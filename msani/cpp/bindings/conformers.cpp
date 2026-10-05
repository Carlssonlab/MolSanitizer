#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "stochastic_sampling.h"
#include "common_types.h"
#include "initial_embedder.h"
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <DistGeom/DistGeomUtils.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/ForceFieldHelpers/FFConvenience.h>
#include <ForceField/ForceField.h>
#include <ForceField/MMFF/AngleConstraint.h>
#include <ForceField/MMFF/TorsionConstraint.h>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <algorithm>
#include <cctype>
#include <array>
#include <cmath>
#include <set>
#include <optional>
#include <limits>
#include <sstream>

namespace py = pybind11;
using namespace RDKit;
using namespace StochasticSampling;

// Helper function to extract molecule directly from Python RDKit object (no modifications)
RWMol* extractMolDirectly(py::object mol_obj) {
    if (mol_obj.is_none()) {
        return nullptr;
    }
    
    try {
        py::object binary_data = mol_obj.attr("ToBinary")();
        std::string binary_str = binary_data.cast<std::string>();
        
        std::unique_ptr<ROMol> mol(new ROMol(binary_str));
        if (!mol) {
            throw std::runtime_error("Could not deserialize molecule from binary data");
        }
        
        // Create RWMol copy - pickled molecule should have all data intact
        return new RWMol(*mol);
    } catch (...) {
        throw std::runtime_error("Could not extract molecule from Python object using ToBinary method");
    }
}

// Copy a native RDKit conformer's coordinates into a NumPy array so the
// Python RDKit binding can set every atom position in a single call.
py::array_t<double> conformerPositionsToNumpy(const RDKit::Conformer& conf) {
    py::array_t<double> positions({
        static_cast<py::ssize_t>(conf.getNumAtoms()),
        static_cast<py::ssize_t>(3),
    });
    auto coordinates = positions.mutable_unchecked<2>();

    for (unsigned int atom_idx = 0; atom_idx < conf.getNumAtoms(); ++atom_idx) {
        const RDGeom::Point3D& position = conf.getAtomPos(atom_idx);
        coordinates(atom_idx, 0) = position.x;
        coordinates(atom_idx, 1) = position.y;
        coordinates(atom_idx, 2) = position.z;
    }

    return positions;
}

// Enhanced function that creates RDKit molecule from conformers with direct transfer
py::object createMoleculeWithConformersDirectly(py::object mol_obj, const std::vector<ConformerResult>& products) {
    if (mol_obj.is_none()) {
        return py::none();
    }
    
    if (products.empty()) {
        mol_obj.attr("SetBoolProp")("Failed_sampling", true);
        return mol_obj;  // Return original molecule if no conformers generated
    }
    
    // Check if conformer atom count matches molecule atom count without
    // serializing and deserializing the original Python molecule a second time.
    const auto molecule_num_atoms = mol_obj.attr("GetNumAtoms")().cast<unsigned int>();
    if (products[0].conformer.getNumAtoms() != molecule_num_atoms) {
        py::print("Atom count mismatch! Molecule: " + std::to_string(molecule_num_atoms) +
                  ", Conformer: " + std::to_string(products[0].conformer.getNumAtoms()));
        py::print("   This suggests the stochastic sampling process modified the molecular structure.");
        py::print("   Returning original molecule to preserve structure.");
        return mol_obj;
    }
    
    // Clear existing conformers from the molecule
    mol_obj.attr("RemoveAllConformers")();
    mol_obj.attr("SetBoolProp")("Failed_sampling", false);
    // Add conformers directly to the Python molecule object
    py::object rdkit_chem = py::module::import("rdkit.Chem");
    
    for (size_t i = 0; i < products.size(); ++i) {
        const auto& product = products[i];
        const RDKit::Conformer& conf = product.conformer;
        
        // Create Python conformer object
        py::object py_conformer = rdkit_chem.attr("Conformer")(conf.getNumAtoms());
        
        // Set conformer properties
        py_conformer.attr("SetId")(static_cast<int>(i));
        py_conformer.attr("Set3D")(true);
        py_conformer.attr("SetProp")("Energy", std::to_string(product.energy));
        
        // Copy all 3D coordinates in one Python API call.
        py_conformer.attr("SetPositions")(conformerPositionsToNumpy(conf));
        
        // Add conformer directly to the molecule
        mol_obj.attr("AddConformer")(py_conformer, false);  // false = don't assign ID automatically
    }
    
    return mol_obj;
}

// Helper function to convert Python angle map to C++ format for discrete sampling
AngleMap convertDiscreteAngleMap(const py::dict& py_angle_map) {
    AngleMap angle_map;
    
    for (auto item : py_angle_map) {
        int key = item.first.cast<int>();
        py::tuple value = item.second.cast<py::tuple>();
        
        AngleMapEntry entry;
        if (value.size() != 3) throw std::invalid_argument("Angle entries require pattern, four atoms, and angles");
        // Extract dihedral atoms (4 integers)
        py::list dihedral_list = value[1].cast<py::list>();
        if (dihedral_list.size() != 4) throw std::invalid_argument("A torsion requires exactly four atom indices");
        for (int i = 0; i < 4; ++i) {
            entry.dihedral_atoms.push_back(dihedral_list[i].cast<int>());
        }
        
        // Extract possible angles
        py::list angles_list = value[2].cast<py::list>();
        for (auto angle : angles_list) {
            entry.possible_angles.push_back(angle.cast<double>());
        }
        
        angle_map[key] = entry;
    }
    
    return angle_map;
}

// Helper function to convert Python score map to C++ format
ScoreMap convertScoreMap(const py::dict& py_score_map) {
    ScoreMap score_map;
    
    for (auto item : py_score_map) {
        int key = item.first.cast<int>();
        py::list scores_list = item.second.cast<py::list>();
        std::vector<double> scores;
        for (auto score : scores_list) {
            scores.push_back(score.cast<double>());
        }
        score_map[key] = scores;
    }
    
    return score_map;
}

// Helper function to convert Python torsion library to C++ format for continuous sampling
ContinuousTorsionMap convertContinuousTorsionLibrary(const py::list& match_torlib, int tolerance_level) {
    ContinuousTorsionMap torsion_map;
    
    for (size_t i = 0; i < match_torlib.size(); ++i) {
        py::list rule = match_torlib[i];  // Each rule is a list, not tuple
        if (rule.size() != 3) throw std::invalid_argument("Torsion rules require pattern, four atoms, and peaks");
        py::tuple dihedral_atoms_py = rule[1];  // Dihedral atoms tuple
        py::list peaks_py = rule[2];            // Peaks list
        
        // Extract dihedral atoms
        std::vector<int> dihedral_atoms;
        for (size_t j = 0; j < dihedral_atoms_py.size(); ++j) {
            dihedral_atoms.push_back(dihedral_atoms_py[j].cast<int>());
        }
        
        // Extract peaks
        std::vector<TorsionPeak> peaks;
        for (size_t j = 0; j < peaks_py.size(); ++j) {
            py::tuple peak_tuple = peaks_py[j];
            
            if (peak_tuple.size() != 4) throw std::invalid_argument("Peaks require center, two tolerances, and weight");
            double center = peak_tuple[0].cast<double>();
            double tolerance_1 = peak_tuple[1].cast<double>();
            double tolerance_2 = peak_tuple[2].cast<double>();
            double weight = peak_tuple[3].cast<double>();

            std::vector<double> tolerance;
            tolerance.push_back(tolerance_1);
            tolerance.push_back(tolerance_2);
            
            peaks.emplace_back(center, tolerance, weight);
        }
        
        // Create bond info
        ContinuousTorsionBond bond_info;
        bond_info.dihedral_atoms = dihedral_atoms;
        bond_info.peaks = peaks;
        
        torsion_map[static_cast<int>(i)] = bond_info;
    }
    
    return torsion_map;
}

struct TorsionConstraintData {
    std::vector<std::array<unsigned int, 6>> conjugated_substituted_nitrogen_5aro;
    std::vector<std::array<unsigned int, 7>> conjugated_substituted_nitrogen_6aro;
    std::vector<std::vector<unsigned int>> barbiturate_matches;
    std::vector<std::vector<unsigned int>> hydantoin_matches;
    std::vector<std::array<unsigned int, 4>> substituted_N_barbi_hydan_like;
    std::vector<std::vector<unsigned int>> planar_rings;
    std::vector<std::array<unsigned int, 4>> alkyne;

    bool empty() const {
        return conjugated_substituted_nitrogen_5aro.empty() &&
               conjugated_substituted_nitrogen_6aro.empty() &&
               barbiturate_matches.empty() &&
               hydantoin_matches.empty() &&
               substituted_N_barbi_hydan_like.empty() &&
               planar_rings.empty() &&
               alkyne.empty();
    }
};

std::vector<unsigned int> convertIndexSequence(const py::handle &seq_handle) {
    py::sequence seq = seq_handle.cast<py::sequence>();
    const Py_ssize_t sequence_size = py::len(seq);
    std::vector<unsigned int> result;
    result.reserve(static_cast<std::size_t>(sequence_size));
    for (Py_ssize_t i = 0; i < sequence_size; ++i) {
        result.push_back(seq[i].cast<unsigned int>());
    }
    return result;
}

template <size_t N>
std::optional<std::array<unsigned int, N>> convertFixedSizeArray(const py::handle &seq_handle) {
    py::sequence seq = seq_handle.cast<py::sequence>();
    if (py::len(seq) != static_cast<Py_ssize_t>(N)) {
        return std::nullopt;
    }
    std::array<unsigned int, N> arr{};
    for (size_t i = 0; i < N; ++i) {
        arr[i] = seq[i].cast<unsigned int>();
    }
    return arr;
}

TorsionConstraintData parseTorsionConstraintData(const py::object& constraints_obj) {
    TorsionConstraintData data;
    if (constraints_obj.is_none()) {
        return data;
    }

    py::dict constraint_dict;
    bool is_dict = py::isinstance<py::dict>(constraints_obj);
    if (is_dict) {
        constraint_dict = constraints_obj.cast<py::dict>();
    }

    auto fetch = [&](const char *name) -> py::object {
        if (is_dict) {
            py::str key(name);
            if (constraint_dict.contains(key)) {
                return constraint_dict[key];
            }
        } else if (py::hasattr(constraints_obj, name)) {
            return constraints_obj.attr(name);
        }
        return py::none();
    };

    auto conj5_obj = fetch("conjugated_substituted_nitrogen_5aro");
    if (!conj5_obj.is_none() && py::len(conj5_obj)) {
        py::list conj5_list = conj5_obj.cast<py::list>();
        for (auto entry : conj5_list) {
            auto arr = convertFixedSizeArray<6>(entry);
            if (arr) {
                data.conjugated_substituted_nitrogen_5aro.push_back(*arr);
            }
        }
    }

    auto conj6_obj = fetch("conjugated_substituted_nitrogen_6aro");
    if (!conj6_obj.is_none() && py::len(conj6_obj)) {
        py::list conj6_list = conj6_obj.cast<py::list>();
        for (auto entry : conj6_list) {
            auto arr = convertFixedSizeArray<7>(entry);
            if (arr) {
                data.conjugated_substituted_nitrogen_6aro.push_back(*arr);
            }
        }
    }

    auto barbiturate_obj = fetch("barbiturate_matches");
    if (!barbiturate_obj.is_none() && py::len(barbiturate_obj)) {
        py::list matches = barbiturate_obj.cast<py::list>();
        for (auto entry : matches) {
            auto vec = convertIndexSequence(entry);
            if (vec.size() >= 4) {
                data.barbiturate_matches.push_back(std::move(vec));
            }
        }
    }

    auto hydantoin_obj = fetch("hydantoin_matches");
    if (!hydantoin_obj.is_none() && py::len(hydantoin_obj)) {
        py::list matches = hydantoin_obj.cast<py::list>();
        for (auto entry : matches) {
            auto vec = convertIndexSequence(entry);
            if (vec.size() >= 4) {
                data.hydantoin_matches.push_back(std::move(vec));
            }
        }
    }

    auto substituted_obj = fetch("substituted_N_barbi_hydan_like");
    if (!substituted_obj.is_none() && py::len(substituted_obj)) {
        py::list entries = substituted_obj.cast<py::list>();
        for (auto entry : entries) {
            auto arr = convertFixedSizeArray<4>(entry);
            if (arr) {
                data.substituted_N_barbi_hydan_like.push_back(*arr);
            }
        }
    }

    auto planar_obj = fetch("planar_rings");
    if (!planar_obj.is_none() && py::len(planar_obj)) {
        py::list rings = planar_obj.cast<py::list>();
        for (auto entry : rings) {
            auto vec = convertIndexSequence(entry);
            if (vec.size() >= 4) {
                data.planar_rings.push_back(std::move(vec));
            }
        }
    }

    auto alkyne_obj = fetch("alkyne");
    if (!alkyne_obj.is_none() && py::len(alkyne_obj)) {
        py::list entries = alkyne_obj.cast<py::list>();
        for (auto entry : entries) {
            auto arr = convertFixedSizeArray<4>(entry);
            if (arr) {
                data.alkyne.push_back(*arr);
            }
        }
    }

    return data;
}

void addTorsionConstraint(ForceFields::ForceField &ff,
                          unsigned int a,
                          unsigned int b,
                          unsigned int c,
                          unsigned int d,
                          double minDeg,
                          double maxDeg,
                          double forceConstant) {
    auto *constraint = new ForceFields::MMFF::TorsionConstraintContrib(
        &ff, a, b, c, d, false, minDeg, maxDeg, forceConstant);
    ff.contribs().push_back(ForceFields::ContribPtr(constraint));
}

void addAngleConstraint(ForceFields::ForceField &ff,
                        unsigned int a,
                        unsigned int b,
                        unsigned int c,
                        double minDeg,
                        double maxDeg,
                        double forceConstant) {
    auto *constraint = new ForceFields::MMFF::AngleConstraintContrib(
        &ff, a, b, c, false, minDeg, maxDeg, forceConstant);
    ff.contribs().push_back(ForceFields::ContribPtr(constraint));
}

void applyTorsionConstraints(ForceFields::ForceField &ff, const TorsionConstraintData &data) {
    constexpr double tightMin = -2.0;
    constexpr double tightMax = 2.0;
    constexpr double transMin = 178.0;
    constexpr double transMax = 182.0;
    constexpr double alkyneMin = 179.5;
    constexpr double alkyneMax = 180.5;
    constexpr double alkyneForce = 5.0;
    constexpr double alkyneBendMin = 170;
    // RDKit bond angles cannot exceed 180 degrees, so 180.0 is the
    // realizable upper bound of the requested 179.5-180.5 degree interval.
    constexpr double alkyneBendMax = 180.0;
    constexpr double alkyneBendForce = 5.0;
    constexpr double defaultForce = 1.0;

    for (const auto &entry : data.conjugated_substituted_nitrogen_5aro) {
        unsigned int a = entry[0];
        unsigned int b = entry[1];
        unsigned int c = entry[2];
        unsigned int d = entry[3];
        unsigned int e = entry[4];
        unsigned int f = entry[5];
        addTorsionConstraint(ff, a, b, c, d, transMin, transMax, defaultForce);
        addTorsionConstraint(ff, a, b, f, e, transMin, transMax, defaultForce);
        addTorsionConstraint(ff, b, c, d, e, tightMin, tightMax, defaultForce);
        addTorsionConstraint(ff, d, e, f, b, tightMin, tightMax, defaultForce);
    }

    for (const auto &entry : data.conjugated_substituted_nitrogen_6aro) {
        unsigned int a = entry[0];
        unsigned int b = entry[1];
        unsigned int c = entry[2];
        unsigned int d = entry[3];
        unsigned int e = entry[4];
        unsigned int f = entry[5];
        unsigned int g = entry[6];
        addTorsionConstraint(ff, a, b, c, d, transMin, transMax, defaultForce);
        addTorsionConstraint(ff, a, b, g, f, transMin, transMax, defaultForce);
        addTorsionConstraint(ff, b, c, d, e, tightMin, tightMax, defaultForce);
        addTorsionConstraint(ff, e, f, g, b, tightMin, tightMax, defaultForce);
    }

    for (const auto &match : data.barbiturate_matches) {
        const size_t n = match.size();
        if (n < 4) {
            continue;
        }
        for (size_t i = 0; i < n; ++i) {
            unsigned int a = match[i % n];
            unsigned int b = match[(i + 1) % n];
            unsigned int c = match[(i + 2) % n];
            unsigned int d = match[(i + 3) % n];
            addTorsionConstraint(ff, a, b, c, d, tightMin, tightMax, defaultForce);
        }
    }

    for (const auto &match : data.hydantoin_matches) {
        const size_t n = match.size();
        if (n < 4) {
            continue;
        }
        for (size_t i = 0; i < n; ++i) {
            unsigned int a = match[i % n];
            unsigned int b = match[(i + 1) % n];
            unsigned int c = match[(i + 2) % n];
            unsigned int d = match[(i + 3) % n];
            addTorsionConstraint(ff, a, b, c, d, tightMin, tightMax, defaultForce);
        }
    }

    for (const auto &entry : data.substituted_N_barbi_hydan_like) {
        addTorsionConstraint(ff, entry[0], entry[1], entry[2], entry[3], transMin, transMax, defaultForce);
    }

    for (const auto &ring : data.planar_rings) {
        const size_t n = ring.size();
        if (n < 4) {
            continue;
        }
        for (size_t i = 0; i + 1 < n; ++i) {
            unsigned int a = ring[i];
            unsigned int b = ring[(i + 1) % n];
            unsigned int c = ring[(i + 2) % n];
            unsigned int d = ring[(i + 3) % n];
            addTorsionConstraint(ff, a, b, c, d, tightMin, tightMax, defaultForce);
        }
    }

    for (const auto &entry : data.alkyne) {
        addTorsionConstraint(ff, entry[0], entry[1], entry[2], entry[3],
                             alkyneMin, alkyneMax, alkyneForce);
        addAngleConstraint(ff, entry[0], entry[1], entry[2],
                           alkyneBendMin, alkyneBendMax, alkyneBendForce);
        addAngleConstraint(ff, entry[1], entry[2], entry[3],
                           alkyneBendMin, alkyneBendMax, alkyneBendForce);
    }
}

// Main discrete sampling wrapper
namespace {
void requireFinite(double value, const char* name, bool positive = false) {
    if (!std::isfinite(value) || (positive ? value <= 0.0 : value < 0.0)) {
        throw std::invalid_argument(std::string(name) + " must be finite and " +
                                    (positive ? "positive" : "nonnegative"));
    }
}

void validateAtom(int index, const ROMol& mol) {
    if (index < 0 || static_cast<unsigned int>(index) >= mol.getNumAtoms()) {
        throw std::invalid_argument("Torsion atom index is outside the molecule");
    }
}

void validateTorsionAtoms(const std::vector<int>& atoms, const ROMol& mol) {
    if (atoms.size() != 4 || std::set<int>(atoms.begin(), atoms.end()).size() != 4) {
        throw std::invalid_argument("A torsion requires four distinct atom indices");
    }
    for (int atom : atoms) validateAtom(atom, mol);
    for (std::size_t i = 1; i < atoms.size(); ++i) {
        if (!mol.getBondBetweenAtoms(atoms[i - 1], atoms[i])) {
            throw std::invalid_argument("Torsion atoms must form a bonded path");
        }
    }
}

void validateWeights(const std::vector<double>& weights, std::size_t expected) {
    if (weights.size() != expected || weights.empty()) {
        throw std::invalid_argument("Sampling weights must match the nonempty choices");
    }
    double sum = 0.0;
    for (double weight : weights) { requireFinite(weight, "Weight"); sum += weight; }
    requireFinite(sum, "Weight sum", true);
}

void validateSamplingParameters(const ROMol& mol, std::size_t torsions,
                                const HeteroBonds& bonds, int count, int attempts,
                                int timeout, double window, double rmsd,
                                double clash_scale, double eps) {
    if (count <= 0 || attempts <= 0 || timeout < 0) {
        throw std::invalid_argument("Conformer and attempt counts must be positive; timeout must be nonnegative");
    }
    requireFinite(window, "Energy window");
    requireFinite(rmsd, "RMSD");
    requireFinite(clash_scale, "Clash scale", true);
    requireFinite(eps, "Dielectric constant", true);
    if (bonds.size() > torsions) {
        throw std::invalid_argument("More hydroxyl bonds than torsion entries");
    }
    for (const auto& bond : bonds) {
        validateAtom(bond.first, mol); validateAtom(bond.second, mol);
        if (!mol.getBondBetweenAtoms(bond.first, bond.second)) {
            throw std::invalid_argument("Hydroxyl atom indices must describe a bond");
        }
    }
}
}  // namespace

// Main discrete sampling wrapper (fixed argument ordering)
py::object stochasticSamplingDiscreteWrapper(py::object mol_obj,
                                             const py::dict& py_angle_map,
                                             const py::dict& py_score_map,
                                             long long possible_numConfs,
                                             const py::list& py_importance_order,
                                             double window,
                                             int max_attempts,
                                             const py::list& py_hetero_H_bonds,
                                             int timeout_conf,
                                             double rmsd,
                                             int numConfs,
                                             double clash_scale,
                                             bool verbose = false,
                                             const std::string& mmff_variant = "MMFF94s",
                                             double eps = 1.0,
                                             int randomSeed = 42) {
    try {
        // Extract molecule
        std::unique_ptr<RWMol> mol(extractMolDirectly(mol_obj));
        if (!mol) {
            throw std::runtime_error("Failed to extract molecule");
        }

        // Convert Python data structures to C++
        AngleMap angle_map = convertDiscreteAngleMap(py_angle_map);
        ScoreMap score_map = convertScoreMap(py_score_map);

        // Convert importance order
        ImportanceOrder importance_order;
        for (auto item : py_importance_order) {
            importance_order.push_back(item.cast<double>());
        }

        // Convert hetero H bonds
        HeteroBonds hetero_H_bonds;
        for (auto item : py_hetero_H_bonds) {
            py::tuple bond_tuple = item.cast<py::tuple>();
            if (bond_tuple.size() != 2) throw std::invalid_argument("Hydroxyl bonds require two atom indices");
            int atom1 = bond_tuple[0].cast<int>();
            int atom2 = bond_tuple[1].cast<int>();
            hetero_H_bonds.emplace_back(atom1, atom2);
        }

        validateSamplingParameters(*mol, angle_map.size(), hetero_H_bonds,
                                   numConfs, max_attempts, timeout_conf, window, rmsd, clash_scale, eps);
        if (possible_numConfs < 0) throw std::invalid_argument("Possible conformer count must be nonnegative");
        if (!angle_map.empty()) validateWeights(importance_order, angle_map.size());
        for (const auto& [key, entry] : angle_map) {
            validateTorsionAtoms(entry.dihedral_atoms, *mol);
            if (entry.possible_angles.empty()) throw std::invalid_argument("Torsion angles cannot be empty");
            for (double angle : entry.possible_angles) {
                if (!std::isfinite(angle)) throw std::invalid_argument("Torsion angles must be finite");
            }
            const auto weights = score_map.find(key);
            if (weights != score_map.end()) validateWeights(weights->second, entry.possible_angles.size());
        }
        // Call the discrete sampling function
        ProductList products = stochasticSamplingDiscrete(*mol,
                                                         angle_map,
                                                         score_map,
                                                         possible_numConfs,
                                                         importance_order,
                                                         window,
                                                         max_attempts,
                                                         hetero_H_bonds,
                                                         timeout_conf,
                                                         rmsd,
                                                         numConfs,
                                                         clash_scale,
                                                         verbose,
                                                         mmff_variant,
                                                         eps,
                                                         randomSeed);

        // Convert results back to Python molecule
        return createMoleculeWithConformersDirectly(mol_obj, products);

    } catch (const std::invalid_argument&) {
        throw;
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("Error in discrete sampling: ") + e.what());
    }
}

// Main continuous sampling wrapper
py::object stochasticSamplingContinuousWrapper(py::object mol_obj,
                                               const py::list& match_torlib,
                                               int tolerance_level,
                                               int numConfs,
                                               double window,
                                               int max_attempts,
                                               int timeout_conf,
                                               double rmsd,
                                               double clash_scale,
                                               const py::list& py_hetero_H_bonds,
                                               bool verbose = false,
                                               const std::string& mmff_variant = "MMFF94s",
                                               double eps = 1.0,
                                               const std::string& random_method = "uniform",
                                               int randomSeed = 42) {
    
    try {
        // Extract molecule
        std::unique_ptr<RWMol> mol(extractMolDirectly(mol_obj));
        if (!mol) {
            throw std::runtime_error("Failed to extract molecule");
        }
        // Convert torsion library
        ContinuousTorsionMap torsion_library = convertContinuousTorsionLibrary(match_torlib, tolerance_level);
        // Convert hetero H bonds
        HeteroBonds hetero_H_bonds;
        for (auto item : py_hetero_H_bonds) {
            py::tuple bond_tuple = item.cast<py::tuple>();
            if (bond_tuple.size() != 2) throw std::invalid_argument("Hydroxyl bonds require two atom indices");
            int atom1 = bond_tuple[0].cast<int>();
            int atom2 = bond_tuple[1].cast<int>();
            hetero_H_bonds.emplace_back(atom1, atom2);
        }
        validateSamplingParameters(*mol, torsion_library.size(), hetero_H_bonds,
                                   numConfs, max_attempts, timeout_conf, window, rmsd, clash_scale, eps);
        if (tolerance_level != 1 && tolerance_level != 2) {
            throw std::invalid_argument("Tolerance level must be 1 or 2");
        }
        for (const auto& [key, bond] : torsion_library) {
            validateTorsionAtoms(bond.dihedral_atoms, *mol);
            std::vector<double> weights;
            for (const auto& peak : bond.peaks) {
                if (!std::isfinite(peak.center)) throw std::invalid_argument("Peak centers must be finite");
                for (double tolerance : peak.tolerance) requireFinite(tolerance, "Peak tolerance");
                weights.push_back(peak.weight);
            }
            validateWeights(weights, bond.peaks.size());
        }
        // Call the continuous sampling function
    ProductList products = stochasticSamplingContinuous(*mol,
                               torsion_library,
                               tolerance_level,
                               numConfs,
                               window,
                               max_attempts,
                               timeout_conf,
                               rmsd,
                               clash_scale,
                               hetero_H_bonds,
                               verbose,
                               mmff_variant,
                               eps,
                               random_method,
                               randomSeed);
        
        // Convert results back to Python molecule
        return createMoleculeWithConformersDirectly(mol_obj, products);
        
    } catch (const std::invalid_argument&) {
        throw;
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("Error in continuous sampling: ") + e.what());
    }
}

py::object embedMultipleConfsWrapper(py::object mol_obj,
                                     unsigned int numConfs,
                                     py::object params_obj = py::none(),
                                     py::object constraints_obj = py::none()) {
    try {
        std::unique_ptr<RWMol> mol(extractMolDirectly(mol_obj));
        if (!mol) {
            throw std::runtime_error("Could not extract molecule");
        }

        static const std::unordered_map<std::string, const DGeomHelpers::EmbedParameters *> presetMap = {
            {"kdg", &InitialEmbedder::KDG},
            {"etdg", &InitialEmbedder::ETDG},
            {"etdgv2", &InitialEmbedder::ETDGv2},
            {"etkdg", &InitialEmbedder::ETKDG},
            {"etkdgv2", &InitialEmbedder::ETKDGv2},
            {"etkdgv3", &InitialEmbedder::ETKDGv3},
            {"sretkdgv3", &InitialEmbedder::srETKDGv3}
        };

        DGeomHelpers::EmbedParameters params = InitialEmbedder::ETKDGv3;
        if (!params_obj.is_none()) {
            if (py::isinstance<py::str>(params_obj)) {
                std::string preset = params_obj.cast<std::string>();
                std::string lowered;
                lowered.resize(preset.size());
                std::transform(preset.begin(), preset.end(), lowered.begin(),
                               [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
                auto it = presetMap.find(lowered);
                if (it != presetMap.end()) {
                    params = *(it->second);
                } else {
                    throw std::invalid_argument("Unknown embedding preset: " + preset);
                }
            } else {
                auto copyUnsigned = [&](const char *name, unsigned int &target) {
                    if (py::hasattr(params_obj, name)) {
                        target = params_obj.attr(name).cast<unsigned int>();
                    }
                };
                auto copyInt = [&](const char *name, int &target) {
                    if (py::hasattr(params_obj, name)) {
                        target = params_obj.attr(name).cast<int>();
                    }
                };
                auto copyDouble = [&](const char *name, double &target) {
                    if (py::hasattr(params_obj, name)) {
                        target = params_obj.attr(name).cast<double>();
                    }
                };
                auto copyBool = [&](const char *name, bool &target) {
                    if (py::hasattr(params_obj, name)) {
                        target = params_obj.attr(name).cast<bool>();
                    }
                };

                copyUnsigned("maxIterations", params.maxIterations);
                copyInt("numThreads", params.numThreads);
                copyInt("randomSeed", params.randomSeed);
                copyBool("clearConfs", params.clearConfs);
                copyBool("useRandomCoords", params.useRandomCoords);
                copyDouble("boxSizeMult", params.boxSizeMult);
                copyBool("randNegEig", params.randNegEig);
                copyUnsigned("numZeroFail", params.numZeroFail);
                copyDouble("optimizerForceTol", params.optimizerForceTol);
                copyBool("ignoreSmoothingFailures", params.ignoreSmoothingFailures);
                copyBool("enforceChirality", params.enforceChirality);
                copyBool("useExpTorsionAnglePrefs", params.useExpTorsionAnglePrefs);
                copyBool("useBasicKnowledge", params.useBasicKnowledge);
                copyBool("verbose", params.verbose);
                copyDouble("basinThresh", params.basinThresh);
                copyDouble("pruneRmsThresh", params.pruneRmsThresh);
                copyBool("onlyHeavyAtomsForRMS", params.onlyHeavyAtomsForRMS);
                copyUnsigned("ETversion", params.ETversion);
                copyBool("embedFragmentsSeparately", params.embedFragmentsSeparately);
                copyBool("useSmallRingTorsions", params.useSmallRingTorsions);
                copyBool("useMacrocycleTorsions", params.useMacrocycleTorsions);
                copyBool("useMacrocycle14config", params.useMacrocycle14config);
                copyUnsigned("timeout", params.timeout);
                copyBool("trackFailures", params.trackFailures);
                copyBool("useSymmetryForPruning", params.useSymmetryForPruning);
                copyBool("symmetrizeConjugatedTerminalGroupsForPruning",
                         params.symmetrizeConjugatedTerminalGroupsForPruning);
                copyBool("enableSequentialRandomSeeds",
                         params.enableSequentialRandomSeeds);
                copyBool("forceTransAmides", params.forceTransAmides);
                copyDouble("boundsMatForceScaling", params.boundsMatForceScaling);
            }
        }

        TorsionConstraintData constraint_data = parseTorsionConstraintData(constraints_obj);

        // Extract MMFF variant from constraints if available
        std::string mmff_variant = "MMFF94s"; // default
        if (!constraints_obj.is_none()) {
            py::dict constraint_dict;
            bool is_dict = py::isinstance<py::dict>(constraints_obj);
            if (is_dict) {
                constraint_dict = constraints_obj.cast<py::dict>();
                py::str key("forcefield");
                if (constraint_dict.contains(key)) {
                    mmff_variant = constraint_dict[key].cast<std::string>();
                }
            } else if (py::hasattr(constraints_obj, "forcefield")) {
                mmff_variant = constraints_obj.attr("forcefield").cast<std::string>();
            }
        }

        INT_VECT res;
        InitialEmbedder::EmbedMultipleConfs(*mol, res, numConfs, params);

        // Create MMFF properties using the same pattern as RDKit's MMFFOptimizeMoleculeConfs
        RDKit::MMFF::MMFFMolProperties mmff_properties(*mol, mmff_variant);
        if (!mmff_properties.isValid()) {
            throw std::runtime_error("Failed to compute valid MMFF properties for molecule");
        }
        
        // Turn off electrostatic term
        mmff_properties.setMMFFEleTerm(false);

        // Pre-construct the force field once instead of creating it for each conformer
        std::unique_ptr<ForceFields::ForceField> base_ff;
        if (!res.empty() && res[0] >= 0) {
            base_ff.reset(RDKit::MMFF::constructForceField(*mol, &mmff_properties, 100.0, -1, true));
            if (base_ff) {
                // Apply torsion constraints to the base force field
                applyTorsionConstraints(*base_ff, constraint_data);
            }
        }

        mol_obj.attr("RemoveAllConformers")();
        py::object rdkit_chem = py::module::import("rdkit.Chem");

        for (size_t i = 0; i < res.size(); ++i) {
            int confId = res[i];
            if (confId < 0 || !base_ff) {
                continue;
            }

            // Update force field positions to point to current conformer
            Conformer& mutable_conf = mol->getConformer(confId);
            for (unsigned int atom_idx = 0; atom_idx < mol->getNumAtoms(); ++atom_idx) {
                base_ff->positions()[atom_idx] = &mutable_conf.getAtomPos(atom_idx);
            }
            
            // Initialize and minimize using the shared force field
            base_ff->initialize();
            base_ff->minimize();
            double energy = base_ff->calcEnergy();

            py::object py_conformer = rdkit_chem.attr("Conformer")(mutable_conf.getNumAtoms());
            py_conformer.attr("SetId")(confId);
            py_conformer.attr("Set3D")(true);
            py_conformer.attr("SetDoubleProp")("MMFF_Energy", energy);

            py_conformer.attr("SetPositions")(conformerPositionsToNumpy(mutable_conf));

            mol_obj.attr("AddConformer")(py_conformer, false);
        }


        return mol_obj;
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("Error in embed_multiple_confs: ") + e.what());
    }
}

// PyBind11 module definition
PYBIND11_MODULE(msani_confgen_cpp, m) {
    m.doc() = "The MolSanitizer C++ Conformer Generation Module\n\n"
              "This module provides high-performance conformer generation using stochastic sampling\n"
              "and initial embedding with torsion constraints for molecular structures.";
    
    // Version information
    m.attr("__author__") = "Phong Lam, Uppsala University (2025)";

    // Discrete sampling function
    m.def("stochastic_sampling_discrete", &stochasticSamplingDiscreteWrapper,
        R"pbdoc(
        Discrete stochastic sampling using predefined angle values.
        
        Generates conformers by sampling from discrete torsion angle sets with
        importance-weighted selection and MMFF energy evaluation.
        
        Args:
            mol (rdkit.Chem.Mol): Input molecule object
            angle_map (dict): Dictionary mapping bond indices to angle data
                Format: {bond_id: (pattern, [atom1, atom2, atom3, atom4], [angle1, angle2, ...])}
            score_map (dict): Dictionary mapping bond indices to score arrays
                Format: {bond_id: [score1, score2, ...]}
            possible_numConfs (int): Maximum number of conformers to attempt
            importance_order (list): List of importance weights for bond selection
            window (float, optional): Energy window for conformer acceptance. Defaults to 25.0.
            max_attempts (int, optional): Maximum sampling attempts. Defaults to 50000.
            hetero_H_bonds (list): List of heterogeneous hydrogen bond pairs
                Format: [(atom1, atom2), ...]
            timeout_conf (int): Timeout per conformer in seconds
            rmsd (float, optional): RMSD threshold for conformer pruning. Defaults to 0.5.
            numConfs (int, optional): Target number of output conformers. Defaults to 600.
            clash_scale (float, optional): Scale applied to the sum of atomic van der Waals radii for clash detection. Defaults to 0.7.
            verbose (bool, optional): Enable verbose output. Defaults to False.
            mmff_variant (str, optional): MMFF variant to use. Defaults to "MMFF94s".
            eps (float, optional): Numerical epsilon for calculations. Defaults to 1.0.
            randomSeed (int, optional): Random seed for reproducibility. Defaults to 42.
            
        Returns:
            rdkit.Chem.Mol: Molecule with generated conformers attached
            
        Raises:
            RuntimeError: If molecule extraction or sampling fails
        )pbdoc",
        py::arg("mol"),
        py::arg("angle_map"),
        py::arg("score_map"),
        py::arg("possible_numConfs"),
        py::arg("importance_order"),
        py::arg("window") = 25.0,
        py::arg("max_attempts") = 50000,
        py::arg("hetero_H_bonds"),
        py::arg("timeout_conf"),
        py::arg("rmsd") = 0.5,
        py::arg("numConfs") = 600,
        py::arg("clash_scale") = 0.6,
        py::arg("verbose") = false,
        py::arg("mmff_variant") = "MMFF94s",
        py::arg("eps") = 1.0,
        py::arg("randomSeed") = 42);
    
    // Continuous sampling function
    m.def("stochastic_sampling_continuous", &stochasticSamplingContinuousWrapper,
        R"pbdoc(
        Continuous stochastic sampling using random angles around torsional peaks.
        
        Generates conformers by sampling continuous angle distributions around
        known torsional energy minima with Gaussian or uniform sampling.
        
        Args:
            mol (rdkit.Chem.Mol): Input molecule object
            match_torlib (list): Torsion library matching data
                Format: [rule1, rule2, ...] where each rule is 
                [pattern, (atom1, atom2, atom3, atom4), [(peak1_center, tol1, tol2, weight1), ...]]
            tolerance_level (int): Tolerance level for peak sampling (1 or 2, selecting tolerance1 or tolerance2)
            numConfs (int): Target number of output conformers
            window (float, optional): Energy window for conformer acceptance. Defaults to 25.0.
            max_attempts (int, optional): Maximum sampling attempts. Defaults to 50000.
            timeout_conf (int): Timeout per conformer in seconds
            rmsd (float, optional): RMSD threshold for conformer pruning. Defaults to 0.5.
            clash_scale (float, optional): Scale applied to the sum of atomic van der Waals radii for clash detection. Defaults to 0.7.
            hetero_H_bonds (list): List of heterogeneous hydrogen bond pairs
                Format: [(atom1, atom2), ...]
            verbose (bool, optional): Enable verbose output. Defaults to False.
            mmff_variant (str, optional): MMFF variant to use. Defaults to "MMFF94s".
            eps (float, optional): Numerical epsilon for calculations. Defaults to 1.0.
            random_method (str, optional): Random sampling method ("uniform" or "gaussian"). Defaults to "uniform".
            randomSeed (int, optional): Random seed for reproducibility. Defaults to 42.
            
        Returns:
            rdkit.Chem.Mol: Molecule with generated conformers attached
            
        Raises:
            RuntimeError: If molecule extraction or sampling fails
        )pbdoc",
        py::arg("mol"),
        py::arg("match_torlib"),
        py::arg("tolerance_level"),
        py::arg("numConfs"),
        py::arg("window") = 25.0,
        py::arg("max_attempts") = 50000,
        py::arg("timeout_conf"),
        py::arg("rmsd") = 0.5,
        py::arg("clash_scale") = 0.6,
        py::arg("hetero_H_bonds"),
        py::arg("verbose") = false,
        py::arg("mmff_variant") = "MMFF94s",
        py::arg("eps") = 1.0,
        py::arg("random_method") = "uniform",
        py::arg("randomSeed") = 42);
    
    // Initial embedder function
    m.def("embed_multiple_confs", &embedMultipleConfsWrapper,
        R"pbdoc(
        Embed multiple conformers using RDKit with torsion constraints and MMFF minimization.
        
        Generates initial conformers using RDKit's distance geometry embedding,
        applies optional torsion constraints, performs MMFF minimization, and
        calculates final energies.
        
        Args:
            mol (rdkit.Chem.Mol): Input molecule object
            numConfs (int): Number of conformers to generate
            params (str or object, optional): Embedding parameters. Can be:
                - String preset: "kdg", "etdg", "etkdg", "etkdgv2", "etkdgv3", "sretkdgv3"
                - RDKit EmbedParameters object with custom settings
                - None for default ETKDGv3 parameters
            constraints (dict or object, optional): Torsion constraint specifications.
                Can be dictionary or object with any of these attributes:
                - forcefield (str): MMFF variant to use (e.g., "MMFF94", "MMFF94s"). Defaults to "MMFF94s".
                - conjugated_substituted_nitrogen_5aro: List of 6-atom index arrays
                - conjugated_substituted_nitrogen_6aro: List of 7-atom index arrays  
                - barbiturate_matches: List of variable-length atom index arrays
                - hydantoin_matches: List of variable-length atom index arrays
                - substituted_N_barbi_hydan_like: List of 4-atom index arrays
                - planar_rings: List of variable-length atom index arrays (≥4 atoms)
                - alkyne: List of 4-atom index arrays with a linear torsion constraint and
                  179.5–180 degree bending constraints on atoms 0-1-2 and 1-2-3
                
        Returns:
            rdkit.Chem.Mol: Input molecule with conformers added, each having:
                - 3D coordinates from embedding + minimization
                - "MMFF_Energy" property with final energy value
                
        Raises:
            RuntimeError: If molecule extraction, embedding, or minimization fails
            ValueError: If constraint data format is invalid
            
        Example:
            >>> import msani_confgen_cpp as mcc
            >>> from rdkit import Chem
            >>> mol = Chem.MolFromSmiles("CCO")
            >>> mol = Chem.AddHs(mol)
            >>> constraints = {
            ...     "forcefield": "MMFF94s",
            ...     "planar_rings": [[0, 1, 2, 3]]
            ... }
            >>> result_mol = mcc.embed_multiple_confs(mol, 10, "etkdgv3", constraints)
            >>> conf = result_mol.GetConformer(0)
            >>> energy = conf.GetDoubleProp("MMFF_Energy")
        )pbdoc",
        py::arg("mol"),
        py::arg("numConfs"),
        py::arg("params") = py::none(),
        py::arg("constraints") = py::none());
}
