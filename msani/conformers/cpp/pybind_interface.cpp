#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "stochastic_sampling.h"
#include "common_types.h"
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <DistGeom/DistGeomUtils.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <GraphMol/ForceFieldHelpers/FFConvenience.h>
#include <memory>
#include <stdexcept>

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

// Enhanced function that creates RDKit molecule from conformers with direct transfer
py::object createMoleculeWithConformersDirectly(py::object mol_obj, const std::vector<ConformerResult>& products) {
    std::unique_ptr<RWMol> mol(extractMolDirectly(mol_obj));
    if (!mol) {
        return py::none();
    }
    
    if (products.empty()) {
        // py::print("No conformers generated from stochastic sampling");
        mol_obj.attr("SetBoolProp")("Failed_sampling", true);
        return mol_obj;  // Return original molecule if no conformers generated
    }
    
    // Check if conformer atom count matches molecule atom count
    if (!products.empty() && products[0].conformer.getNumAtoms() != mol->getNumAtoms()) {
        py::print("Atom count mismatch! Molecule: " + std::to_string(mol->getNumAtoms()) + 
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
    py::object rdgeom = py::module::import("rdkit.Geometry");
    
    for (size_t i = 0; i < products.size(); ++i) {
        const auto& product = products[i];
        const RDKit::Conformer& conf = product.conformer;
        
        // Create Python conformer object
        py::object py_conformer = rdkit_chem.attr("Conformer")(conf.getNumAtoms());
        
        // Set conformer properties
        py_conformer.attr("SetId")(static_cast<int>(i));
        py_conformer.attr("Set3D")(true);
        py_conformer.attr("SetProp")("Energy", std::to_string(product.energy));
        
        // Copy 3D coordinates from C++ conformer to Python conformer
        for (unsigned int atom_idx = 0; atom_idx < conf.getNumAtoms(); ++atom_idx) {
            const RDGeom::Point3D& pos = conf.getAtomPos(atom_idx);
            
            // Set position using Python API
            py::object point3d = rdgeom.attr("Point3D")(pos.x, pos.y, pos.z);
            py_conformer.attr("SetAtomPosition")(atom_idx, point3d);
        }
        
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
        // Extract dihedral atoms (4 integers)
        py::list dihedral_list = value[1].cast<py::list>();
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

// Main discrete sampling wrapper
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
                                             double clash_threshold,
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
            int atom1 = bond_tuple[0].cast<int>();
            int atom2 = bond_tuple[1].cast<int>();
            hetero_H_bonds.emplace_back(atom1, atom2);
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
                                                         clash_threshold,
                                                         verbose,
                                                         mmff_variant,
                                                         eps,
                                                         randomSeed);

        // Convert results back to Python molecule
        return createMoleculeWithConformersDirectly(mol_obj, products);

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
                                               double clash_threshold,
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
            int atom1 = bond_tuple[0].cast<int>();
            int atom2 = bond_tuple[1].cast<int>();
            hetero_H_bonds.emplace_back(atom1, atom2);
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
                               clash_threshold,
                               hetero_H_bonds,
                               verbose,
                               mmff_variant,
                               eps,
                               random_method,
                               randomSeed);
        
        // Convert results back to Python molecule
        return createMoleculeWithConformersDirectly(mol_obj, products);
        
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("Error in continuous sampling: ") + e.what());
    }
}

// PyBind11 module definition
PYBIND11_MODULE(stochastic_sampling_combined, m) {
    m.doc() = "Combined Stochastic Sampling Module - Discrete and Continuous methods";
    
    // Discrete sampling function
    m.def("stochastic_sampling_discrete", &stochasticSamplingDiscreteWrapper,
        "Discrete stochastic sampling using predefined angle values",
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
        py::arg("clash_threshold") = 1.6,
        py::arg("verbose") = false,
        py::arg("mmff_variant") = "MMFF94s",
        py::arg("eps") = 1.0,
        py::arg("randomSeed") = 42);
    
    // Continuous sampling function
    m.def("stochastic_sampling_continuous", &stochasticSamplingContinuousWrapper,
        "Continuous stochastic sampling using random angles around torsional peaks",
        py::arg("mol"),
        py::arg("match_torlib"),
        py::arg("tolerance_level"),
        py::arg("numConfs"),
        py::arg("window") = 25.0,
        py::arg("max_attempts") = 50000,
        py::arg("timeout_conf"),
        py::arg("rmsd") = 0.5,
        py::arg("clash_threshold") = 1.6,
        py::arg("hetero_H_bonds"),
        py::arg("verbose") = false,
        py::arg("mmff_variant") = "MMFF94s",
        py::arg("eps") = 1.0,
        py::arg("random_method") = "uniform",
        py::arg("randomSeed") = 42);
    
    // Version information
    m.attr("__version__") = "0.4.0";
    m.attr("__author__") = "Phong Lam, Uppsala University (2025)";
}