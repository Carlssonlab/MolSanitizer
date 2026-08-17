#include "common_types.h"
#include "accelerated_rmsd.h"
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <GraphMol/ForceFieldHelpers/FFConvenience.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/PeriodicTable.h>
#include <ForceField/ForceField.h>
#include <Geometry/point.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <functional>
#include <unordered_set>

namespace StochasticSampling {

// ========================
// CONFORMER RESULT IMPLEMENTATION
// ========================

RDKit::Conformer ConformerResult::createStandaloneConformer(const RDKit::Conformer& source) {
    RDKit::Conformer result(source.getNumAtoms());
    
    // Copy all atom positions
    for (unsigned int i = 0; i < source.getNumAtoms(); ++i) {
        result.setAtomPos(i, source.getAtomPos(i));
    }
    
    // Copy conformer properties if any
    result.setId(source.getId());
    
    return result;
}

ConformerResult::ConformerResult(const RDKit::Conformer& conf, double e) 
    : conformer(createStandaloneConformer(conf)), energy(e) {}

// ========================
// CONFORMER CACHE IMPLEMENTATION
// ========================

void ConformerCache::initialize(const RDKit::ROMol& reference_mol) {
    // Create H-stripped molecule template for efficient RDKit RMSD comparison
    cache_mol_no_h = std::make_unique<RDKit::RWMol>(reference_mol);
    RDKit::MolOps::removeAllHs(*cache_mol_no_h);
    cache_mol_no_h->clearConformers();
    energies.clear();
    cache_ref_radii.clear();
    cache_ref_centroids.clear();
    
    // Initialize accelerated RMSD calculator
    if (!rmsd_calculator) {
        rmsd_calculator = std::make_unique<AcceleratedRMSD::SameMoleculeRMSDCalculator>();
    }
    rmsd_calculator->initialize(*cache_mol_no_h);
}

ConformerCache::PreparedHeavyConformer ConformerCache::prepareHeavyConformer(
        const RDKit::Conformer& conf,
        const std::vector<int>& heavy_atom_mapping) const {
    if (!cache_mol_no_h) {
        throw std::runtime_error("ConformerCache not initialized");
    }

    PreparedHeavyConformer prepared;
    prepared.conformer = std::make_unique<RDKit::Conformer>(cache_mol_no_h->getNumAtoms());
    for (unsigned int i = 0; i < conf.getNumAtoms(); ++i) {
        const int heavy_idx = heavy_atom_mapping[i];
        if (heavy_idx >= 0) {
            prepared.conformer->setAtomPos(heavy_idx, conf.getAtomPos(i));
        }
    }

    const size_t n_heavy = cache_mol_no_h->getNumAtoms();
    prepared.radii.resize(n_heavy);
    for (size_t i = 0; i < n_heavy; ++i) {
        prepared.centroid += prepared.conformer->getAtomPos(static_cast<unsigned int>(i));
    }
    prepared.centroid *= 1.0 / static_cast<double>(n_heavy);
    for (size_t i = 0; i < n_heavy; ++i) {
        const auto& p = prepared.conformer->getAtomPos(static_cast<unsigned int>(i));
        const double dx = p.x - prepared.centroid.x;
        const double dy = p.y - prepared.centroid.y;
        const double dz = p.z - prepared.centroid.z;
        prepared.radii[i] = std::sqrt(dx * dx + dy * dy + dz * dz);
    }
    return prepared;
}

bool ConformerCache::isSimilarFast(const PreparedHeavyConformer& probe,
                                   double rmsd_threshold) const {
    if (rmsd_threshold == 0.0 || !cache_mol_no_h || cache_mol_no_h->getNumConformers() == 0) {
        return false;
    }
    try {
        return rmsd_calculator && rmsd_calculator->isSimilarToAny(
            *cache_mol_no_h, *probe.conformer, probe.radii, rmsd_threshold,
            &cache_ref_radii, &cache_ref_centroids);
    } catch (const std::exception& e) {
        std::cerr << "❌ Error in ConformerCache::isSimilarFast: " << e.what() << std::endl;
    }
    return false;
}

void ConformerCache::addPreparedConformer(PreparedHeavyConformer&& prepared,
                                          double energy) {
    if (!cache_mol_no_h) throw std::runtime_error("ConformerCache not initialized");
    if (!prepared.conformer) throw std::invalid_argument("Prepared conformer has no coordinates");

    cache_mol_no_h->addConformer(prepared.conformer.release(), true);
    energies.push_back(energy);
    cache_ref_radii.push_back(std::move(prepared.radii));
    cache_ref_centroids.push_back(prepared.centroid);
}

size_t ConformerCache::size() const {
    return cache_mol_no_h ? cache_mol_no_h->getNumConformers() : 0;
}

void ConformerCache::clear() {
    if (cache_mol_no_h) {
        cache_mol_no_h->clearConformers();
    }
    energies.clear();
    cache_ref_radii.clear();
    cache_ref_centroids.clear();
}

// ========================
// RANDOM GENERATOR IMPLEMENTATION
// ========================

int RandomGenerator::weightedChoice(const std::vector<double>& weights) {
    // OPTIMIZATION: Check if weights are uniform (all equal)
    if (!weights.empty()) {
        double first_weight = weights[0];
        bool all_equal = std::all_of(weights.begin(), weights.end(), 
                                   [first_weight](double w) { return std::abs(w - first_weight) < 1e-9; });
        
        if (all_equal) {
            // Use fast uniform_int_distribution for uniform weights
            size_t size = weights.size();
            auto it = uniform_distributions.find(size);
            if (it == uniform_distributions.end()) {
                auto dist = std::make_unique<std::uniform_int_distribution<>>(0, static_cast<int>(size - 1));
                auto result = dist->operator()(gen);
                uniform_distributions[size] = std::move(dist);
                return result;
            } else {
                return it->second->operator()(gen);
            }
        }
    }
    
    // Fall back to discrete_distribution for non-uniform weights
    auto it = cached_distributions.find(weights);
    if (it == cached_distributions.end()) {
        auto dist = std::make_unique<std::discrete_distribution<>>(weights.begin(), weights.end());
        auto result = dist->operator()(gen);
        cached_distributions[weights] = std::move(dist);
        return result;
    } else {
        return it->second->operator()(gen);
    }
}

std::vector<int> RandomGenerator::weightedChoices(const std::vector<double>& weights, int k) {
    std::vector<int> result;
    result.reserve(k);
    
    // Use the optimized weightedChoice for each selection
    for (int i = 0; i < k; ++i) {
        result.push_back(weightedChoice(weights));
    }
    
    return result;
}

double RandomGenerator::getRandomAngle(double center, double tolerance, const std::string& method) {
    if (tolerance == 0.0) return center;
    
    double angle;
    do {
        if (method == "gauss") {
            angle = normal_dist(gen) * tolerance + center;
        } else { // uniform
            angle = uniform_real_dist(gen) * 2 * tolerance - tolerance + center;
        }
    } while (angle < center - tolerance || angle > center + tolerance);
    
    // Normalize to [-180, 180] range
    return fmod(angle + 180.0, 360.0) - 180.0;
}

// ========================
// SAMPLING UTILS IMPLEMENTATION
// ========================

std::pair<AtomPairs, AtomPairs> SamplingUtils::precomputeBondedAndSameParentPairs(const RDKit::ROMol& mol) {
    AtomPairs bonded_pairs;
    AtomPairs same_parent_pairs;
    
    // Iterate over all atoms in the molecule
    for (const auto& atom : mol.atoms()) {
        int atom_idx = atom->getIdx();
        std::vector<int> neighbor_indices;
        
        // Get neighbors and collect their indices
        for (const auto& neighbor : mol.atomNeighbors(atom)) {
            int neighbor_idx = neighbor->getIdx();
            neighbor_indices.push_back(neighbor_idx);
            
            // Add bonded pairs in both directions (like Python version)
            bonded_pairs.emplace_back(atom_idx, neighbor_idx);
            bonded_pairs.emplace_back(neighbor_idx, atom_idx);
        }
        
        // Get atoms that share the same parent (common neighbors)
        if (neighbor_indices.size() > 1) {
            for (size_t i = 0; i < neighbor_indices.size(); ++i) {
                for (size_t j = i + 1; j < neighbor_indices.size(); ++j) {
                    // Add same parent pairs in both directions (like Python version)
                    same_parent_pairs.emplace_back(neighbor_indices[i], neighbor_indices[j]);
                    same_parent_pairs.emplace_back(neighbor_indices[j], neighbor_indices[i]);
                }
            }
        }
    }
    
    return std::make_pair(bonded_pairs, same_parent_pairs);
}

NonbondedClashPairs SamplingUtils::precomputeNonbondedPairs(const RDKit::ROMol& mol,
                                                            const AtomPairs& bonded_pairs,
                                                            const AtomPairs& same_parent_pairs,
                                                            double clash_scale) {
    NonbondedClashPairs candidate_pairs;
    std::set<std::pair<int, int>> excluded_pairs;
    const auto* periodic_table = RDKit::PeriodicTable::getTable();
    const auto canonical_pair = [](int atom1, int atom2) {
        return std::make_pair(std::min(atom1, atom2), std::max(atom1, atom2));
    };
    
    // Add bonded and same parent pairs to exclusion set
    for (const auto& pair : bonded_pairs) {
        excluded_pairs.insert(canonical_pair(pair.first, pair.second));
    }
    for (const auto& pair : same_parent_pairs) {
        excluded_pairs.insert(canonical_pair(pair.first, pair.second));
    }
    
    // Generate all possible pairs and exclude bonded/same parent pairs
    int num_atoms = mol.getNumAtoms();
    for (int i = 0; i < num_atoms; ++i) {
        for (int j = i + 1; j < num_atoms; ++j) {
            std::pair<int, int> pair_ij = {i, j};
            if (excluded_pairs.find(pair_ij) == excluded_pairs.end()) {
                const double vdw_sum = periodic_table->getRvdw(mol.getAtomWithIdx(i)->getAtomicNum()) +
                                       periodic_table->getRvdw(mol.getAtomWithIdx(j)->getAtomicNum());
                const double cutoff = clash_scale * vdw_sum;
                candidate_pairs.push_back({i, j, cutoff * cutoff});
            }
        }
    }
    
    return candidate_pairs;
}

bool SamplingUtils::checkTooCloseNonbondedAtoms(const RDKit::Conformer& conf,
                                                const NonbondedClashPairs& candidate_pairs) {
    for (const auto& pair : candidate_pairs) {
        const RDGeom::Point3D& pos1 = conf.getAtomPos(pair.atom1);
        const RDGeom::Point3D& pos2 = conf.getAtomPos(pair.atom2);
        
        // Calculate squared distance (avoid sqrt)
        double dx = pos1.x - pos2.x;
        double dy = pos1.y - pos2.y;
        double dz = pos1.z - pos2.z;
        double distance_sq = dx*dx + dy*dy + dz*dz;
        
        if (distance_sq < pair.cutoff_sq) {
            return true;  // Found clash
        }
    }
    return false;
}

bool SamplingUtils::checkTimeout(std::chrono::steady_clock::time_point start_time,
                                int timeout_seconds) {
    auto current_time = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time);
    return elapsed.count() >= timeout_seconds;
}

std::vector<double> SamplingUtils::extractCoreAngles(const std::vector<double>& full_angles,
                                                     int num_hetero_H_bonds) {
    std::vector<double> core_angles = full_angles;
    
    if (num_hetero_H_bonds > 0 && core_angles.size() > static_cast<size_t>(num_hetero_H_bonds)) {
        core_angles.resize(core_angles.size() - num_hetero_H_bonds);
    } else if (num_hetero_H_bonds > 0 && core_angles.size() == static_cast<size_t>(num_hetero_H_bonds)) {
        // Edge case: core_angles would become empty after resize
        // This represents a valid empty core state that should allow multiple hetero-H orientations
        core_angles.clear();
    }
    
    return core_angles;
}

// Check if a conformer (vector of angles) is similar to any in a list, within a given tolerance (deg)
bool SamplingUtils::isSimilarConformer(const std::vector<double>& angles1,
                                       const std::vector<std::vector<double>>& visited,
                                       double tolerance) {
    for (const auto& angles2 : visited) {
        if (angles1.size() != angles2.size()) continue;
        bool similar = true;
        for (size_t i = 0; i < angles1.size(); ++i) {
            double diff = std::abs(angles1[i] - angles2[i]);
            // Handle angle wrapping
            if (diff > 180.0) diff = 360.0 - diff;
            if (diff > tolerance) {
                similar = false;
                break;
            }
        }
        if (similar) return true;
    }
    return false;
}
} // namespace StochasticSampling
