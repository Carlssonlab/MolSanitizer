//
// Accelerated RMSD Calculation for Same Molecular Graph
// Optimized for conformer comparison within the same molecule
//

#include "accelerated_rmsd.h"
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#if __cplusplus >= 202002L
#include <span>
#endif

namespace StochasticSampling {
namespace AcceleratedRMSD {

/**
 * Symmetrize conjugated terminal groups (like COO-, NO2-) to handle resonance structures
 * This implements the same logic as RDKit's details::symmetrizeTerminalAtoms
 */
void symmetrizeTerminalAtoms(RDKit::RWMol &mol) {
    // SMARTS pattern to identify terminal O,N atoms in conjugated systems
    // Matches patterns like C(=O)[O-] and C([O-])=O allowing flexible matching
    const std::string qsmarts = 
        "[O,N;D1;$([O,N;D1]-[*]=[O,N;D1]),$([O,N;D1]=[*]-[O,N;D1])]~[*;!$(P([*X{2-}])([*X{2-}])(=O)[O-])]";
    // Exclude non-terminal P(=O)[O-] groups to avoid over-symmetrization.
    
    try {
        std::unique_ptr<RDKit::ROMol> qry(RDKit::SmartsToMol(qsmarts));
        if (!qry) {
            return; // Failed to parse SMARTS, skip symmetrization
        }
        
        auto matches = RDKit::SubstructMatch(mol, *qry);
        if (matches.empty()) {
            return; // No conjugated terminal groups found
        }

        // Create a query bond that matches both single and double bonds
        RDKit::QueryBond qb;
        qb.setBondType(RDKit::Bond::SINGLE);
        // Add both single and double bond possibilities
        qb.setQuery(RDKit::makeBondOrderEqualsQuery(RDKit::Bond::SINGLE));
        qb.expandQuery(RDKit::makeBondOrderEqualsQuery(RDKit::Bond::DOUBLE), Queries::COMPOSITE_OR);
        
        for (const auto &match : matches) {
            // Neutralize formal charge on terminal atom
            mol.getAtomWithIdx(match[0].second)->setFormalCharge(0);
            
            // Replace the bond with a single-or-double query bond
            auto bond = mol.getBondBetweenAtoms(match[0].second, match[1].second);
            if (bond) {
                mol.replaceBond(bond->getIdx(), &qb);
            }
        }
        
    } catch (const std::exception& e) {
        // If symmetrization fails, continue without it
        // This maintains robustness while providing the feature when possible
    }
}

void SameMoleculeRMSDCalculator::generateSymmetricMappings(const RDKit::ROMol& mol) {
    symmetric_mappings_.clear();
    
    // Compute SO2 pairs
    std::vector<std::pair<int, int>> so2_pairs;
    for (const auto& atom : mol.atoms()) {
        if (atom->getAtomicNum() == 16) {  // Sulfur
            std::vector<int> oxygens;
            for (const auto& bond : mol.atomBonds(atom)) {
                const auto& neighbor = bond->getOtherAtom(atom);
                if (neighbor->getAtomicNum() == 8 && bond->getBondType() == RDKit::Bond::DOUBLE) {
                    oxygens.push_back(neighbor->getIdx());
                }
            }
            if (oxygens.size() == 2) {
                int o1 = std::min(oxygens[0], oxygens[1]);
                int o2 = std::max(oxygens[0], oxygens[1]);
                so2_pairs.emplace_back(o1, o2);
            }
        }
    }
    
    // Create a copy of the molecule for potential symmetrization
    std::unique_ptr<RDKit::RWMol> mol_for_match;
    const RDKit::ROMol* mol_to_use = &mol;
    
    if (symmetrize_conjugated_terminal_groups_) {
        mol_for_match.reset(new RDKit::RWMol(mol));
        symmetrizeTerminalAtoms(*mol_for_match);
        mol_to_use = mol_for_match.get();
    }
    
    // Get all substructure matches of the molecule against itself
    std::vector<RDKit::MatchVectType> all_matches;
    
    RDKit::SubstructMatchParameters params;
    params.uniquify = false;
    params.recursionPossible = true;
    params.useChirality = false;
    params.useQueryQueryMatches = false;
    params.maxMatches = 1000; // Limit to prevent excessive computation
    
    // Set extraFinalCheck to filter out invalid SO2 mappings
#if __cplusplus >= 202002L
    params.extraFinalCheck = [&so2_pairs](const RDKit::ROMol &mol, std::span<const unsigned int> match) -> bool {
#else
    params.extraFinalCheck = [&so2_pairs](const RDKit::ROMol &mol, const std::vector<unsigned int> &match) -> bool {
#endif
        for (const auto& pair : so2_pairs) {
            int o1 = pair.first;
            int o2 = pair.second;
            if (match[o1] == static_cast<unsigned int>(o2) && match[o2] == static_cast<unsigned int>(o1)) {
                return false;  // Invalid mapping
            }
        }
        return true;
    };
    
    try {
        all_matches = RDKit::SubstructMatch(mol, *mol_to_use, params);
    } catch (const std::exception& e) {
        // If SubstructMatch fails, fall back to identity mapping
        all_matches.clear();
        RDKit::MatchVectType identity_match;
        for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
            identity_match.emplace_back(i, i);
        }
        all_matches.push_back(identity_match);
    }
    
    // Convert matches to heavy-atom-only mappings
    for (const auto& match : all_matches) {
        std::vector<std::pair<int, int>> heavy_mapping;
        
        // Create mapping for heavy atoms only
        for (const auto& pair : match) {
            int probe_idx = pair.first;
            int ref_idx = pair.second;
            
            // Check if both atoms are heavy atoms
            if (mol.getAtomWithIdx(probe_idx)->getAtomicNum() != 1 && 
                mol.getAtomWithIdx(ref_idx)->getAtomicNum() != 1) {
                
                // Find the indices in the heavy_atom_indices_ array
                auto probe_it = std::find(heavy_atom_indices_.begin(), heavy_atom_indices_.end(), probe_idx);
                auto ref_it = std::find(heavy_atom_indices_.begin(), heavy_atom_indices_.end(), ref_idx);
                
                if (probe_it != heavy_atom_indices_.end() && ref_it != heavy_atom_indices_.end()) {
                    int probe_heavy_idx = std::distance(heavy_atom_indices_.begin(), probe_it);
                    int ref_heavy_idx = std::distance(heavy_atom_indices_.begin(), ref_it);
                    heavy_mapping.emplace_back(probe_heavy_idx, ref_heavy_idx);
                }
            }
        }
        
        if (!heavy_mapping.empty() && heavy_mapping.size() == num_heavy_atoms_) {
            symmetric_mappings_.push_back(heavy_mapping);
        }
    }
    
    // If no valid mappings found, create identity mapping
    if (symmetric_mappings_.empty()) {
        std::vector<std::pair<int, int>> identity_mapping;
        for (size_t i = 0; i < num_heavy_atoms_; ++i) {
            identity_mapping.emplace_back(i, i);
        }
        symmetric_mappings_.push_back(identity_mapping);
    }
    
} 

SameMoleculeRMSDCalculator::SameMoleculeRMSDCalculator(bool use_symmetry, bool symmetrize_conjugated_terminal_groups) 
    : initialized_(false), use_symmetry_(use_symmetry), symmetrize_conjugated_terminal_groups_(symmetrize_conjugated_terminal_groups), num_heavy_atoms_(0) {}

void SameMoleculeRMSDCalculator::initialize(const RDKit::ROMol& mol) {
    initialized_ = false;
    scored_mappings_.clear();
    candidates_.clear();
    heavy_atom_indices_.clear();
    symmetric_mappings_.clear();
    
    // Pre-compute heavy atom indices (exclude hydrogens)
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
        const RDKit::Atom* atom = mol.getAtomWithIdx(i);
        if (atom->getAtomicNum() != 1) { // Not hydrogen
            heavy_atom_indices_.push_back(static_cast<int>(i));
        }
    }
    
    num_heavy_atoms_ = heavy_atom_indices_.size();
    if (num_heavy_atoms_ == 0) {
        throw std::invalid_argument(
            "RMSD calculation requires at least one heavy atom");
    }
    
    // Size/reserve every calculateAlignedRMSD() scratch buffer up front.
    heavy_probe_positions_.resize(num_heavy_atoms_);
    heavy_ref_positions_.resize(num_heavy_atoms_);
    ref_points_.clear();
    probe_points_.clear();
    ref_points_.reserve(num_heavy_atoms_);
    probe_points_.reserve(num_heavy_atoms_);
    
    // Generate symmetric mappings if enabled
    if (use_symmetry_) {
        generateSymmetricMappings(mol);
    } else {
        // Create identity mapping only
        std::vector<std::pair<int, int>> identity_mapping;
        for (size_t i = 0; i < num_heavy_atoms_; ++i) {
            identity_mapping.emplace_back(i, i);
        }
        symmetric_mappings_.push_back(identity_mapping);
    }

    scored_mappings_.clear();
    scored_mappings_.reserve(symmetric_mappings_.size());
    
    initialized_ = true;
}


double SameMoleculeRMSDCalculator::calculateAlignedRMSD(
    const RDKit::Conformer& probe_conf,
    const RDKit::Conformer& ref_conf,
    RDGeom::Transform3D* transform,
    double thres2,
    const std::vector<double>* cached_ref_radii,
    const RDGeom::Point3D* cached_ref_centroid,
    const std::vector<double>* cached_probe_radii) const {
    if (!initialized_) {
        throw std::runtime_error("Calculator not initialized");
    }

    double best_msd = std::numeric_limits<double>::max();

    // Cache heavy-atom position pointers once per reference. Storage is owned
    // by the calculator so this does not allocate on the hot path.
    const size_t n_heavy = num_heavy_atoms_;
    const auto& probe_positions = probe_conf.getPositions();
    const auto& ref_positions = ref_conf.getPositions();
    if (cached_probe_radii == nullptr || cached_ref_radii == nullptr) {
        throw std::invalid_argument(
            "Precomputed probe and reference radii are required");
    }
    if (cached_probe_radii->size() < n_heavy ||
        cached_ref_radii->size() < n_heavy) {
        throw std::invalid_argument(
            "Precomputed radii do not cover every heavy atom");
    }
    if (!heavy_atom_indices_.empty()) {
        const size_t last_atom_index =
            static_cast<size_t>(heavy_atom_indices_.back());
        if (probe_positions.size() <= last_atom_index ||
            ref_positions.size() <= last_atom_index) {
            throw std::invalid_argument(
                "Conformer atom count does not match the initialized molecule");
        }
    }
    for (size_t i = 0; i < n_heavy; ++i) {
        int atom_idx = heavy_atom_indices_[static_cast<int>(i)];
        heavy_probe_positions_[i] = &probe_positions[atom_idx];
        heavy_ref_positions_[i] = &ref_positions[atom_idx];
    }

    // The caller owns descriptors for both conformers; alias them instead of
    // allocating and copying two radius arrays for every reference.
    const auto& probe_radii = *cached_probe_radii;
    const auto& ref_radii = *cached_ref_radii;

    // Compute heuristic score (radial lower bound) for each mapping
    scored_mappings_.clear();

    for (size_t mapping_index = 0;
         mapping_index < symmetric_mappings_.size();
         ++mapping_index) {
        const auto& mapping = symmetric_mappings_[mapping_index];
        double sum_sq = 0.0;
        for (auto& pair : mapping) {
            double d = probe_radii[pair.first] - ref_radii[pair.second];
            sum_sq += d * d;
        }
        double heuristic = (mapping.empty() ? 0.0 : sum_sq / mapping.size());
        scored_mappings_.push_back({heuristic, mapping_index});
    }

    // Sort by heuristic (smaller = more similar)
    std::sort(scored_mappings_.begin(), scored_mappings_.end(),
              [](const ScoredMapping& a, const ScoredMapping& b) {
                  return a.heuristic < b.heuristic;
              });

    // Try each symmetric mapping
    for (const auto& sm : scored_mappings_) {
        const auto& mapping = symmetric_mappings_[sm.mapping_index];
        ref_points_.clear();
        probe_points_.clear();

        double lower_bound_msd = sm.heuristic; // use heuristic directly

        if (lower_bound_msd >= best_msd) { continue; }
        if (thres2 >= 0.0 && lower_bound_msd > thres2) { continue; }

        for (auto& pair : mapping) {
            probe_points_.push_back(heavy_probe_positions_[pair.first]);
            ref_points_.push_back(heavy_ref_positions_[pair.second]);
        }

        double ssr = RDNumeric::Alignments::AlignPoints(
            ref_points_, probe_points_, alignment_transform_, nullptr, false, 25);

        double msd = ssr / static_cast<double>(num_heavy_atoms_);

        if (msd < best_msd) {
            best_msd = msd;
        }

        if (thres2 >= 0.0 && msd <= thres2) {
            return msd;
        }
    }

    return best_msd;
}

bool SameMoleculeRMSDCalculator::isSimilarToAny(const RDKit::ROMol& mol,
                                               const RDKit::Conformer& probe_conf,
                                               const std::vector<double>& probe_radii,
                                               double rmsd_threshold,
                                               const std::vector<std::vector<double>>* cached_ref_radii,
                                               const std::vector<RDGeom::Point3D>* cached_ref_centroids) const {
    if (!initialized_) {
        throw std::runtime_error("Calculator not initialized");
    }

    if (cached_ref_radii == nullptr) {
        throw std::invalid_argument("Cached reference radii are required");
    }
    if (probe_radii.size() < num_heavy_atoms_) {
        throw std::invalid_argument(
            "Probe radii do not cover every heavy atom");
    }
    if (cached_ref_radii->size() != mol.getNumConformers()) {
        throw std::invalid_argument(
            "Cached reference radii do not match the conformer cache");
    }
    // Centroids are retained in the API for compatibility, but radial
    // descriptors contain everything needed by the current lower bound.
    (void)cached_ref_centroids;
    
    double thres2 = (rmsd_threshold >= 0.0) ? rmsd_threshold * rmsd_threshold : -1.0;
    const size_t n_heavy = num_heavy_atoms_;

    // Build candidate list with heuristic (mean squared difference of radii)
    candidates_.clear();
    size_t cache_index = 0;
    for (auto conformer = mol.beginConformers();
         conformer != mol.endConformers();
         ++conformer, ++cache_index) {
        // Cached radii are indexed identically to the cache molecule's conformers.
        const auto& ref_r = (*cached_ref_radii)[cache_index];
        if (ref_r.size() < n_heavy) {
            throw std::invalid_argument(
                "Cached reference radii do not cover every heavy atom");
        }
        double sum_sq = 0.0;
        for (size_t i = 0; i < n_heavy; ++i) {
            double d = probe_radii[i] - ref_r[i];
            sum_sq += d * d;
        }
        double heuristic = sum_sq / static_cast<double>(n_heavy);
        candidates_.push_back({conformer->get(), cache_index, heuristic});
    }

    std::sort(candidates_.begin(), candidates_.end(),
              [](const Candidate& a, const Candidate& b) {
                  return a.heuristic < b.heuristic;
              });

    // Evaluate candidates in heuristic order (most promising first)
    for (const auto& candidate : candidates_) {
        const std::vector<double>* ref_radii_ptr =
            &((*cached_ref_radii)[candidate.cache_index]);

        // Forward both conformers' precomputed radial descriptors. Centroids
        // are not needed once those descriptors have been built.
        double msd = calculateAlignedRMSD(probe_conf, *candidate.conformer,
                                          nullptr, thres2,
                                          ref_radii_ptr, nullptr,
                                          &probe_radii);
        if (thres2 >= 0.0) {
            if (msd <= thres2) return true;
        } else {
            // If no threshold requested, any successful alignment is considered (msd finite)
            if (msd < std::numeric_limits<double>::infinity()) return true;
        }
    }
    
    return false;
}


size_t SameMoleculeRMSDCalculator::getNumHeavyAtoms() const {
    return num_heavy_atoms_;
}

const std::vector<int>& SameMoleculeRMSDCalculator::getHeavyAtomIndices() const {
    return heavy_atom_indices_;
}

size_t SameMoleculeRMSDCalculator::getNumSymmetricMappings() const {
    return symmetric_mappings_.size();
}

bool SameMoleculeRMSDCalculator::isSymmetryEnabled() const {
    return use_symmetry_;
}

} // namespace AcceleratedRMSD
} // namespace StochasticSampling
