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
#include <random>
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
        //"[O,N;D1;$([O,N;D1]-[*]=[O,N;D1]),$([O,N;D1]=[*]-[O,N;D1])]~[*]"; Removed the non-terminal P(=O)[O-] group to avoid over-symmetrization
    
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
    
    // Debug output
    // fprintf(stderr, "DEBUG: Generated %zu symmetric mappings for molecule\n", symmetric_mappings_.size());
    // for (size_t i = 0; i < symmetric_mappings_.size(); ++i) {
    //     fprintf(stderr, "  Mapping %zu: %zu atom pairs\n", i, symmetric_mappings_[i].size());
    // }
} 

SameMoleculeRMSDCalculator::SameMoleculeRMSDCalculator(bool use_symmetry, bool symmetrize_conjugated_terminal_groups) 
    : initialized_(false), use_symmetry_(use_symmetry), symmetrize_conjugated_terminal_groups_(symmetrize_conjugated_terminal_groups), num_heavy_atoms_(0) {}

void SameMoleculeRMSDCalculator::initialize(const RDKit::ROMol& mol) {
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
    
    // Reserve capacity for point vectors to avoid repeated allocations
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

    // Precompute heavy atom positions and radii once
    const size_t n_heavy = num_heavy_atoms_;
    std::vector<const RDGeom::Point3D*> heavy_probe_pos(n_heavy), heavy_ref_pos(n_heavy);
    // Populate pointers to heavy atom positions (no centroid accumulation here — we assume cached centroids/radii exist)
    for (size_t i = 0; i < n_heavy; ++i) {
        int atom_idx = heavy_atom_indices_[static_cast<int>(i)];
        heavy_probe_pos[i] = &probe_conf.getAtomPos(atom_idx);
        heavy_ref_pos[i]   = &ref_conf.getAtomPos(atom_idx);
    }
    // Precompute radii (distance from centroid)
    std::vector<double> probe_radii(n_heavy), ref_radii(n_heavy);
    // Use provided probe radii (distance from centroid) — caller guarantees this is populated
    for (size_t i = 0; i < n_heavy; ++i) {
        probe_radii[i] = (*cached_probe_radii)[i];
        ref_radii[i] = (*cached_ref_radii)[i];
    }

    // Compute heuristic score (radial lower bound) for each mapping
    struct ScoredMapping {
        double heuristic;
        const std::vector<std::pair<int,int>>* mapping;
    };
    std::vector<ScoredMapping> scored;
    scored.reserve(symmetric_mappings_.size());

    for (const auto& mapping : symmetric_mappings_) {
        double sum_sq = 0.0;
        for (auto& pair : mapping) {
            double d = probe_radii[pair.first] - ref_radii[pair.second];
            sum_sq += d * d;
        }
        double heuristic = (mapping.empty() ? 0.0 : sum_sq / mapping.size());
        scored.push_back({heuristic, &mapping});
    }

    // Sort by heuristic (smaller = more similar)
    std::sort(scored.begin(), scored.end(),
              [](const ScoredMapping& a, const ScoredMapping& b) {
                  return a.heuristic < b.heuristic;
              });

    // Try each symmetric mapping
    for (const auto& sm : scored) {
        const auto& mapping = *sm.mapping;
        ref_points_.clear();
        probe_points_.clear();

        double lower_bound_msd = sm.heuristic; // use heuristic directly

        if (lower_bound_msd >= best_msd) { continue; }
        if (thres2 >= 0.0 && lower_bound_msd > thres2) { continue; }

        for (auto& pair : mapping) {
            probe_points_.push_back(heavy_probe_pos[pair.first]);
            ref_points_.push_back(heavy_ref_pos[pair.second]);
        }

        RDGeom::Transform3D trans;
        double ssr = RDNumeric::Alignments::AlignPoints(
            ref_points_, probe_points_, trans, nullptr, false, 25);

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
                                               int probe_conf_id,
                                               const std::vector<int>& ref_conf_ids,
                                               double rmsd_threshold,
                                               const std::vector<std::vector<double>>* cached_ref_radii,
                                               const std::vector<RDGeom::Point3D>* cached_ref_centroids) const {
    if (!initialized_) {
        throw std::runtime_error("Calculator not initialized");
    }
    
    const RDKit::Conformer& probe_conf = mol.getConformer(probe_conf_id);
    double thres2 = (rmsd_threshold >= 0.0) ? rmsd_threshold * rmsd_threshold : -1.0;
    // Precompute probe centroid and squared radii once for the probe and reuse for heuristics and alignments
    const size_t n_heavy = num_heavy_atoms_;
    std::vector<double> probe_radii(n_heavy);
    RDGeom::Point3D probe_centroid_all(0.0, 0.0, 0.0);
    for (size_t i = 0; i < n_heavy; ++i) {
        int atom_idx = heavy_atom_indices_[static_cast<int>(i)];
        const RDGeom::Point3D &p = probe_conf.getAtomPos(atom_idx);
        probe_centroid_all += p;
    }
    probe_centroid_all *= (1.0 / static_cast<double>(n_heavy));
    for (size_t i = 0; i < n_heavy; ++i) {
        int atom_idx = heavy_atom_indices_[static_cast<int>(i)];
        const RDGeom::Point3D &p = probe_conf.getAtomPos(atom_idx);
        double dx = p.x - probe_centroid_all.x;
        double dy = p.y - probe_centroid_all.y;
        double dz = p.z - probe_centroid_all.z;
        probe_radii[i] = std::sqrt(dx*dx + dy*dy + dz*dz);
    }

    // Build candidate list with heuristic (mean squared difference of radii)
    struct Candidate { int conf_id; double heuristic; };
    std::vector<Candidate> candidates;
    candidates.reserve(ref_conf_ids.size());

    for (int ref_id : ref_conf_ids) {
        // Assume cached_ref_radii is provided and contains squared radii for each reference conformer
        const auto &ref_r = (*cached_ref_radii)[static_cast<size_t>(ref_id)];
        double sum_sq = 0.0;
        for (size_t i = 0; i < n_heavy; ++i) {
            double d = probe_radii[i] - ref_r[i];
            sum_sq += d * d;
        }
        double heuristic = sum_sq / static_cast<double>(n_heavy);
        candidates.push_back({ref_id, heuristic});
    }

    std::sort(candidates.begin(), candidates.end(), [](const Candidate &a, const Candidate &b){ return a.heuristic < b.heuristic; });

    // Evaluate candidates in heuristic order (most promising first)
    for (const auto &cand : candidates) {
        int ref_id = cand.conf_id;
        const RDKit::Conformer& ref_conf = mol.getConformer(ref_id);
    // Directly get references into the cached per-ref arrays (assumed present)
    const std::vector<double>* ref_radii_ptr = &((*cached_ref_radii)[static_cast<size_t>(ref_id)]);
    const RDGeom::Point3D* ref_centroid_ptr = &((*cached_ref_centroids)[static_cast<size_t>(ref_id)]);

    // Forward the precomputed probe squared radii to avoid recomputation inside calculateAlignedRMSD
    double msd = calculateAlignedRMSD(probe_conf, ref_conf, nullptr, thres2, ref_radii_ptr, ref_centroid_ptr, &probe_radii);
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
