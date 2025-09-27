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
    fprintf(stderr, "DEBUG: Generated %zu symmetric mappings for molecule\n", symmetric_mappings_.size());
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

double SameMoleculeRMSDCalculator::calculateAlignedRMSD(const RDKit::Conformer& probe_conf, 
                                                       const RDKit::Conformer& ref_conf,
                                                       RDGeom::Transform3D* transform,
                                                       double thres2) const {
    if (!initialized_) {
        throw std::runtime_error("Calculator not initialized");
    }
    
    double best_msd = std::numeric_limits<double>::max();
    
    // Create a shuffled copy of symmetric mappings for stochastic early exit
    std::vector<std::vector<std::pair<int, int>>> shuffled_mappings = symmetric_mappings_;
    std::shuffle(shuffled_mappings.begin(), shuffled_mappings.end(), std::mt19937{std::random_device{}()});
    
    // Precompute heavy-atom positions, centroids and radii once to avoid repeated work
    const size_t n_heavy = num_heavy_atoms_;
    std::vector<const RDGeom::Point3D*> heavy_probe_pos(n_heavy), heavy_ref_pos(n_heavy);
    RDGeom::Point3D probe_centroid_all(0.0, 0.0, 0.0), ref_centroid_all(0.0, 0.0, 0.0);
    for (size_t i = 0; i < n_heavy; ++i) {
        int atom_idx = heavy_atom_indices_[static_cast<int>(i)];
        heavy_probe_pos[i] = &probe_conf.getAtomPos(atom_idx);
        heavy_ref_pos[i] = &ref_conf.getAtomPos(atom_idx);
        probe_centroid_all += *heavy_probe_pos[i];
        ref_centroid_all += *heavy_ref_pos[i];
    }
    double inv_n_heavy = 1.0 / static_cast<double>(n_heavy);
    probe_centroid_all *= inv_n_heavy;
    ref_centroid_all *= inv_n_heavy;

    // Precompute radii (distance to centroid) once per heavy atom
    std::vector<double> probe_radii(n_heavy), ref_radii(n_heavy);
    for (size_t i = 0; i < n_heavy; ++i) {
        const RDGeom::Point3D &pp = *heavy_probe_pos[i];
        const RDGeom::Point3D &rp = *heavy_ref_pos[i];
        double px = pp.x - probe_centroid_all.x;
        double py = pp.y - probe_centroid_all.y;
        double pz = pp.z - probe_centroid_all.z;
        double rx = rp.x - ref_centroid_all.x;
        double ry = rp.y - ref_centroid_all.y;
        double rz = rp.z - ref_centroid_all.z;
        probe_radii[i] = std::sqrt(px*px + py*py + pz*pz);
        ref_radii[i] = std::sqrt(rx*rx + ry*ry + rz*rz);
    }

    // Try each symmetric mapping and find the one with lowest RMSD
    for (const auto& mapping : shuffled_mappings) {
        // Reset point arrays without deallocating capacity
        ref_points_.resize(0);
        probe_points_.resize(0);
    
        const double npts = static_cast<double>(mapping.size());
        const double inv_n = (npts > 0.0) ? 1.0 / npts : inv_n_heavy;

        // Compute radial lower bound using precomputed radii (no sqrt per mapping)
        double radial_sq_sum = 0.0;
        for (const auto& pair : mapping) {
            int probe_heavy_idx = pair.first;
            int ref_heavy_idx = pair.second;

            probe_points_.push_back(heavy_probe_pos[probe_heavy_idx]);
            ref_points_.push_back(heavy_ref_pos[ref_heavy_idx]);

            double d = probe_radii[probe_heavy_idx] - ref_radii[ref_heavy_idx];
            radial_sq_sum += d * d;
        }

        double lower_bound_msd = radial_sq_sum * inv_n;

        // If radial lower bound already worse than best or threshold, skip expensive alignment
        if (lower_bound_msd >= best_msd) { continue; }
        if (thres2 >= 0.0 && lower_bound_msd > thres2) { continue; }

        // Perform alignment for this mapping
        RDGeom::Transform3D trans;
        double ssr = RDNumeric::Alignments::AlignPoints(
            ref_points_, probe_points_, trans, nullptr, false, 25);

        // Compare squared RMSD directly to avoid sqrt in hot loop
        double msd = ssr / static_cast<double>(num_heavy_atoms_);


        if (msd < best_msd) {
            best_msd = msd;
        }

        // Early exit if threshold is set and met
        if (thres2 >= 0.0 && msd <= thres2) {
            return msd;
        }
    }
    
    return best_msd;
}

bool SameMoleculeRMSDCalculator::isSimilarToAny(const RDKit::ROMol& mol,
                                               int probe_conf_id,
                                               const std::vector<int>& ref_conf_ids,
                                               double rmsd_threshold) const {
    if (!initialized_) {
        throw std::runtime_error("Calculator not initialized");
    }
    
    const RDKit::Conformer& probe_conf = mol.getConformer(probe_conf_id);
    double thres2 = (rmsd_threshold >= 0.0) ? rmsd_threshold * rmsd_threshold : -1.0;
    for (int ref_id : ref_conf_ids) {
        const RDKit::Conformer& ref_conf = mol.getConformer(ref_id);
        double msd = calculateAlignedRMSD(probe_conf, ref_conf, nullptr, thres2);
        if (msd <= thres2) {
            return true;
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