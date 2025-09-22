//
// Accelerated RMSD Calculation for Same Molecular Graph
// Optimized for conformer comparison within the same molecule
//

#include "accelerated_rmsd.h"
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/QueryBond.h>
#include <stdexcept>

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
        "[O,N;D1;$([O,N;D1]-[*]=[O,N;D1]),$([O,N;D1]=[*]-[O,N;D1])]~[*]";
    
    try {
        std::unique_ptr<RDKit::ROMol> qry(RDKit::SmartsToMol(qsmarts));
        if (!qry) {
            return; // Failed to parse SMARTS, skip symmetrization
        }
        
        auto matches = RDKit::SubstructMatch(mol, *qry);
        if (matches.empty()) {
            return; // No conjugated terminal groups found
        }
        
        // For each match, we need to identify the parent heavy atom and terminal atoms
        for (const auto& match : matches) {
            if (match.size() >= 2) { // Should have at least terminal and parent atom
                int terminal_idx = match[0].second;  // Terminal O/N atom
                int parent_idx = match[1].second;    // Parent atom
                
                // Find all terminal O/N atoms connected to this parent
                std::vector<int> terminal_atoms;
                for (const auto& bond : mol.atomBonds(mol.getAtomWithIdx(parent_idx))) {
                    int neighbor_idx = bond->getOtherAtomIdx(parent_idx);
                    const RDKit::Atom* neighbor = mol.getAtomWithIdx(neighbor_idx);
                    
                    // Check if it's a terminal O or N (degree 1)
                    if ((neighbor->getAtomicNum() == 8 || neighbor->getAtomicNum() == 7) && 
                        neighbor->getDegree() == 1) {
                        terminal_atoms.push_back(neighbor_idx);
                    }
                }
                
                // If we have multiple terminal atoms, mark them as equivalent
                if (terminal_atoms.size() > 1) {
                    // Set the same atom map number for equivalent atoms
                    int map_num = terminal_atoms[0] + 1000; // Arbitrary offset to avoid conflicts
                    for (int atom_idx : terminal_atoms) {
                        mol.getAtomWithIdx(atom_idx)->setAtomMapNum(map_num);
                    }
                }
            }
        }
    } catch (const std::exception& e) {
        // If symmetrization fails, continue without it
        return;
    }
}

void SameMoleculeRMSDCalculator::generateSymmetricMappings(const RDKit::ROMol& mol) {
    symmetric_mappings_.clear();
    
    // Create a copy of the molecule for symmetrization
    RDKit::RWMol mol_copy(mol);
    symmetrizeTerminalAtoms(mol_copy);
    
    // Get all substructure matches of the molecule against itself
    std::vector<RDKit::MatchVectType> all_matches;
    bool uniquify = false;
    bool recursionPossible = true;
    bool useChirality = false;
    bool useQueryQueryMatches = false;
    int maxMatches = 1000; // Limit to prevent excessive computation
    
    try {
        RDKit::SubstructMatch(mol_copy, mol_copy, all_matches, uniquify, recursionPossible,
                             useChirality, useQueryQueryMatches, maxMatches);
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
        
        if (!heavy_mapping.empty() && heavy_mapping.size() == heavy_atom_indices_.size()) {
            symmetric_mappings_.push_back(heavy_mapping);
        }
    }
    
    // If no valid mappings found, create identity mapping
    if (symmetric_mappings_.empty()) {
        std::vector<std::pair<int, int>> identity_mapping;
        for (size_t i = 0; i < heavy_atom_indices_.size(); ++i) {
            identity_mapping.emplace_back(i, i);
        }
        symmetric_mappings_.push_back(identity_mapping);
    }
}

SameMoleculeRMSDCalculator::SameMoleculeRMSDCalculator(bool use_symmetry) 
    : initialized_(false), use_symmetry_(use_symmetry) {}

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
    
    // Generate symmetric mappings if enabled
    if (use_symmetry_) {
        generateSymmetricMappings(mol);
    } else {
        // Create identity mapping only
        std::vector<std::pair<int, int>> identity_mapping;
        for (size_t i = 0; i < heavy_atom_indices_.size(); ++i) {
            identity_mapping.emplace_back(i, i);
        }
        symmetric_mappings_.push_back(identity_mapping);
    }
    
    initialized_ = true;
}

double SameMoleculeRMSDCalculator::calculateAlignedRMSD(const RDKit::Conformer& probe_conf, 
                                                       const RDKit::Conformer& ref_conf,
                                                       RDGeom::Transform3D* transform) const {
    if (!initialized_) {
        throw std::runtime_error("Calculator not initialized");
    }
    
    double best_rmsd = std::numeric_limits<double>::max();
    RDGeom::Transform3D best_transform;
    
    // Try each symmetric mapping and find the one with lowest RMSD
    for (const auto& mapping : symmetric_mappings_) {
        // Prepare point arrays for this mapping
        RDGeom::Point3DConstPtrVect ref_points, probe_points;
        
        for (const auto& pair : mapping) {
            int probe_heavy_idx = pair.first;
            int ref_heavy_idx = pair.second;
            
            // Get actual atom indices
            int probe_atom_idx = heavy_atom_indices_[probe_heavy_idx];
            int ref_atom_idx = heavy_atom_indices_[ref_heavy_idx];
            
            probe_points.push_back(&probe_conf.getAtomPos(probe_atom_idx));
            ref_points.push_back(&ref_conf.getAtomPos(ref_atom_idx));
        }
        
        // Perform alignment for this mapping
        RDGeom::Transform3D trans;
        double ssr = RDNumeric::Alignments::AlignPoints(
            ref_points, probe_points, trans, nullptr, false, 50);
        
        double rmsd = std::sqrt(ssr / heavy_atom_indices_.size());
        
        if (rmsd < best_rmsd) {
            best_rmsd = rmsd;
            // Copy transform data manually since assignment operator is deleted
            for (unsigned int i = 0; i < 4; ++i) {
                for (unsigned int j = 0; j < 4; ++j) {
                    best_transform.setVal(i, j, trans.getVal(i, j));
                }
            }
        }
    }
    
    if (transform) {
        // Copy best transform data manually since assignment operator is deleted
        for (unsigned int i = 0; i < 4; ++i) {
            for (unsigned int j = 0; j < 4; ++j) {
                transform->setVal(i, j, best_transform.getVal(i, j));
            }
        }
    }
    
    return best_rmsd;
}

std::vector<double> SameMoleculeRMSDCalculator::calculateBatchRMSD(const RDKit::ROMol& mol,
                                                                   int ref_conf_id,
                                                                   const std::vector<int>& probe_conf_ids) const {
    if (!initialized_) {
        throw std::runtime_error("Calculator not initialized");
    }
    
    std::vector<double> rmsd_values;
    rmsd_values.reserve(probe_conf_ids.size());
    
    const RDKit::Conformer& ref_conf = mol.getConformer(ref_conf_id);
    
    for (int probe_id : probe_conf_ids) {
        const RDKit::Conformer& probe_conf = mol.getConformer(probe_id);
        double rmsd = calculateAlignedRMSD(probe_conf, ref_conf);
        rmsd_values.push_back(rmsd);
    }
    
    return rmsd_values;
}

double SameMoleculeRMSDCalculator::getBestRMSDSameMolecule(const RDKit::ROMol& mol,
                                                          int probe_conf_id,
                                                          const std::vector<int>& ref_conf_ids,
                                                          int* best_ref_id) const {
    if (!initialized_) {
        throw std::runtime_error("Calculator not initialized");
    }
    
    double best_rmsd = std::numeric_limits<double>::max();
    int best_id = -1;
    
    const RDKit::Conformer& probe_conf = mol.getConformer(probe_conf_id);
    
    for (int ref_id : ref_conf_ids) {
        const RDKit::Conformer& ref_conf = mol.getConformer(ref_id);
        double rmsd = calculateAlignedRMSD(probe_conf, ref_conf);
        
        if (rmsd < best_rmsd) {
            best_rmsd = rmsd;
            best_id = ref_id;
        }
    }
    
    if (best_ref_id) {
        *best_ref_id = best_id;
    }
    
    return best_rmsd;
}

bool SameMoleculeRMSDCalculator::isSimilarToAny(const RDKit::ROMol& mol,
                                               int probe_conf_id,
                                               const std::vector<int>& ref_conf_ids,
                                               double rmsd_threshold) const {
    if (!initialized_) {
        throw std::runtime_error("Calculator not initialized");
    }
    
    const RDKit::Conformer& probe_conf = mol.getConformer(probe_conf_id);
    
    for (int ref_id : ref_conf_ids) {
        const RDKit::Conformer& ref_conf = mol.getConformer(ref_id);
        double rmsd = calculateAlignedRMSD(probe_conf, ref_conf);
        
        if (rmsd <= rmsd_threshold) {
            return true;
        }
    }
    
    return false;
}

size_t SameMoleculeRMSDCalculator::getNumHeavyAtoms() const {
    return heavy_atom_indices_.size();
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

double calculateFastAlignedRMSD(const RDKit::ROMol& mol,
                               int conf_id1,
                               int conf_id2,
                               bool heavy_atoms_only) {
    const RDKit::Conformer& conf1 = mol.getConformer(conf_id1);
    const RDKit::Conformer& conf2 = mol.getConformer(conf_id2);
    
    // Prepare point arrays for alignment
    RDGeom::Point3DConstPtrVect ref_points, probe_points;
    
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
        // Skip hydrogens if heavy_atoms_only is true
        if (heavy_atoms_only && mol.getAtomWithIdx(i)->getAtomicNum() == 1) {
            continue;
        }
        
        probe_points.push_back(&conf1.getAtomPos(i));
        ref_points.push_back(&conf2.getAtomPos(i));
    }
    
    // Perform alignment
    RDGeom::Transform3D trans;
    double ssr = RDNumeric::Alignments::AlignPoints(
        ref_points, probe_points, trans, nullptr, false, 50);
    
    return std::sqrt(ssr / probe_points.size());
}

OptimizedConformerCache::OptimizedConformerCache() 
    : rmsd_calc_(std::make_unique<SameMoleculeRMSDCalculator>()) {}

void OptimizedConformerCache::initialize(const RDKit::ROMol& mol) {
    rmsd_calc_->initialize(mol);
    cached_conf_ids_.clear();
}

bool OptimizedConformerCache::isSimilar(const RDKit::ROMol& mol,
                                       int probe_conf_id,
                                       double rmsd_threshold) const {
    if (cached_conf_ids_.empty()) {
        return false;
    }
    
    return rmsd_calc_->isSimilarToAny(mol, probe_conf_id, cached_conf_ids_, 
                                     rmsd_threshold);
}

void OptimizedConformerCache::addConformer(int conf_id) {
    cached_conf_ids_.push_back(conf_id);
}

size_t OptimizedConformerCache::getNumCachedConformers() const {
    return cached_conf_ids_.size();
}

void OptimizedConformerCache::clear() {
    cached_conf_ids_.clear();
}

} // namespace AcceleratedRMSD
} // namespace StochasticSampling