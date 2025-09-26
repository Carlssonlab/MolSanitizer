//
// Accelerated RMSD Calculation for Same Molecular Graph
// Header file for optimized conformer comparison within the same molecule
//

#pragma once

#include <GraphMol/ROMol.h>
#include <GraphMol/Conformer.h>
#include <Geometry/point.h>
#include <Numerics/Alignment/AlignPoints.h>
#include <Numerics/Vector.h>
#include <Geometry/Transform3D.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <vector>
#include <memory>
#include <algorithm>
#include <cmath>
#include <limits>

namespace StochasticSampling {
namespace AcceleratedRMSD {

/**
 * Optimized RMSD calculator for conformers of the same molecular graph
 * Pre-computes mappings and atom indices to avoid repeated calculations
 */
class SameMoleculeRMSDCalculator {
private:
    std::vector<int> heavy_atom_indices_;
    std::vector<std::vector<std::pair<int, int>>> symmetric_mappings_;
    bool initialized_;
    bool use_symmetry_;
    bool symmetrize_conjugated_terminal_groups_;
    size_t num_heavy_atoms_;
    
    // Reusable point vectors to avoid repeated allocations
    mutable RDGeom::Point3DConstPtrVect ref_points_;
    mutable RDGeom::Point3DConstPtrVect probe_points_;
    
    /**
     * Generate all symmetric mappings for the same molecule
     * Uses RDKit's SubstructMatch to find equivalent atoms, with optional 
     * symmetrization of conjugated terminal groups (like COO-, NO2-)
     */
    void generateSymmetricMappings(const RDKit::ROMol& mol);
    
public:
    SameMoleculeRMSDCalculator(bool use_symmetry = true, bool symmetrize_conjugated_terminal_groups = true);
    
    /**
     * Initialize the calculator with a reference molecule
     * This pre-computes heavy atom indices and symmetric mappings
     */
    void initialize(const RDKit::ROMol& mol);
    
    /**
     * Calculate RMSD with optimal alignment considering symmetry
     * Tries all symmetric mappings and returns the best (lowest) RMSD
     */
    double calculateAlignedRMSD(const RDKit::Conformer& probe_conf,
                                   const RDKit::Conformer& ref_conf,
                                   RDGeom::Transform3D* transform = nullptr,
                                   double rmsd_threshold = -1.0) const;
     
    /**
     * Check if a conformer is similar to any in a set (below threshold)
     * Optimized version of the conformer filtering logic
     */
    bool isSimilarToAny(const RDKit::ROMol& mol,
                       int probe_conf_id,
                       const std::vector<int>& ref_conf_ids,
                       double rmsd_threshold) const;
    
    /**
     * Get the number of heavy atoms being used for RMSD calculation
     */
    size_t getNumHeavyAtoms() const;
    
    /**
     * Get the heavy atom indices
     */
    const std::vector<int>& getHeavyAtomIndices() const;
    
    /**
     * Get the number of symmetric mappings found
     */
    size_t getNumSymmetricMappings() const;
    
    /**
     * Check if symmetry handling is enabled
     */
    bool isSymmetryEnabled() const;
};

/**
 * Optimized conformer similarity checker
 * Replacement for the ConformerCache::isSimilarFast function
 */
class OptimizedConformerCache {
private:
    std::unique_ptr<SameMoleculeRMSDCalculator> rmsd_calc_;
    std::vector<int> cached_conf_ids_;
    
public:
    OptimizedConformerCache();
    
    void initialize(const RDKit::ROMol& mol);
    
    bool isSimilar(const RDKit::ROMol& mol,
                  int probe_conf_id,
                  double rmsd_threshold) const;
    
    void addConformer(int conf_id);
    
    size_t getNumCachedConformers() const;
    
    void clear();
};

} // namespace AcceleratedRMSD
} // namespace StochasticSampling