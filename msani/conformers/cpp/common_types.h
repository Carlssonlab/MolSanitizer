#pragma once

#include <vector>
#include <map>
#include <set>
#include <tuple>
#include <random>
#include <chrono>
#include <memory>
#include <array>
#include <unordered_set>
#include <functional>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>

namespace StochasticSampling {

// Forward declarations
namespace AcceleratedRMSD {
    class SameMoleculeRMSDCalculator;
}

// ========================
// COMMON DATA STRUCTURES
// ========================

// Type aliases
using AtomPairs = std::vector<std::pair<int, int>>;
using ProductList = std::vector<struct ConformerResult>;
using HeteroBonds = std::vector<std::pair<int, int>>;

// Result structure for conformers
struct ConformerResult {
    RDKit::Conformer conformer;
    double energy;
    
    // Default constructor
    ConformerResult() : conformer(), energy(0.0) {}
    
    // Constructor that creates a standalone conformer (no molecule reference)
    ConformerResult(const RDKit::Conformer& conf, double e);

private:
    // Helper to create standalone conformer without molecule reference
    static RDKit::Conformer createStandaloneConformer(const RDKit::Conformer& source);
};

// ========================
// DISCRETE SAMPLING DATA STRUCTURES
// ========================

// Data structures for discrete stochastic sampling
struct AngleMapEntry {
    std::vector<int> dihedral_atoms;  // 4 atom indices defining the dihedral
    std::vector<double> possible_angles;  // Possible angle values
};

using AngleMap = std::map<int, AngleMapEntry>;
using ScoreMap = std::map<int, std::vector<double>>;
using ImportanceOrder = std::vector<double>;

// ========================
// CONTINUOUS SAMPLING DATA STRUCTURES
// ========================

// Simplified peak structure for continuous approach
struct TorsionPeak {
    double center;                    // Peak center angle
    std::vector<double> tolerance;    // Tolerance at different levels [strict, medium, relaxed]
    double weight;                    // Peak weight for selection
    
    TorsionPeak(double c, const std::vector<double>& tol, double w) 
        : center(c), tolerance(tol), weight(w) {}
};

// Bond information for continuous sampling
struct ContinuousTorsionBond {
    std::vector<int> dihedral_atoms;     // [atom1, atom2, atom3, atom4]
    std::vector<TorsionPeak> peaks;      // Available peaks for this bond
};

using ContinuousTorsionMap = std::map<int, ContinuousTorsionBond>;

// ========================
// CONFORMER CACHE
// ========================

// Cache for efficient RMSD checking with H-stripped conformers
class ConformerCache {
private:
    std::unique_ptr<RDKit::RWMol> cache_mol_no_h;  // Molecule without hydrogens for RDKit RMSD
    std::vector<double> energies;  // Track energies for each cached conformer
    std::unique_ptr<AcceleratedRMSD::SameMoleculeRMSDCalculator> rmsd_calculator;  // Accelerated RMSD calculator
    // Per-conformer cached data for heavy-atom-only conformers (per-atom radial distances from centroid)
    // cache_ref_radii[conf_id][atom_index] == radius (distance from that conformer's centroid)
    std::vector<std::vector<double>> cache_ref_radii; // radii per conformer (indexed by conformer id)
    std::vector<RDGeom::Point3D> cache_ref_centroids;  // centroid per conformer
    
public:
    ConformerCache() = default;
    
    // Initialize cache with reference molecule (will be stripped of Hs once)
    void initialize(const RDKit::ROMol& reference_mol);
    
    // FAST methods that avoid creating new molecules
    bool isSimilarFast(const RDKit::Conformer& conf, 
                       const RDKit::ROMol& mol_with_h,
                       const std::vector<int>& heavy_atom_mapping,
                       double rmsd_threshold) const;
    
    void addConformerFast(const RDKit::Conformer& conf,
                          const RDKit::ROMol& mol_with_h,
                          const std::vector<int>& heavy_atom_mapping,
                          double energy);
    
    // Get number of cached conformers
    size_t size() const;
    
    // Clear cache
    void clear();
};

// ========================
// RANDOM GENERATOR
// ========================

// Optimized random generator with caching for uniform vs non-uniform weights
class RandomGenerator {
private:
    std::mt19937 gen;
    
    // Cache for uniform distributions (size -> distribution)
    std::unordered_map<size_t, std::unique_ptr<std::uniform_int_distribution<>>> uniform_distributions;
    
    // Cache for weighted distributions (weights -> distribution)
    std::map<std::vector<double>, std::unique_ptr<std::discrete_distribution<>>> cached_distributions;
    
public:
    explicit RandomGenerator(unsigned seed = std::chrono::steady_clock::now().time_since_epoch().count()) 
        : gen(seed) {}
    
    // Optimized weighted choice with caching
    int weightedChoice(const std::vector<double>& weights);
    
    // Multiple weighted choices
    std::vector<int> weightedChoices(const std::vector<double>& weights, int k);
    
    // Shuffle vector
    template<typename T>
    void shuffle(std::vector<T>& vec) {
        std::shuffle(vec.begin(), vec.end(), gen);
    }
    
    // Random integer in range [min, max]
    int randint(int min, int max) {
        std::uniform_int_distribution<int> dist(min, max);
        return dist(gen);
    }
    
    // Random angle for continuous sampling
    double getRandomAngle(double center, double tolerance, const std::string& method = "uniform");
    
private:
    std::normal_distribution<double> normal_dist{0.0, 1.0};
    std::uniform_real_distribution<double> uniform_real_dist{0.0, 1.0};
};

// ========================
// SAMPLING UTILITIES
// ========================

class SamplingUtils {
    // Check if a conformer (vector of angles) is similar to any in a list, within a given tolerance (deg)
    public:
        // Check if a conformer (vector of angles) is similar to any in a list, within a given tolerance (deg)
        static bool isSimilarConformer(const std::vector<double>& angles1,
                                       const std::vector<std::vector<double>>& visited,
                                       double tolerance);
public:
    // Precompute bonded and same parent pairs for clash detection
    static std::pair<AtomPairs, AtomPairs> precomputeBondedAndSameParentPairs(const RDKit::ROMol& mol);
    
    // Precompute non-bonded pairs for clash detection
    static AtomPairs precomputeNonbondedPairs(const RDKit::ROMol& mol,
                                              const AtomPairs& bonded_pairs,
                                              const AtomPairs& same_parent_pairs);
    
    // Check for clashes between non-bonded atoms
    static bool checkTooCloseNonbondedAtoms(const RDKit::Conformer& conf,
                                            const AtomPairs& candidate_pairs,
                                            double threshold);
    
    // Check timeout condition
    static bool checkTimeout(std::chrono::steady_clock::time_point start_time,
                             int timeout_seconds);
    
    // Extract core angles (remove hetero-H bond angles)
    static std::vector<double> extractCoreAngles(const std::vector<double>& full_angles,
                                                  int num_hetero_H_bonds);
};

// ========================
// ENERGY TRACKER
// ========================

// Simple energy tracker for conformer counting
class ConformerEnergyTracker {
private:
    std::vector<double> energies;
    
public:
    void addEnergy(double energy) {
        energies.push_back(energy);
    }
    
    int getValidConformerCount(double min_energy, double window) const {
        int count = 0;
        for (double energy : energies) {
            if (energy <= min_energy + window) {
                count++;
            }
        }
        return count;
    }
    
    void clear() {
        energies.clear();
    }
};

} // namespace StochasticSampling