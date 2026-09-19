#include "stochastic_sampling.h"
#include "common_types.h"
#include "hydroxyl_sampling.h"
#include "accelerated_rmsd.h"

#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <GraphMol/ForceFieldHelpers/FFConvenience.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/MolOps.h>
#include <ForceField/ForceField.h>
#include <Geometry/point.h>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <cstdio>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <functional>
#include <unordered_set>
#include <utility>

namespace StochasticSampling {

// Simplified hash for better performance (avoids rounding operations)
struct VectorDoubleHash {
    std::size_t operator()(const std::vector<double>& v) const {
        std::size_t seed = v.size();
        for (auto& angle : v) {
            // Direct hash without rounding for better performance
            seed ^= std::hash<double>{}(angle) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

struct VectorDoubleEqual {
    bool operator()(const std::vector<double>& lhs, const std::vector<double>& rhs) const {
        if (lhs.size() != rhs.size()) return false;
        // Direct comparison since angles are discrete values from the angle map
        return lhs == rhs;
    }
};

// Type aliases for the hash sets
using VisitedStateSet = std::unordered_set<std::vector<double>, VectorDoubleHash, VectorDoubleEqual>;

// ========================
// FORCE FIELD CACHE FOR ENERGY CALCULATIONS
// ========================

// Force field cache for efficient energy calculations
class ForceFieldCache {
private:
    mutable std::unique_ptr<ForceFields::ForceField> cached_ff;
    mutable bool ff_valid = false;
    
public:
    ForceFieldCache() = default;
    
    double calcEnergyFast(RDKit::RWMol& working_mol, 
                         RDKit::MMFF::MMFFMolProperties* mmffMolProperties) const {
        // Construct force field only once
        if (!ff_valid || !cached_ff) {
            cached_ff.reset(RDKit::MMFF::constructForceField(working_mol, mmffMolProperties, 100.0));
            ff_valid = true;
            if (!cached_ff) {
                throw std::runtime_error("Failed to construct force field");
            }
        }
        
        // Rebind the force field to the current conformer and refresh its internal state.
        // This mirrors RDKit's ForceFieldHelpers::OptimizeMoleculeConfs workflow.
        auto& conf = working_mol.getConformer(0);
        auto& positions = cached_ff->positions();
        for (unsigned int i = 0; i < working_mol.getNumAtoms(); ++i) {
            positions[i] = &conf.getAtomPos(i);
        }
        cached_ff->initialize();
        
        return cached_ff->calcEnergy();
    }
    
    void invalidate() {
        ff_valid = false;
        cached_ff.reset();
    }
};

// ========================
// DISCRETE STOCHASTIC SAMPLING IMPLEMENTATION
// ========================

ProductList stochasticSamplingDiscrete(RDKit::ROMol& mol, 
                                       const AngleMap& angle_map,
                                       const ScoreMap& score_map,
                                       long long possible_numConfs,
                                       const ImportanceOrder& importance_order,
                                       double window,
                                       int max_attempts,
                                       const HeteroBonds& hetero_H_bonds,
                                       int timeout_conf,
                                       double rmsd,
                                       int numConfs,
                                       double clash_scale,
                                       bool verbose,
                                       const std::string& mmff_variant,
                                       double eps,
                                       int randomSeed) {

    ProductList products;
    RandomGenerator rand_gen(randomSeed);  // Use provided seed
    
    // Initialize the conformer cache for efficient RMSD checking
    ConformerCache conformer_cache;
    conformer_cache.initialize(mol);
    
    // SIMPLIFIED: Build MMFF properties once (topology-dependent, doesn't change)
    std::unique_ptr<RDKit::MMFF::MMFFMolProperties> mmffMolProperties;
    try {
        mmffMolProperties = std::make_unique<RDKit::MMFF::MMFFMolProperties>(mol, mmff_variant);
        if (!mmffMolProperties || !mmffMolProperties->isValid()) {
            return products;
        }
        
        // // Set MMFF variant (MMFF94 or MMFF94s)
        // mmffMolProperties->setMMFFVariant(mmff_variant);
        
        // Set dielectric constant (eps)
        mmffMolProperties->setMMFFDielectricConstant(eps);
        if (verbose) {
            fprintf(stderr, "MMFF properties created with variant: %s, dielectric constant: %f (will be reused for all conformers)\n",
                    mmff_variant.c_str(), eps);
        }
    } catch (const std::exception& e) {
        std::cerr << "❌ Error creating MMFF properties: " << e.what() << std::endl;
        return products;
    }
    
    // Initialize force field cache for efficient energy calculations
    ForceFieldCache ff_cache;
    
    // OPTIMIZATION: Create a single working molecule to avoid deep copies
    RDKit::RWMol working_mol(mol);
    if (working_mol.getNumConformers() == 0) {
        return products;
    }
    
    // Create mapping from original atom indices to heavy atom indices
    const unsigned int numAtoms = working_mol.getNumAtoms();
    std::vector<int> heavy_atom_mapping;
    heavy_atom_mapping.reserve(numAtoms);
    int heavy_idx = 0;
    for (unsigned int i = 0; i < numAtoms; ++i) {
        const RDKit::Atom* atom = mol.getAtomWithIdx(i);
        if (atom->getAtomicNum() != 1) { // Not hydrogen
            heavy_atom_mapping.push_back(heavy_idx++);
        } else {
            heavy_atom_mapping.push_back(-1); // Mark as hydrogen
        }
    }
    
    // Precompute bonded and same parent pairs
    auto [bonded_pairs, same_parent_pairs] = SamplingUtils::precomputeBondedAndSameParentPairs(mol);
    NonbondedClashPairs nonbonded_pairs = SamplingUtils::precomputeNonbondedPairs(mol, bonded_pairs, same_parent_pairs, clash_scale);
    
    // Track visited states and energy statistics
    VisitedStateSet visited_full;
    VisitedStateSet visited_core;
    int num_hetero_H_bonds = hetero_H_bonds.size();
    double min_energy = 1e6;
    
    // Initialize energy tracker for real-time conformer count monitoring
    ConformerEnergyTracker energy_tracker;
    int current_valid_conformer_count = 0;
    
    // Store angle map keys in vector for indexing
    std::vector<int> angle_keys;
    for (const auto& [key, entry] : angle_map) {
        angle_keys.push_back(key);
    }
    
    // Calculate hydroxyl combinations and scaling
    std::uint64_t num_hydroxyl_combinations = 1;
    int f = 1; // scaling factor
    int core_allocation = numConfs;

    if (num_hetero_H_bonds > 0) {

        for (int i = 0; i < num_hetero_H_bonds; ++i) {
            if (i < static_cast<int>(angle_keys.size())) {
                const auto& entry = angle_map.at(angle_keys[angle_keys.size() - 1 - i]);
                num_hydroxyl_combinations = HydroxylCombinationSampler::cappedProduct(
                    num_hydroxyl_combinations, entry.possible_angles.size());
            }
        }
        
        // Calculate scaling factor from hydroxyl combinations
        if (num_hydroxyl_combinations > 30) {
            f = 30;
        } else {
            f = static_cast<int>(std::min<std::uint64_t>(num_hydroxyl_combinations, 3));
        }
        core_allocation = (f < 30) ? core_allocation / f : core_allocation / 30;
        
        if (verbose) {
            fprintf(stderr, "Hydroxyl combinations: %llu, f: %d, core allocation: %d\n",
                    static_cast<unsigned long long>(num_hydroxyl_combinations), f, core_allocation);
        }
        possible_numConfs /= static_cast<long long>(num_hydroxyl_combinations);
        if (verbose) {
            fprintf(stderr, "Adjusted possible_numConfs: %lld\n", possible_numConfs);
        }
    }
    
    // Main sampling logic; avoid over-enumerate for compounds with too many hydroxyls (sugars). For these, use stochastic sampling.
    std::uint64_t full_combination_count = 1;
    for (const auto& [key, entry] : angle_map) {
        full_combination_count = HydroxylCombinationSampler::cappedProduct(
            full_combination_count, entry.possible_angles.size());
    }
    if (possible_numConfs <= max_attempts &&
        full_combination_count <= static_cast<std::uint64_t>(max_attempts) * 5 &&
        num_hydroxyl_combinations <= static_cast<std::uint64_t>(max_attempts) * 5 /
            static_cast<std::uint64_t>(std::max(1LL, possible_numConfs))) {
        // ENUMERATION PATH WITH RANDOMNESS
        if (verbose) {
            fprintf(stderr, "Using enumeration approach with randomness (combinations: %lld)\n", possible_numConfs);
        }
        
        // Pregenerate all combinations
        std::vector<std::vector<double>> all_combinations;
        std::function<void(std::size_t, std::vector<double>&)> generate_combinations;
        generate_combinations = [&](std::size_t depth, std::vector<double>& current_angles) {
            if (depth == angle_keys.size()) {
                all_combinations.push_back(current_angles);
                return;
            }
            
            const auto& entry = angle_map.at(angle_keys[depth]);
            for (double angle : entry.possible_angles) {
                current_angles.push_back(angle);
                generate_combinations(depth + 1, current_angles);
                current_angles.pop_back();
            }
        };
        
        std::vector<double> current_angles;
        generate_combinations(0, current_angles);
        
        // Shuffle combinations for randomness
        rand_gen.shuffle(all_combinations);
        if (verbose) {
            fprintf(stderr, "Generated and shuffled %zu angle combinations\n", all_combinations.size());
        }

        // Process combinations
        int attempts = 0;
        int max_stagnation = std::max(std::min(max_attempts / 10, 3000), 50);
        int stagnation_counter = 0;
        int last_product_size = 0;
        auto start_time = std::chrono::steady_clock::now();
        int timeout_check_counter = 0;
        
        RDKit::Conformer& work_conf = working_mol.getConformer(0);
        
        size_t combination_index = 0;
        while (combination_index < all_combinations.size() && 
               current_valid_conformer_count < core_allocation && 
               attempts < max_attempts) {

            // Check timeout every 10 iterations
            if (timeout_conf > 0 && timeout_check_counter % 10 == 0 && 
                SamplingUtils::checkTimeout(start_time, timeout_conf)) {
                if (verbose) {
                    fprintf(stderr, "Timeout criteria met (attempted %d). Generated %zu conformers (estimated valid within window: %d).\n",
                            attempts, products.size(), current_valid_conformer_count);
                }
                break;
            }
            timeout_check_counter++;
            
            if (stagnation_counter >= max_stagnation) {
                if (verbose) {
                    fprintf(stderr, "Early stopping criteria met (attempted %d). Generated %zu conformers (estimated valid within window: %d).\n",
                            attempts, products.size(), current_valid_conformer_count);
                }
                break;
            }
            const auto& angle_combination = all_combinations[combination_index];
            std::vector<double> core_angles = SamplingUtils::extractCoreAngles(angle_combination, num_hetero_H_bonds);
            
            // Early exit: Check if this core has already been processed
            if (!core_angles.empty() && visited_core.find(core_angles) != visited_core.end()) {
                attempts++;
                combination_index++;
                continue;
            }
            
            // Apply current combination of angles to working molecule
            bool success = true;
            
            for (size_t i = 0; i < angle_keys.size() && success; ++i) {
                const auto& entry = angle_map.at(angle_keys[i]);
                
                MolTransforms::setDihedralDeg(work_conf, 
                                            entry.dihedral_atoms[0],
                                            entry.dihedral_atoms[1], 
                                            entry.dihedral_atoms[2],
                                            entry.dihedral_atoms[3],
                                            angle_combination[i]);
            }

            const bool has_clash = SamplingUtils::checkTooCloseNonbondedAtoms(work_conf, nonbonded_pairs);
            if (has_clash) {
                stagnation_counter++;
                attempts++;
                combination_index++;
                continue;
            }

            if (success) {
                try {
                    double energy = ff_cache.calcEnergyFast(working_mol, const_cast<RDKit::MMFF::MMFFMolProperties*>(mmffMolProperties.get()));
                    
                    if (energy < min_energy) {
                        min_energy = energy;
                    }
                    
                    if (energy <= min_energy + window) {
                        bool should_add = false;
                        
                        auto prepared = conformer_cache.prepareHeavyConformer(work_conf, heavy_atom_mapping);
                        bool is_similar = conformer_cache.isSimilarFast(prepared, rmsd);
                        
                        if (!is_similar) {
                            should_add = true;
                            visited_core.insert(core_angles);
                            conformer_cache.addPreparedConformer(std::move(prepared), energy);
                        }
                        
                        if (should_add) {
                            energy_tracker.addEnergy(energy);
                            
                            if (energy == min_energy) {
                                current_valid_conformer_count = energy_tracker.getValidConformerCount(min_energy, window);
                            } else {
                                current_valid_conformer_count++;
                            }
                            
                            ConformerResult result(work_conf, energy);
                            products.push_back(result);
                        }
                    }
                } catch (const std::exception& e) {
                    std::cerr << "Error calculating energy: " << e.what() << std::endl;
                }
            }
            
            // Check for stagnation
            if (products.size() == static_cast<size_t>(last_product_size)) {
                stagnation_counter++;
            } else {
                stagnation_counter = 0;
                last_product_size = products.size();
            }
            
            attempts++;
            combination_index++;
        }
        
    } else {
        // STOCHASTIC SAMPLING PATH WITH IMPORTANCE SAMPLING
        if (verbose) {
            fprintf(stderr, "Using stochastic sampling approach with importance-based weights (combinations: %lld)\n", possible_numConfs);
        }
        
        // Initialize visiting state tracking
        std::vector<double> visiting(angle_map.size(), 0.0);
        
        // Reset initial angle to the first possible angle for each bond
        int idx = 0;
        RDKit::Conformer& work_conf = working_mol.getConformer(0);
        
        for (const auto& [key, entry] : angle_map) {
            if (!entry.possible_angles.empty()) {
                visiting[idx] = entry.possible_angles[0];
                
                MolTransforms::setDihedralDeg(work_conf,
                                            entry.dihedral_atoms[0],
                                            entry.dihedral_atoms[1], 
                                            entry.dihedral_atoms[2],
                                            entry.dihedral_atoms[3],
                                            visiting[idx]);
            }
            idx++;
        }
        
        // Adaptive sampling parameters
        int max_stagnation = std::min(max_attempts / 10, 2000);
        int stagnation_counter = 0;
        int last_product_size = 0;
        int k = angle_map.size();  // Number of bonds to rotate each iteration
        
        int attempts = 0;
        auto start_time = std::chrono::steady_clock::now();
        int timeout_check_counter = 0;
        
        while (current_valid_conformer_count < core_allocation && attempts < max_attempts) {
            // Check timeout every 10 iterations
            if (timeout_conf > 0 && timeout_check_counter % 10 == 0 && 
                SamplingUtils::checkTimeout(start_time, timeout_conf)) {
                if (verbose) {
                    fprintf(stderr, "Timeout criteria met (attempted %d). Generated %zu conformers (estimated valid within window: %d).\n",
                            attempts, products.size(), current_valid_conformer_count);
                }
                break;
            }
            timeout_check_counter++;
            
            if (stagnation_counter >= max_stagnation) {
                if (verbose) {
                    fprintf(stderr, "Early stopping criteria met (attempted %d). Generated %zu conformers (estimated valid within window: %d).\n",
                            attempts, products.size(), current_valid_conformer_count);
                }
                break;
            }
            // Use importance-based weights for rotation selection
            std::vector<int> to_rotate_raw = rand_gen.weightedChoices(importance_order, k);
            
            // Deduplicate selected bond indices
            std::unordered_set<int> unique_bonds(to_rotate_raw.begin(), to_rotate_raw.end());
            std::vector<int> to_rotate(unique_bonds.begin(), unique_bonds.end());
            
            // For each selected bond, choose a random angle
            for (int bond_idx : to_rotate) {
                if (bond_idx >= static_cast<int>(angle_map.size())) continue;
                
                auto map_it = angle_map.begin();
                std::advance(map_it, bond_idx);
                const auto& entry = map_it->second;
                
                // Get score weights for this bond
                auto score_it = score_map.find(map_it->first);
                std::vector<double> angle_scores;
                if (score_it != score_map.end()) {
                    angle_scores = score_it->second;
                } else {
                    angle_scores.assign(entry.possible_angles.size(), 1.0);
                }
                
                int angle_idx = rand_gen.weightedChoice(angle_scores);
                double angle = entry.possible_angles[angle_idx];
                visiting[bond_idx] = angle;
                
                MolTransforms::setDihedralDeg(work_conf,
                                            entry.dihedral_atoms[0],
                                            entry.dihedral_atoms[1], 
                                            entry.dihedral_atoms[2],
                                            entry.dihedral_atoms[3],
                                            angle);
            }

            // Check if this state has been visited
            std::vector<double> state_tuple = visiting;
            std::vector<double> core_angles = SamplingUtils::extractCoreAngles(state_tuple, num_hetero_H_bonds);
            
            if (visited_full.find(state_tuple) != visited_full.end()) {
                stagnation_counter++;
                attempts++;
                continue;
            }
            visited_full.insert(state_tuple);
            
            // Early exit: Check if this core has already been processed
            if (num_hetero_H_bonds > 0 && (core_angles.empty() || visited_core.find(core_angles) != visited_core.end())) {
                stagnation_counter++;
                attempts++;
                continue;
            }
            // Check for clashes
            const bool has_clash = SamplingUtils::checkTooCloseNonbondedAtoms(work_conf, nonbonded_pairs);
            if (has_clash) {
                stagnation_counter++;
                attempts++;
                continue;
            }
            
            // Calculate energy and apply filtering
            try {
                double energy = ff_cache.calcEnergyFast(working_mol, const_cast<RDKit::MMFF::MMFFMolProperties*>(mmffMolProperties.get()));
                
                if (energy < min_energy) {
                    min_energy = energy;
                }
                
                if (energy <= min_energy + window) {
                    bool should_add = false;
                    
                    auto prepared = conformer_cache.prepareHeavyConformer(work_conf, heavy_atom_mapping);
                    bool is_similar = conformer_cache.isSimilarFast(prepared, rmsd);
                    
                    if (!is_similar) {
                        should_add = true;
                        visited_core.insert(core_angles);
                        conformer_cache.addPreparedConformer(std::move(prepared), energy);
                    }
                    
                    if (should_add) {
                        energy_tracker.addEnergy(energy);
                        
                        if (energy == min_energy) {
                            current_valid_conformer_count = energy_tracker.getValidConformerCount(min_energy, window);
                        } else {
                            current_valid_conformer_count++;
                        }
                        
                        ConformerResult result(work_conf, energy);
                        products.push_back(result);
                    }
                }
                
                // Check for stagnation
                if (products.size() == static_cast<size_t>(last_product_size)) {
                    stagnation_counter++;
                } else {
                    stagnation_counter = 0;
                    last_product_size = products.size();
                }
                
            } catch (const std::exception& e) {
                std::cerr << "Error in force field calculation: " << e.what() << std::endl;
            }
            
            attempts++;
        }
    }
    
    // Sort by energy (lowest first)
    std::sort(products.begin(), products.end(), 
              [](const ConformerResult& a, const ConformerResult& b) {
                  return a.energy < b.energy;
              });
        
    // Apply energy window filtering
    if (!products.empty() && window > 0) {
        double min_energy_final = products[0].energy;
        std::vector<ConformerResult> filtered_products;
        for (const auto& product : products) {
            if (product.energy - min_energy_final <= window) {
                filtered_products.push_back(product);
            }
        }
        products = std::move(filtered_products);
        if (verbose) {
            fprintf(stderr, "Energy window filter (%.2f kcal/mol): %zu conformers within window\n", window, products.size());
        }
    }
    
    // Handle hydroxyl rotations if present
    if (num_hetero_H_bonds > 0 && !products.empty()) {
        if (verbose) {
            fprintf(stderr, "Starting hydroxyl rotation enumeration...\n");
        }
        
        // Resize products to core_allocation size
            if (products.size() > static_cast<size_t>(core_allocation)) {
            products.resize(core_allocation);
            if (verbose) {
                fprintf(stderr, "Resized to core allocation: %d core conformers\n", core_allocation);
            }
        }
        
        std::vector<std::vector<double>> hydroxyl_choices;
        for (int depth = 0; depth < num_hetero_H_bonds; ++depth) {
            const auto& entry = angle_map.at(angle_keys[angle_keys.size() - 1 - depth]);
            hydroxyl_choices.push_back(entry.possible_angles);
        }
        HydroxylCombinationSampler hydroxyl_sampler(std::move(hydroxyl_choices));

        // For each core conformer, apply random f variations of hydroxyl orientations
        std::vector<ConformerResult> final_products;
        final_products.reserve(products.size() * f);
        bool reached_limit = false;
        
        for (const auto& core_product : products) {
            auto shuffled_combinations = hydroxyl_sampler.sample(f + 2, rand_gen);
            
            int variations_to_try = std::min(f + 2, static_cast<int>(shuffled_combinations.size()));
            
            for (int var = 0; var < variations_to_try; ++var) {
                try {
                    RDKit::Conformer& work_conf = working_mol.getConformer(0);
                    // Copy coordinates from core conformer
                    for (unsigned int i = 0; i < core_product.conformer.getNumAtoms(); ++i) {
                        work_conf.setAtomPos(i, core_product.conformer.getAtomPos(i));
                    }
                    
                    // Apply hydroxyl angles only
                    const auto& hydroxyl_angles = shuffled_combinations[var];
                    for (int h = 0; h < num_hetero_H_bonds; ++h) {
                        int bond_idx = angle_keys.size() - 1 - h;
                        const auto& entry = angle_map.at(angle_keys[bond_idx]);
                        MolTransforms::setDihedralDeg(work_conf,
                                                    entry.dihedral_atoms[0],
                                                    entry.dihedral_atoms[1], 
                                                    entry.dihedral_atoms[2],
                                                    entry.dihedral_atoms[3],
                                                    hydroxyl_angles[h]);
                    }

                    // Check for clashes
                    const bool has_clash = SamplingUtils::checkTooCloseNonbondedAtoms(work_conf, nonbonded_pairs);
                    if (has_clash) {
                        continue;
                    }

                    // Calculate energy for this variation
                    double energy = ff_cache.calcEnergyFast(working_mol, const_cast<RDKit::MMFF::MMFFMolProperties*>(mmffMolProperties.get()));

                    if (energy > min_energy + window) {
                        continue;
                    }
                    ConformerResult hydroxyl_result(work_conf, energy);
                    final_products.push_back(hydroxyl_result);
                    
                        if (final_products.size() >= static_cast<size_t>(numConfs)) {
                        if (verbose) {
                            fprintf(stderr, "Early termination: reached numConfs of %d conformers\n", numConfs);
                        }
                        reached_limit = true;
                        break;
                    }
                    
                } catch (const std::exception& e) {
                    std::cerr << "Error in hydroxyl variation: " << e.what() << std::endl;
                    continue;
                }
            }
            if (reached_limit) break;
        }
        
        if (verbose) {
            fprintf(stderr, "Generated %zu conformers with hydroxyl variations\n", final_products.size());
        }
        
        products = std::move(final_products);
    }
    
    // Apply final number filtering to numConfs
        if (products.size() > static_cast<size_t>(numConfs)) {
        products.resize(numConfs);
        if (verbose) {
            fprintf(stderr, "Final count filter: Limited to %d conformers\n", numConfs);
        }
    }
    
    return products;
}


// ========================
// CONTINUOUS STOCHASTIC SAMPLING (FALLBACK/ANGLE-BASED) IMPLEMENTATION
// ========================

ProductList stochasticSamplingContinuous(RDKit::ROMol& mol,
                                         const ContinuousTorsionMap& torsion_library,
                                         int tolerance_level,
                                         int numConfs,
                                         double window,
                                         int max_attempts,
                                         int timeout_conf,
                                         double rmsd,
                                         double clash_scale,
                                         const HeteroBonds& hetero_H_bonds,
                                         bool verbose,
                                         const std::string& mmff_variant,
                                         double eps,
                                         const std::string& random_method,
                                         int randomSeed) {
    ProductList products;
    RandomGenerator rand_gen(randomSeed);  // Use provided seed

    // Build MMFF properties once
    std::unique_ptr<RDKit::MMFF::MMFFMolProperties> mmffMolProperties;
    try {
        mmffMolProperties = std::make_unique<RDKit::MMFF::MMFFMolProperties>(mol);
        if (!mmffMolProperties || !mmffMolProperties->isValid()) {
            return products;
        }
        mmffMolProperties->setMMFFVariant(mmff_variant);
        mmffMolProperties->setMMFFDielectricConstant(eps);
        if (verbose) {
            fprintf(stderr, "[Fallback] MMFF properties created with variant: %s, dielectric constant: %f\n",
                    mmff_variant.c_str(), eps);
        }
    } catch (const std::exception& e) {
        std::cerr << "❌ Error creating MMFF properties: " << e.what() << std::endl;
        return products;
    }

    // Initialize force field cache
    ForceFieldCache ff_cache;

    // Create working molecule
    RDKit::RWMol working_mol(mol);
    if (working_mol.getNumConformers() == 0) {
        return products;
    }
    
    // Prepare for RMSD-based deduplication
    // Create heavy atom mapping
    const unsigned int numAtoms = working_mol.getNumAtoms();
    std::vector<int> heavy_atom_mapping;
    heavy_atom_mapping.reserve(numAtoms);
    int heavy_idx = 0;
    for (unsigned int i = 0; i < numAtoms; ++i) {
        const RDKit::Atom* atom = mol.getAtomWithIdx(i);
        if (atom->getAtomicNum() != 1) {
            heavy_atom_mapping.push_back(heavy_idx++);
        } else {
            heavy_atom_mapping.push_back(-1);
        }
    }
    ConformerCache conformer_cache;
    conformer_cache.initialize(mol);

    // Precompute pairs for clash detection
    auto [bonded_pairs, same_parent_pairs] = SamplingUtils::precomputeBondedAndSameParentPairs(mol);
    NonbondedClashPairs nonbonded_pairs = SamplingUtils::precomputeNonbondedPairs(mol, bonded_pairs, same_parent_pairs, clash_scale);

    

    RDKit::Conformer& work_conf = working_mol.getConformer(0);

    // Main sampling loop
    std::vector<std::vector<double>> visited_angles;
    double min_energy = 1e6;
    int attempts = 0;
    int stagnation_counter = 0;
    int last_product_size = 0;
    int max_stagnation = std::min(max_attempts / 10, 5000);
    int n_transform = torsion_library.size();
    std::vector<double> visiting(n_transform, 0.0);
    // Add energy tracker and conformer count as in discrete
    ConformerEnergyTracker energy_tracker;
    int current_valid_conformer_count = 0;
    int num_hetero_H_bonds = hetero_H_bonds.size();
    int f = 1;
    int core_allocation = numConfs;
    std::uint64_t num_hydroxyl_combinations = 1;
    auto start_time = std::chrono::steady_clock::now();
    int timeout_check_counter = 0;
    

    // Initialize visiting state with first angles from each bond
    int idx = 0;
    for (const auto& [bond_id, bond_info] : torsion_library) {
        if (!bond_info.peaks.empty()) {
            double angle = rand_gen.getRandomAngle(
                bond_info.peaks[0].center,
                bond_info.peaks[0].tolerance[tolerance_level-1],
                random_method
            );
            visiting[idx] = angle;
            MolTransforms::setDihedralDeg(work_conf,
                                         bond_info.dihedral_atoms[0],
                                         bond_info.dihedral_atoms[1],
                                         bond_info.dihedral_atoms[2],
                                         bond_info.dihedral_atoms[3],
                                         angle);
        }
        idx++;
    }




    // Build bond_indices once outside the loop
    std::vector<int> bond_indices;
    std::vector<std::vector<double>> bond_weights;
    for (const auto& [bond_id, bond_info] : torsion_library) {
        bond_indices.push_back(bond_id);
        std::vector<double> weights;
        weights.reserve(bond_info.peaks.size());
        for (const auto& peak : bond_info.peaks) weights.push_back(peak.weight);
        bond_weights.push_back(std::move(weights));
    }

    if (num_hetero_H_bonds > 0) {
        for (int i = 0; i < num_hetero_H_bonds; ++i) {
            // Assuming hetero_H_bonds correspond to the last entries in torsion_library
            if (static_cast<std::size_t>(i) < bond_indices.size()) {
                const std::size_t bond_idx_in_torsion_library =
                    bond_indices.size() - 1 - static_cast<std::size_t>(i);
                int actual_bond_id = bond_indices[bond_idx_in_torsion_library];
                const auto& bond_info = torsion_library.at(actual_bond_id);
                num_hydroxyl_combinations = HydroxylCombinationSampler::cappedProduct(
                    num_hydroxyl_combinations, bond_info.peaks.size());
            }
        }
        
        // Calculate scaling factor from hydroxyl combinations
        if (num_hydroxyl_combinations > 30) {
            f = 30;
        } else {
            f = static_cast<int>(std::min<std::uint64_t>(num_hydroxyl_combinations, 3));
        }
        core_allocation = (f < 30) ? core_allocation / f : core_allocation / 30;
        if (verbose) {
            fprintf(stderr, "[Fallback] Hydroxyl combinations: %llu, f: %d, core allocation: %d\n",
                    static_cast<unsigned long long>(num_hydroxyl_combinations), f, core_allocation);
        }
    }

    // Only check isSimilarConformer to the "core" (i.e., angles except hetero_H_bonds)
    std::vector<std::vector<double>> visited_core_angles;
    while (current_valid_conformer_count < core_allocation && attempts < max_attempts) {
        // Timeout check (every 10 iterations)
        if (timeout_conf > 0 && timeout_check_counter % 10 == 0 && 
            SamplingUtils::checkTimeout(start_time, timeout_conf)) {
            if (verbose) {
                fprintf(stderr, "Timeout criteria met (attempted %d). Generated %zu conformers (estimated valid within window: %d).\n",
                        attempts, products.size(), current_valid_conformer_count);
            }
            break;
        }
        timeout_check_counter++;
        // Rotate n_transform bonds randomly
        for (int i = 0; i < n_transform; ++i) {
            int bond_idx = rand_gen.randint(0, n_transform - 1);
            int actual_bond_id = bond_indices[bond_idx];
            const auto& bond_info = torsion_library.at(actual_bond_id);
            int peak_idx = rand_gen.weightedChoice(bond_weights[bond_idx]);
            if (peak_idx < 0 || peak_idx >= (int)bond_info.peaks.size()) peak_idx = 0;
            double angle = rand_gen.getRandomAngle(
                bond_info.peaks[peak_idx].center,
                bond_info.peaks[peak_idx].tolerance[tolerance_level-1],
                random_method
            );
            visiting[bond_idx] = angle;
            MolTransforms::setDihedralDeg(work_conf,
                                         bond_info.dihedral_atoms[0],
                                         bond_info.dihedral_atoms[1],
                                         bond_info.dihedral_atoms[2],
                                         bond_info.dihedral_atoms[3],
                                         angle);
        }


        // Extract core angles (excluding hetero_H_bonds)
        std::vector<double> core_angles = SamplingUtils::extractCoreAngles(visiting, num_hetero_H_bonds);
        bool is_similar_core = SamplingUtils::isSimilarConformer(core_angles, visited_core_angles, 30.0);
        if (is_similar_core) {
            attempts++;
            stagnation_counter++;
            continue;
        }
        visited_core_angles.push_back(core_angles);

        // Check for clashes
        bool has_clash = SamplingUtils::checkTooCloseNonbondedAtoms(work_conf, nonbonded_pairs);
        if (has_clash) {
            attempts++;
            stagnation_counter++;
            continue;
        }

        // Calculate energy
        try {
            double energy = ff_cache.calcEnergyFast(working_mol, const_cast<RDKit::MMFF::MMFFMolProperties*>(mmffMolProperties.get()));
            if (energy < min_energy) {
                min_energy = energy;
            }
            if (energy <= min_energy + window) {
                // Only check RMSD similarity for core conformers
                auto prepared = conformer_cache.prepareHeavyConformer(work_conf, heavy_atom_mapping);
                bool is_similar_rmsd = conformer_cache.isSimilarFast(prepared, rmsd);
                if (!is_similar_rmsd) {
                    energy_tracker.addEnergy(energy);
                    if (energy == min_energy) {
                        current_valid_conformer_count = energy_tracker.getValidConformerCount(min_energy, window);
                    } else {
                        current_valid_conformer_count++;
                    }
                    ConformerResult result(work_conf, energy);
                    products.push_back(result);
                    conformer_cache.addPreparedConformer(std::move(prepared), energy);
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "Error calculating energy: " << e.what() << std::endl;
        }

        // Check for stagnation
        if ((int)products.size() == last_product_size) {
            stagnation_counter++;
        } else {
            stagnation_counter = 0;
            last_product_size = products.size();
        }
        attempts++;
        if (stagnation_counter > max_stagnation) {
            if (verbose) {
                fprintf(stderr, "Early stopping criteria met (attempted %d). Generated %zu conformers.\n",
                        attempts, products.size());
            }
            break;
        }
    }

    // Downsample by f, sample f conformers per core as in Discrete method (to be implemented after this main loop)
    // ...existing code for sorting, filtering, and hydroxyl sampling as in Discrete method should be added here...
    std::sort(products.begin(), products.end(),
              [](const ConformerResult& a, const ConformerResult& b) {
                  return a.energy < b.energy;
              });
    
    // Handle hydroxyl rotations if present
    if (num_hetero_H_bonds > 0 && !products.empty()) {
        if (verbose) {
            fprintf(stderr, "[Fallback] Starting hydroxyl rotation enumeration...\n");
        }
        
        // Resize products to core_allocation size
        if (products.size() > static_cast<size_t>(core_allocation)) {
            products.resize(core_allocation);
            if (verbose) {
                fprintf(stderr, "[Fallback] Resized to core allocation: %d core conformers\n", core_allocation);
            }
        }
        
        std::vector<std::vector<double>> hydroxyl_choices;
        for (int depth = 0; depth < num_hetero_H_bonds; ++depth) {
            const auto& bond = torsion_library.at(bond_indices[bond_indices.size() - 1 - depth]);
            std::vector<double> choices;
            choices.reserve(bond.peaks.size());
            for (const auto& peak : bond.peaks) choices.push_back(peak.center);
            hydroxyl_choices.push_back(std::move(choices));
        }
        HydroxylCombinationSampler hydroxyl_sampler(std::move(hydroxyl_choices));

        // For each core conformer, apply random f variations of hydroxyl orientations
        std::vector<ConformerResult> final_products;
        final_products.reserve(products.size() * f);
        bool reached_limit = false;
        
        for (const auto& core_product : products) {
            auto shuffled_combinations = hydroxyl_sampler.sample(f + 2, rand_gen);
            
            int variations_to_try = std::min(f + 2, static_cast<int>(shuffled_combinations.size()));
            
            for (int var = 0; var < variations_to_try; ++var) {
                try {
                    RDKit::Conformer& work_conf_hydroxyl = working_mol.getConformer(0);
                    
                    // Copy coordinates from core conformer
                    for (unsigned int i = 0; i < core_product.conformer.getNumAtoms(); ++i) {
                        work_conf_hydroxyl.setAtomPos(i, core_product.conformer.getAtomPos(i));
                    }
                    
                    // Apply hydroxyl angles only
                    const auto& hydroxyl_angles = shuffled_combinations[var];
                    for (int h = 0; h < num_hetero_H_bonds; ++h) {
                        int bond_idx_in_torsion_library = bond_indices.size() - 1 - h;
                        int actual_bond_id = bond_indices[bond_idx_in_torsion_library];
                        const auto& bond_info = torsion_library.at(actual_bond_id);
                        MolTransforms::setDihedralDeg(work_conf_hydroxyl,
                                                    bond_info.dihedral_atoms[0],
                                                    bond_info.dihedral_atoms[1],
                                                    bond_info.dihedral_atoms[2],
                                                    bond_info.dihedral_atoms[3],
                                                    hydroxyl_angles[h]);
                    }
                    
                    // Check for clashes
                    if (SamplingUtils::checkTooCloseNonbondedAtoms(work_conf_hydroxyl, nonbonded_pairs)) {
                        continue;
                    }
                    
                    // Calculate energy for this variation
                    double energy = ff_cache.calcEnergyFast(working_mol, const_cast<RDKit::MMFF::MMFFMolProperties*>(mmffMolProperties.get()));
                    
                    if (energy > min_energy + window) {
                        continue;
                    }
                    
                    ConformerResult hydroxyl_result(work_conf_hydroxyl, energy);
                    final_products.push_back(hydroxyl_result);
                    
                    if (final_products.size() >= static_cast<size_t>(numConfs)) {
                        if (verbose) {
                            fprintf(stderr, "[Fallback] Early termination: reached numConfs of %d conformers\n", numConfs);
                        }
                        reached_limit = true;
                        break;
                    }
                    
                } catch (const std::exception& e) {
                    std::cerr << "[Fallback] Error in hydroxyl variation: " << e.what() << std::endl;
                    continue;
                }
            }
            if (reached_limit) break;
        }
        
        if (verbose) {
            fprintf(stderr, "[Fallback] Generated %zu conformers with hydroxyl variations\n", final_products.size());
        }
        
        products = std::move(final_products);
    }
    

    if (products.size() > static_cast<size_t>(numConfs)) {
        products.resize(numConfs);
    }
    if (verbose) {
        fprintf(stderr, "[Fallback] Random angle sampling completed: %zu conformers generated\n", products.size());
    }
    return products;
}

} // namespace StochasticSampling
