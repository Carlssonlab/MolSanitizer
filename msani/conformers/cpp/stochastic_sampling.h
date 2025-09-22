#pragma once

#include "common_types.h"
#include "accelerated_rmsd.h"
#include <string>

namespace StochasticSampling {

/**
 * Discrete Stochastic Sampling (previously stochasticSamplingV2)
 * Uses predefined discrete angle values for torsional sampling
 */
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
                                       int requested_numConfs,
                                       double clash_threshold,
                                       bool verbose = false,
                                       const std::string& mmff_variant = "MMFF94s",
                                       double eps = 1.0,
                                       int randomSeed = 42);

/**
 * Continuous Stochastic Sampling (previously fallback sampling)
 * Uses continuous random angle sampling around torsional peaks
 */
ProductList stochasticSamplingContinuous(RDKit::ROMol& mol,
                                         const ContinuousTorsionMap& torsion_library,
                                         int tolerance_level,
                                         int requested_numConfs,
                                         double window,
                                         int max_attempts,
                                         double rmsd,
                                         double clash_threshold,
                                         const HeteroBonds& hetero_H_bonds,
                                         bool verbose = false,
                                         const std::string& mmff_variant = "MMFF94s",
                                         double eps = 1.0,
                                         const std::string& random_method = "uniform",
                                         int randomSeed = 42);

} // namespace StochasticSampling