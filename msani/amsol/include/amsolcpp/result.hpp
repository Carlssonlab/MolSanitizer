#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace amsolcpp {

struct AtomicSolvationResult {
    double cm2_charge{};
    double polarization_kcal_mol{};
    double area_angstrom2{};
    double surface_tension_cal_mol_angstrom2{};
    double cds_kcal_mol{};
    double subtotal_kcal_mol{};
    int quadrature_points{};
};

struct ScfIterationRecord {
    std::size_t iteration{};
    double electronic_energy_ev{};
    double energy_residual_kcal_mol{};
    double density_residual{};
    double commutator_residual{};
    double damping_factor{};
};

struct CalculationResult {
    bool converged = false;
    std::size_t scf_iterations = 0;
    int output_decimal_precision = -1;  // mirrors CalculationOptions::output_decimal_precision

    double electronic_energy = 0.0;
    double nuclear_repulsion_energy = 0.0;
    double gas_phase_energy = 0.0;
    double solvation_free_energy = 0.0;
    double solution_phase_energy = 0.0;

    double polarization_free_energy = 0.0;
    double cds_free_energy = 0.0;
    double total_surface_area = 0.0;
    double energy_residual = 0.0;
    double density_residual = 0.0;
    double commutator_residual = 0.0;

    bool damping_used = false;
    bool rescue_used = false;

    // NDDO Mulliken populations and the CM2-corrected charges are both part
    // of the supported scientific result. atomic_charges remains the
    // backwards-compatible CM2 spelling.
    std::vector<double> mulliken_charges;
    std::vector<double> atomic_charges;
    std::vector<double> orbital_energies;
    std::vector<AtomicSolvationResult> atomic_solvation;
    double large_surface_contribution = 0.0;
    std::vector<ScfIterationRecord> iteration_records;
    std::vector<std::string> warnings;
};

struct DualSolventResult {
    CalculationResult water;
    CalculationResult hexadecane;
};

}  // namespace amsolcpp
