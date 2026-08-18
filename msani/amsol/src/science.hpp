#pragma once

#include <cstddef>
#include <span>
#include <vector>

#include <amsolcpp/atom.hpp>
#include <amsolcpp/elements.hpp>
#include <amsolcpp/options.hpp>
#include <amsolcpp/result.hpp>

namespace amsolcpp::detail {

struct IntegralBlock {
    std::size_t atom_i{};
    std::size_t atom_j{};
    std::size_t orbital_count_i{};
    std::size_t orbital_count_j{};
    std::vector<double> values;
};

struct ElectronicSystem {
    std::vector<Atom> atoms;
    std::vector<const ElementParameters*> parameters;
    std::vector<std::size_t> first_orbital;
    std::vector<std::size_t> ao_to_atom;
    std::size_t orbital_count{};
    std::size_t electron_count{};
    std::size_t occupied_orbitals{};
    std::vector<double> hcore;
    std::vector<IntegralBlock> integral_blocks;
    double nuclear_repulsion_ev{};
};

struct ScfState {
    bool converged{};
    bool rescue_used{};
    bool damping_used{};
    std::size_t iterations{};
    double electronic_energy_ev{};
    double energy_residual_kcal{};
    double density_residual{};
    double commutator_residual{};
    std::vector<double> density;
    std::vector<double> orbital_energies;
    std::vector<double> mulliken_charges;
    std::vector<double> cm2_charges;
    std::vector<double> bond_orders;
    std::vector<ScfIterationRecord> records;
};

struct SolvationModel {
    Solvent solvent{Solvent::Water};
    double dielectric{};
    std::vector<double> coulomb_radii;
    std::vector<double> born_radii;
    std::vector<double> fgb;
    std::vector<double> surface_radii;
    std::vector<double> areas;
    std::vector<double> base_sigma;
    std::vector<double> special_sigma;
    double large_surface_contribution{};
};

[[nodiscard]] ElectronicSystem build_electronic_system(
    std::span<const Atom> atoms,
    const CalculationOptions& options
);

[[nodiscard]] ScfState run_scf(
    const ElectronicSystem& system,
    std::span<const double> fgb,
    const CalculationOptions& options
);

[[nodiscard]] SolvationModel build_solvation_model(
    std::span<const Atom> atoms,
    Solvent solvent
);

void complete_cds(
    SolvationModel& model,
    std::span<const Atom> atoms
);

}  // namespace amsolcpp::detail
