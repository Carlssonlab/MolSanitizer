#include <amsolcpp/api.hpp>

#include <amsolcpp/parser.hpp>
#include <amsolcpp/writer.hpp>

#include "science.hpp"
#include "performance_diagnostics.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <string>
#include <utility>

namespace amsolcpp {

namespace {

void validate_options(const CalculationOptions& options) {
    if (options.solvent != Solvent::Water && options.solvent != Solvent::Hexadecane) {
        throw InputError("solvent must be water or hexadecane");
    }
    if (options.max_scf_iterations == 0) throw InputError("max_scf_iterations must be positive");
    if (options.legacy_scf_iterations == 0) throw InputError("legacy_scf_iterations must be positive");
    if (!std::isfinite(options.energy_tolerance) || options.energy_tolerance <= 0.0 ||
        !std::isfinite(options.density_tolerance) || options.density_tolerance <= 0.0 ||
        !std::isfinite(options.commutator_tolerance) || options.commutator_tolerance <= 0.0) {
        throw InputError("SCF tolerances must be finite and positive");
    }
}

CalculationResult assemble_result(std::span<const Atom> atoms, const detail::ElectronicSystem& system,
                                  const detail::SolvationModel& model, const detail::ScfState& scf,
                                  int output_decimal_precision) {
    AMSOLCPP_PERF_SCOPE(detail::PerformanceStage::ResultAssembly);
    const std::size_t n = atoms.size();
    CalculationResult result;
    result.output_decimal_precision = output_decimal_precision;
    result.converged = scf.converged;
    result.scf_iterations = scf.iterations;
    result.electronic_energy = scf.electronic_energy_ev;
    result.nuclear_repulsion_energy = system.nuclear_repulsion_ev;
    result.mulliken_charges = scf.mulliken_charges;
    result.atomic_charges = scf.cm2_charges;
    // Canonicalize values that differ from exact zero by no more than one
    // binary64 epsilon.  This is output normalization only: the SCF and its
    // consistency checks use the unmodified full-precision charges above.
    const auto canonicalize_zero = [](double& value) {
        if (std::abs(value) <= std::numeric_limits<double>::epsilon()) {
            value = 0.0;
        }
    };
    std::for_each(result.mulliken_charges.begin(),
                  result.mulliken_charges.end(), canonicalize_zero);
    std::for_each(result.atomic_charges.begin(),
                  result.atomic_charges.end(), canonicalize_zero);
    result.orbital_energies = scf.orbital_energies;
    result.energy_residual = scf.energy_residual_kcal;
    result.density_residual = scf.density_residual;
    result.commutator_residual = scf.commutator_residual;
    result.damping_used = scf.damping_used;
    result.rescue_used = scf.rescue_used;
    result.iteration_records = scf.records;
    result.atomic_solvation.reserve(n);
    std::vector<double> bp(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            bp[i] += model.fgb[i * n + j] * result.atomic_charges[j];
        }
    }
    const int quadrature_points = std::clamp(static_cast<int>(10.0 + 0.2 * static_cast<double>(n)), 4, 16);
    for (std::size_t i = 0; i < n; ++i) {
        AtomicSolvationResult atom;
        atom.cm2_charge = result.atomic_charges[i];
        atom.polarization_kcal_mol =
            -0.5 * result.atomic_charges[i] * bp[i] * 23.061;
        atom.area_angstrom2 = model.areas[i];
        atom.surface_tension_cal_mol_angstrom2 = model.base_sigma[i] + model.special_sigma[i];
        atom.cds_kcal_mol = 1.0e-3 * atom.surface_tension_cal_mol_angstrom2 * atom.area_angstrom2;
        atom.subtotal_kcal_mol = atom.polarization_kcal_mol + atom.cds_kcal_mol;
        atom.quadrature_points = quadrature_points;
        result.polarization_free_energy += atom.polarization_kcal_mol;
        result.cds_free_energy += atom.cds_kcal_mol;
        result.total_surface_area += atom.area_angstrom2;
        result.atomic_solvation.push_back(atom);
    }
    result.cds_free_energy += model.large_surface_contribution;
    result.large_surface_contribution = model.large_surface_contribution;
    result.solvation_free_energy = result.polarization_free_energy + result.cds_free_energy;
    const double gas_electronic_ev = scf.electronic_energy_ev - result.polarization_free_energy / 23.061;
    result.gas_phase_energy = gas_electronic_ev + system.nuclear_repulsion_ev;
    result.solution_phase_energy = result.gas_phase_energy + result.solvation_free_energy / 23.061;
    return result;
}

}  // namespace

CalculationResult calculate(std::span<const Atom> atoms, const CalculationOptions& options) {
    validate_options(options);
    const auto system = detail::build_electronic_system(atoms, options);
    const auto model = detail::build_solvation_model(atoms, options.solvent);
    const auto scf = detail::run_scf(system, model.fgb, options);
    return assemble_result(atoms, system, model, scf, options.output_decimal_precision);
}

DualSolventResult calculate_water_and_hexadecane(std::span<const Atom> atoms, CalculationOptions options) {
    validate_options(options);
    const auto system = detail::build_electronic_system(atoms, options);
    options.solvent = Solvent::Water;
    const auto water_model = detail::build_solvation_model(atoms, options.solvent);
    const auto water_scf = detail::run_scf(system, water_model.fgb, options);
    auto water = assemble_result(atoms, system, water_model, water_scf, options.output_decimal_precision);
    options.solvent = Solvent::Hexadecane;
    const auto hex_model = detail::build_solvation_model(atoms, options.solvent);
    const auto hex_scf = detail::run_scf(system, hex_model.fgb, options);
    auto hexadecane = assemble_result(atoms, system, hex_model, hex_scf, options.output_decimal_precision);
    return {std::move(water), std::move(hexadecane)};
}

CalculationResult calculate_from_input(const std::string_view input) {
    const auto parsed = parse_legacy_input(input);
    return calculate(parsed.atoms, parsed.options);
}

std::string run_legacy_input(const std::string_view input) {
    const auto parsed = parse_legacy_input(input);
    const auto result = calculate(parsed.atoms, parsed.options);
    return format_legacy_output(parsed.molecule_name, parsed.atoms, result);
}

std::vector<CalculationResult> calculate_batch(
    const std::vector<std::vector<Atom>>& molecules,
    const CalculationOptions& options,
    const bool continue_on_error
) {
    std::vector<CalculationResult> results;
    results.reserve(molecules.size());
    for (std::size_t index = 0; index < molecules.size(); ++index) {
        try {
            results.push_back(calculate(molecules[index], options));
        } catch (const std::exception& error) {
            const std::string message =
                "batch molecule " + std::to_string(index) + " failed: " + error.what();
            if (!continue_on_error) throw InputError(message);
            CalculationResult failed;
            failed.warnings.push_back(message);
            results.push_back(std::move(failed));
        }
    }
    return results;
}

}  // namespace amsolcpp
