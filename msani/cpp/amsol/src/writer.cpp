#include <amsolcpp/writer.hpp>

#include <amsolcpp/elements.hpp>
#include <amsolcpp/exceptions.hpp>

#include "performance_diagnostics.hpp"

#include <array>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace amsolcpp {
namespace {

void json_number(std::ostringstream& output, const double value) {
    if (!std::isfinite(value)) {
        throw NumericalError("cannot serialize a non-finite calculation result");
    }
    output << std::setprecision(17) << value;
}

void json_diagnostic_number(std::ostringstream& output, const double value) {
    if (std::isfinite(value)) {
        output << std::setprecision(17) << value;
    } else {
        output << "null";
    }
}

void json_array(std::ostringstream& output, const std::vector<double>& values) {
    output << '[';
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index != 0) {
            output << ',';
        }
        json_number(output, values[index]);
    }
    output << ']';
}

}  // namespace

std::string result_to_json(
    const std::span<const Atom> atoms,
    const CalculationOptions& options,
    const CalculationResult& result
) {
    AMSOLCPP_PERF_SCOPE(detail::PerformanceStage::OutputFormatting);
    if (result.atomic_charges.size() != atoms.size()
        || result.mulliken_charges.size() != atoms.size()
        || result.atomic_solvation.size() != atoms.size()) {
        throw NumericalError("calculation result atom count mismatch");
    }
    std::ostringstream output;
    output << "{\"converged\":" << (result.converged ? "true" : "false")
           << ",\"scf_iterations\":" << result.scf_iterations
           << ",\"solvent\":\"" << (options.solvent == Solvent::Water ? "water" : "hexadecane") << '\"';
    const std::array<std::pair<const char*, double>, 11> scalars{{
        {"electronic_energy", result.electronic_energy},
        {"nuclear_repulsion_energy", result.nuclear_repulsion_energy},
        {"gas_phase_energy", result.gas_phase_energy},
        {"polarization_free_energy", result.polarization_free_energy},
        {"cds_free_energy", result.cds_free_energy},
        {"solvation_free_energy", result.solvation_free_energy},
        {"solution_phase_energy", result.solution_phase_energy},
        {"total_surface_area", result.total_surface_area},
        {"energy_residual", result.energy_residual},
        {"density_residual", result.density_residual},
        {"commutator_residual", result.commutator_residual},
    }};
    for (const auto& [name, value] : scalars) {
        output << ",\"" << name << "\":";
        json_number(output, value);
    }
    output << ",\"damping_used\":" << (result.damping_used ? "true" : "false")
           << ",\"rescue_used\":" << (result.rescue_used ? "true" : "false")
           << ",\"large_surface_contribution\":";
    json_number(output, result.large_surface_contribution);
    output << ",\"mulliken_charges\":";
    json_array(output, result.mulliken_charges);
    output << ",\"cm2_charges\":";
    json_array(output, result.atomic_charges);
    output << ",\"atomic_charges\":";
    json_array(output, result.atomic_charges);
    output << ",\"orbital_energies\":";
    json_array(output, result.orbital_energies);
    if (!result.iteration_records.empty()) {
        output << ",\"scf_iteration_records\":[";
        for (std::size_t index = 0; index < result.iteration_records.size();
             ++index) {
            if (index != 0) {
                output << ',';
            }
            const auto& record = result.iteration_records[index];
            output << "{\"iteration\":" << record.iteration
                   << ",\"electronic_energy_ev\":";
            json_diagnostic_number(output, record.electronic_energy_ev);
            output << ",\"energy_residual_kcal_mol\":";
            json_diagnostic_number(output, record.energy_residual_kcal_mol);
            output << ",\"density_residual\":";
            json_diagnostic_number(output, record.density_residual);
            output << ",\"commutator_residual\":";
            json_diagnostic_number(output, record.commutator_residual);
            output << ",\"damping_factor\":";
            json_diagnostic_number(output, record.damping_factor);
            output << '}';
        }
        output << ']';
    }
    output << ",\"atoms\":[";
    for (std::size_t index = 0; index < atoms.size(); ++index) {
        if (index != 0) {
            output << ',';
        }
        output << "{\"atomic_number\":" << atoms[index].atomic_number << ",\"symbol\":\""
               << element_symbol(atoms[index].atomic_number) << "\",\"mulliken_charge\":";
        json_number(output, result.mulliken_charges[index]);
        output << ",\"cm2_charge\":";
        json_number(output, result.atomic_charges[index]);
        output << ",\"charge\":";
        json_number(output, result.atomic_charges[index]);
        if (index < result.atomic_solvation.size()) {
            const auto& item = result.atomic_solvation[index];
            output << ",\"polarization_kcal_mol\":";
            json_number(output, item.polarization_kcal_mol);
            output << ",\"area_angstrom2\":";
            json_number(output, item.area_angstrom2);
            output << ",\"surface_tension_cal_mol_angstrom2\":";
            json_number(output, item.surface_tension_cal_mol_angstrom2);
            output << ",\"cds_kcal_mol\":";
            json_number(output, item.cds_kcal_mol);
            output << ",\"subtotal_kcal_mol\":";
            json_number(output, item.subtotal_kcal_mol);
            output << ",\"quadrature_points\":" << item.quadrature_points;
        }
        output << '}';
    }
    output << "]}";
    return output.str();
}

std::string format_legacy_output(
    const std::string_view molecule_name,
    const std::span<const Atom> atoms,
    const CalculationResult& result
) {
    AMSOLCPP_PERF_SCOPE(detail::PerformanceStage::OutputFormatting);
    if (result.atomic_solvation.size() != atoms.size()) {
        throw NumericalError("legacy output requires per-atom solvation results");
    }
    std::ostringstream output;
    output << ' ' << molecule_name << ' ' << atoms.size() << "\n\n"
           << " In the following table subtotal= G_P + SS G_CDS.\n\n"
           << "  Atom   Chem.  CM2      G_P      Area      Sigma k    SS G_CDS Subtotal   M\n"
           << " number symbol   chg.    (kcal) (Ang**2)  cal/(Ang**2)  (kcal)   (kcal)  value\n\n";
    double charge = 0.0;
    for (std::size_t index = 0; index < atoms.size(); ++index) {
        const auto& item = result.atomic_solvation[index];
        output << ' ' << std::setw(3) << index + 1
               << "       " << std::left << std::setw(2) << element_symbol(atoms[index].atomic_number)
               << std::right << std::fixed << std::setprecision(2)
               << std::setw(7) << item.cm2_charge
               << std::setw(9) << item.polarization_kcal_mol << "  "
               << std::setw(7) << item.area_angstrom2 << "    "
               << std::setw(7) << item.surface_tension_cal_mol_angstrom2 << "    "
               << std::setw(7) << item.cds_kcal_mol << "  "
               << std::setw(7) << item.subtotal_kcal_mol << "    "
               << std::setw(3) << item.quadrature_points << '\n';
        charge += item.cm2_charge;
    }
    output << "\n Total:       " << std::fixed << std::setprecision(2)
           << std::setw(7) << charge
           << std::setw(9) << result.polarization_free_energy << "  "
           << std::setw(7) << result.total_surface_area << "               "
           << std::setw(7) << result.cds_free_energy << "  "
           << std::setw(7) << result.solvation_free_energy << '\n';
    return output.str();
}

}  // namespace amsolcpp
