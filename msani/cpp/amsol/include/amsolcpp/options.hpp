#pragma once

#include <cstddef>

namespace amsolcpp {

enum class Solvent {
    Water,
    Hexadecane,
};

struct CalculationOptions {
    Solvent solvent = Solvent::Water;
    int molecular_charge = 0;
    int multiplicity = 1;
    std::size_t max_scf_iterations = 300;
    std::size_t legacy_scf_iterations = 200;
    // AMSOL 7.1 ITER defaults SCFCRT to 1e-6 eV and derives PLTEST from
    // SELCON=23.061*SCFCRT.  Keep the public defaults source-compatible so
    // the accepted electronic root matches the legacy scientific model.
    double energy_tolerance = 1.0e-6;
    // The final consistency residual covers the full dense density, whereas
    // PLTEST covers diagonal packed entries.  Four PLTEST units retain the
    // reference stopping envelope for off-diagonal terms.
    double density_tolerance = 1.5265156583998784e-4;
    // ITER does not gate on a commutator.  This independently calibrated
    // guard rejects grossly non-self-consistent states while accepting the
    // residuals of source-converged AMSOL/pyAMSOL trajectories.
    double commutator_tolerance = 1.0e-3;
    bool enable_scf_rescue = true;
    bool collect_iteration_diagnostics = false;
    // Number of decimal places to round output values (-1 = full precision).
    // cm2_charge is always returned at full precision regardless of this setting.
    int output_decimal_precision = -1;
};

}  // namespace amsolcpp
