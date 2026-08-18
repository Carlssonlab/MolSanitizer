#include <amsolcpp/api.hpp>
#include <amsolcpp/version.hpp>

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#ifdef AMSOLCPP_WITH_RDKIT
// RDKit headers for direct mol object extraction
// (compiled in because AMSOLCPP_WITH_RDKIT is defined)
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/Conformer.h>
#include <Geometry/point.h>
#endif  // AMSOLCPP_WITH_RDKIT

namespace py = pybind11;

// ---------------------------------------------------------------------------
// Helper: extract amsolcpp::Atom list directly from a Python RDKit Mol object.
//
// Only compiled when AMSOLCPP_WITH_RDKIT is defined (i.e. RDKit headers and
// libraries were found at cmake configure time).  When the flag is absent,
// the three public Python functions below are still registered but throw a
// clear RuntimeError so callers get an actionable message instead of an
// obscure ImportError or linker failure.
//
// Accepts any rdkit.Chem.Mol (or RWMol/ROMol) that has already been prepared
// with 3-D coordinates (i.e. it has at least one conformer).  The molecule
// is round-tripped through RDKit's binary pickle format so that we never
// touch the Python-side C++ internals directly – this is safe across all
// RDKit build flavours.
//
// Coordinates are returned in angstroms, exactly as stored by RDKit.
// ---------------------------------------------------------------------------
#ifdef AMSOLCPP_WITH_RDKIT
static std::vector<amsolcpp::Atom> atoms_from_rdkit_mol_impl(py::object mol_obj) {
    if (mol_obj.is_none()) {
        throw std::invalid_argument("rdkit mol object is None");
    }

    // Serialise the Python-side RDKit mol to a binary pickle string.
    py::bytes binary_data;
    try {
        binary_data = mol_obj.attr("ToBinary")();
    } catch (const std::exception& e) {
        throw std::runtime_error(
            std::string("could not call ToBinary() on the rdkit mol object: ") + e.what());
    }
    const std::string binary_str = binary_data.cast<std::string>();

    // Deserialise into an RDKit ROMol using MolPickler.
    RDKit::ROMol mol;
    try {
        RDKit::MolPickler::molFromPickle(binary_str, mol);
    } catch (const std::exception& e) {
        throw std::runtime_error(
            std::string("MolPickler failed to deserialise rdkit mol: ") + e.what());
    }

    // Require at least one 3-D conformer.
    if (mol.getNumConformers() == 0) {
        throw std::invalid_argument(
            "rdkit mol has no conformers; generate 3-D coordinates "
            "(e.g. AllChem.EmbedMolecule) before passing to amsolcpp");
    }

    const RDKit::Conformer& conf = mol.getConformer(0);
    std::vector<amsolcpp::Atom> atoms;
    atoms.reserve(mol.getNumAtoms());

    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
        const RDKit::Atom* atom = mol.getAtomWithIdx(i);
        const RDGeom::Point3D& pos = conf.getAtomPos(i);
        atoms.push_back(amsolcpp::Atom{
            static_cast<int>(atom->getAtomicNum()),
            pos.x,
            pos.y,
            pos.z
        });
    }
    return atoms;
}
#endif  // AMSOLCPP_WITH_RDKIT

// Shared error message used by all three stub registrations below.
static constexpr const char* RDKIT_NOT_COMPILED_MSG =
    "This amsolcpp build was compiled without RDKit support.\n"
    "To enable it, ensure RDKit development headers and libraries are "
    "available in your environment (e.g. install rdkit-dev via conda or "
    "your package manager), then rebuild with:\n"
    "    cmake -DAMSOLCPP_BUILD_PYTHON=ON ...\n"
    "CMake will detect RDKit automatically and set AMSOLCPP_WITH_RDKIT.";

// ---------------------------------------------------------------------------
// SolvDescriptors structs defined at file scope so that compute_solv_descriptors
// (and the RDKit variant) can be referenced before PYBIND11_MODULE.
// ---------------------------------------------------------------------------
struct SolvAtomDescriptor {
  double charge;      // CM2 atomic charge from hexadecane run
  double polar_diff;  // polarization(water) - polarization(hex)
  double surface;     // area_angstrom2 from hexadecane run
  double apolar_diff; // apolar(water) - (apolar(hex) + cs_coeff * area)
  double solv_diff;   // polar_diff + apolar_diff
};

struct SolvDescriptors {
  std::vector<SolvAtomDescriptor> atoms;
  double cs_coeff{};
  double total_surface{};           
  double total_diff_polar{};
  double total_diff_apolar{};
  double total_solv_diff{};
  bool converged_water{};
  bool converged_hexadecane{};
  std::vector<std::string> warnings_water;
  std::vector<std::string> warnings_hexadecane;
};

// ---------------------------------------------------------------------------
// Shared computation: takes a prepared atom list and options, runs dual-solvent
// calculation, and returns SolvDescriptors.  Used by both
// calculate_solv_descriptors (atom-list path) and the RDKit mol-direct path.
// Must be called WITHOUT the GIL held (caller should release before entry).
// ---------------------------------------------------------------------------
static SolvDescriptors compute_solv_descriptors(
    const std::vector<amsolcpp::Atom>& atoms,
    amsolcpp::CalculationOptions options) {
  const auto dual = amsolcpp::calculate_water_and_hexadecane(atoms, options);
  const int prec = options.output_decimal_precision;

  const auto& wat = dual.water;
  const auto& hex = dual.hexadecane;

  const double factor = (prec >= 0) ? std::pow(10.0, prec) : 0.0;
  auto rnd = [prec, factor](double v) -> double {
    if (prec < 0) return v;
    return std::round(v * factor) / factor;
  };

  const std::size_t n = wat.atomic_solvation.size();

  // cs_coeff = (hex total CDS - sum(hex atomic CDS)) / total hex area.
  double sum_apolar_hex = 0.0;
  double total_area = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    sum_apolar_hex += rnd(hex.atomic_solvation[i].cds_kcal_mol);
    total_area += hex.atomic_solvation[i].area_angstrom2;
  }
  total_area = rnd(total_area);

  const double cds_hex_rounded = rnd(hex.cds_free_energy);
  const double cs_coeff =
      (total_area != 0.0) ? (cds_hex_rounded - sum_apolar_hex) / total_area
                          : 0.0;

  SolvDescriptors result;
  result.cs_coeff = cs_coeff;
  result.converged_water = wat.converged;
  result.converged_hexadecane = hex.converged;
  result.warnings_water = wat.warnings;
  result.warnings_hexadecane = hex.warnings;
  result.atoms.reserve(n);

  double total_surface = 0.0;  // sum of per-atom surface areas
  double total_diff_apolar = 0.0; // sum of per-atom apolar differences
  double total_diff_polar = 0.0; // sum of per-atom polar differences
  double total_solv_diff = 0.0; // sum of per-atom solvation differences
  for (std::size_t i = 0; i < n; ++i) {
    SolvAtomDescriptor a;
    const double charge_r    = rnd(hex.atomic_charges[i]);
    const double polar_wat_r = rnd(wat.atomic_solvation[i].polarization_kcal_mol);
    const double polar_hex_r = rnd(hex.atomic_solvation[i].polarization_kcal_mol);
    const double area_r      = rnd(hex.atomic_solvation[i].area_angstrom2);
    const double apolar_wat_r = rnd(wat.atomic_solvation[i].cds_kcal_mol);
    const double apolar_hex_r = rnd(hex.atomic_solvation[i].cds_kcal_mol);

    a.charge      = charge_r;
    a.polar_diff  = polar_wat_r - polar_hex_r;
    a.surface     = area_r;
    a.apolar_diff = apolar_wat_r - (apolar_hex_r + cs_coeff * area_r);
    a.solv_diff   = a.polar_diff + a.apolar_diff;
    result.atoms.push_back(a);
    total_surface += hex.atomic_solvation[i].area_angstrom2;
    total_diff_polar += a.polar_diff;
    total_diff_apolar += a.apolar_diff;
    total_solv_diff += a.solv_diff;
  }
  
  result.total_surface = rnd(total_surface);
  result.total_diff_polar = rnd(total_diff_polar);
  result.total_diff_apolar = rnd(total_diff_apolar);
  result.total_solv_diff = rnd(total_solv_diff);
  return result;
}

PYBIND11_MODULE(_amsolcpp, module) {
  module.doc() = "Deterministic native AM1/CM2/SM5.42R implementation";
  module.attr("__version__") = amsolcpp::version;
  auto error = py::register_exception<amsolcpp::Error>(module, "AmsolError");
  auto input_error = py::register_exception<amsolcpp::InputError>(
      module, "InputError", error.ptr());
  py::register_exception<amsolcpp::UnsupportedCalculationError>(
      module, "UnsupportedCalculationError", input_error.ptr());
  py::register_exception<amsolcpp::UnsupportedElementError>(
      module, "UnsupportedElementError", input_error.ptr());
  auto numerical_error = py::register_exception<amsolcpp::NumericalError>(
      module, "NumericalError", error.ptr());
  py::register_exception<amsolcpp::ScfConvergenceError>(
      module, "ScfConvergenceError", numerical_error.ptr());
  py::register_exception<amsolcpp::DiagonalizationError>(
      module, "DiagonalizationError", numerical_error.ptr());

  py::enum_<amsolcpp::Solvent>(module, "Solvent")
      .value("WATER", amsolcpp::Solvent::Water)
      .value("HEXADECANE", amsolcpp::Solvent::Hexadecane);

  py::class_<amsolcpp::Atom>(module, "Atom")
      .def(py::init<int, double, double, double>(), py::arg("atomic_number"),
           py::arg("x_angstrom"), py::arg("y_angstrom"), py::arg("z_angstrom"))
      .def_readwrite("atomic_number", &amsolcpp::Atom::atomic_number)
      .def_readwrite("x_angstrom", &amsolcpp::Atom::x_angstrom)
      .def_readwrite("y_angstrom", &amsolcpp::Atom::y_angstrom)
      .def_readwrite("z_angstrom", &amsolcpp::Atom::z_angstrom);

  py::class_<amsolcpp::CalculationOptions>(module, "CalculationOptions")
      .def(py::init<>())
      .def_readwrite("solvent", &amsolcpp::CalculationOptions::solvent)
      .def_readwrite("molecular_charge",
                     &amsolcpp::CalculationOptions::molecular_charge)
      .def_readwrite("legacy_scf_iterations",
                     &amsolcpp::CalculationOptions::legacy_scf_iterations)
      .def_readwrite("multiplicity",
                     &amsolcpp::CalculationOptions::multiplicity)
      .def_readwrite("max_scf_iterations",
                     &amsolcpp::CalculationOptions::max_scf_iterations)
      .def_readwrite("energy_tolerance",
                     &amsolcpp::CalculationOptions::energy_tolerance)
      .def_readwrite("density_tolerance",
                     &amsolcpp::CalculationOptions::density_tolerance)
      .def_readwrite("commutator_tolerance",
                     &amsolcpp::CalculationOptions::commutator_tolerance)
      .def_readwrite("enable_scf_rescue",
                     &amsolcpp::CalculationOptions::enable_scf_rescue)
      .def_readwrite(
          "collect_iteration_diagnostics",
          &amsolcpp::CalculationOptions::collect_iteration_diagnostics)
      .def_readwrite("output_decimal_precision",
                     &amsolcpp::CalculationOptions::output_decimal_precision);

  py::class_<amsolcpp::AtomicSolvationResult>(module, "AtomicSolvationResult")
      .def_readonly("cm2_charge", &amsolcpp::AtomicSolvationResult::cm2_charge)
      .def_readonly("polarization_kcal_mol",
                    &amsolcpp::AtomicSolvationResult::polarization_kcal_mol)
      .def_readonly("area_angstrom2",
                    &amsolcpp::AtomicSolvationResult::area_angstrom2)
      .def_readonly(
          "surface_tension_cal_mol_angstrom2",
          &amsolcpp::AtomicSolvationResult::surface_tension_cal_mol_angstrom2)
      .def_readonly("cds_kcal_mol",
                    &amsolcpp::AtomicSolvationResult::cds_kcal_mol)
      .def_readonly("subtotal_kcal_mol",
                    &amsolcpp::AtomicSolvationResult::subtotal_kcal_mol)
      .def_readonly("quadrature_points",
                    &amsolcpp::AtomicSolvationResult::quadrature_points);

  py::class_<amsolcpp::ScfIterationRecord>(module, "ScfIterationRecord")
      .def_readonly("iteration", &amsolcpp::ScfIterationRecord::iteration)
      .def_readonly("electronic_energy_ev",
                    &amsolcpp::ScfIterationRecord::electronic_energy_ev)
      .def_readonly("energy_residual_kcal_mol",
                    &amsolcpp::ScfIterationRecord::energy_residual_kcal_mol)
      .def_readonly("density_residual",
                    &amsolcpp::ScfIterationRecord::density_residual)
      .def_readonly("commutator_residual",
                    &amsolcpp::ScfIterationRecord::commutator_residual)
      .def_readonly("damping_factor",
                    &amsolcpp::ScfIterationRecord::damping_factor);

  // Helper: round v to prec decimal places; returns v unchanged when prec < 0.
  auto dyn_round = [](double v, int prec) -> double {
    if (prec < 0)
      return v;
    const double factor = std::pow(10.0, prec);
    return std::round(v * factor) / factor;
  };

  py::class_<amsolcpp::CalculationResult>(module, "CalculationResult")
      .def(py::init<>())
      .def_readonly("converged", &amsolcpp::CalculationResult::converged)
      .def_readonly("scf_iterations",
                    &amsolcpp::CalculationResult::scf_iterations)
      .def_readonly("output_decimal_precision",
                    &amsolcpp::CalculationResult::output_decimal_precision)
      .def_property_readonly("electronic_energy",
                             [dyn_round](const amsolcpp::CalculationResult &r) {
                               return dyn_round(r.electronic_energy,
                                                r.output_decimal_precision);
                             })
      .def_property_readonly("nuclear_repulsion_energy",
                             [dyn_round](const amsolcpp::CalculationResult &r) {
                               return dyn_round(r.nuclear_repulsion_energy,
                                                r.output_decimal_precision);
                             })
      .def_property_readonly("gas_phase_energy",
                             [dyn_round](const amsolcpp::CalculationResult &r) {
                               return dyn_round(r.gas_phase_energy,
                                                r.output_decimal_precision);
                             })
      .def_property_readonly("solvation_free_energy",
                             [dyn_round](const amsolcpp::CalculationResult &r) {
                               return dyn_round(r.solvation_free_energy,
                                                r.output_decimal_precision);
                             })
      .def_property_readonly("solution_phase_energy",
                             [dyn_round](const amsolcpp::CalculationResult &r) {
                               return dyn_round(r.solution_phase_energy,
                                                r.output_decimal_precision);
                             })
      .def_property_readonly("polarization_free_energy",
                             [dyn_round](const amsolcpp::CalculationResult &r) {
                               return dyn_round(r.polarization_free_energy,
                                                r.output_decimal_precision);
                             })
      .def_property_readonly("cds_free_energy",
                             [dyn_round](const amsolcpp::CalculationResult &r) {
                               return dyn_round(r.cds_free_energy,
                                                r.output_decimal_precision);
                             })
      .def_property_readonly("total_surface_area",
                             [dyn_round](const amsolcpp::CalculationResult &r) {
                               return dyn_round(r.total_surface_area,
                                                r.output_decimal_precision);
                             })
      .def_property_readonly("energy_residual",
                             [dyn_round](const amsolcpp::CalculationResult &r) {
                               return dyn_round(r.energy_residual,
                                                r.output_decimal_precision);
                             })
      .def_property_readonly("density_residual",
                             [dyn_round](const amsolcpp::CalculationResult &r) {
                               return dyn_round(r.density_residual,
                                                r.output_decimal_precision);
                             })
      .def_property_readonly("commutator_residual",
                             [dyn_round](const amsolcpp::CalculationResult &r) {
                               return dyn_round(r.commutator_residual,
                                                r.output_decimal_precision);
                             })
      .def_readonly("damping_used", &amsolcpp::CalculationResult::damping_used)
      .def_readonly("rescue_used", &amsolcpp::CalculationResult::rescue_used)
      .def_property_readonly("mulliken_charges",
                             [dyn_round](const amsolcpp::CalculationResult &r) {
                               if (r.output_decimal_precision < 0)
                                 return r.mulliken_charges;
                               std::vector<double> out;
                               out.reserve(r.mulliken_charges.size());
                               for (double v : r.mulliken_charges)
                                 out.push_back(
                                     dyn_round(v, r.output_decimal_precision));
                               return out;
                             })
      .def_property_readonly("atomic_charges",
                             [dyn_round](const amsolcpp::CalculationResult &r) {
                               if (r.output_decimal_precision < 0)
                                 return r.atomic_charges;
                               std::vector<double> out;
                               out.reserve(r.atomic_charges.size());
                               for (double v : r.atomic_charges)
                                 out.push_back(
                                     dyn_round(v, r.output_decimal_precision));
                               return out;
                             })
      .def_property_readonly("orbital_energies",
                             [dyn_round](const amsolcpp::CalculationResult &r) {
                               if (r.output_decimal_precision < 0)
                                 return r.orbital_energies;
                               std::vector<double> out;
                               out.reserve(r.orbital_energies.size());
                               for (double v : r.orbital_energies)
                                 out.push_back(
                                     dyn_round(v, r.output_decimal_precision));
                               return out;
                             })
      .def_property_readonly(
          "atomic_solvation",
          [dyn_round](const amsolcpp::CalculationResult &r) {
            if (r.output_decimal_precision < 0)
              return r.atomic_solvation;
            std::vector<amsolcpp::AtomicSolvationResult> out;
            out.reserve(r.atomic_solvation.size());
            for (amsolcpp::AtomicSolvationResult a : r.atomic_solvation) {
              // cm2_charge is always returned at full precision
              a.polarization_kcal_mol = dyn_round(a.polarization_kcal_mol,
                                                  r.output_decimal_precision);
              a.area_angstrom2 =
                  dyn_round(a.area_angstrom2, r.output_decimal_precision);
              a.surface_tension_cal_mol_angstrom2 =
                  dyn_round(a.surface_tension_cal_mol_angstrom2,
                            r.output_decimal_precision);
              a.cds_kcal_mol =
                  dyn_round(a.cds_kcal_mol, r.output_decimal_precision);
              a.subtotal_kcal_mol =
                  dyn_round(a.subtotal_kcal_mol, r.output_decimal_precision);
              out.push_back(a);
            }
            return out;
          })
      .def_property_readonly("large_surface_contribution",
                             [dyn_round](const amsolcpp::CalculationResult &r) {
                               return dyn_round(r.large_surface_contribution,
                                                r.output_decimal_precision);
                             })
      .def_property_readonly(
          "iteration_records",
          [dyn_round](const amsolcpp::CalculationResult &r) {
            if (r.output_decimal_precision < 0)
              return r.iteration_records;
            std::vector<amsolcpp::ScfIterationRecord> out;
            out.reserve(r.iteration_records.size());
            for (amsolcpp::ScfIterationRecord rec : r.iteration_records) {
              rec.electronic_energy_ev = dyn_round(rec.electronic_energy_ev,
                                                   r.output_decimal_precision);
              rec.energy_residual_kcal_mol = dyn_round(
                  rec.energy_residual_kcal_mol, r.output_decimal_precision);
              rec.density_residual =
                  dyn_round(rec.density_residual, r.output_decimal_precision);
              rec.commutator_residual = dyn_round(rec.commutator_residual,
                                                  r.output_decimal_precision);
              rec.damping_factor =
                  dyn_round(rec.damping_factor, r.output_decimal_precision);
              out.push_back(rec);
            }
            return out;
          })
      .def_readwrite("warnings", &amsolcpp::CalculationResult::warnings);

  py::class_<amsolcpp::DualSolventResult>(module, "DualSolventResult")
      .def_readonly("water", &amsolcpp::DualSolventResult::water)
      .def_readonly("hexadecane", &amsolcpp::DualSolventResult::hexadecane);

  module.def(
      "calculate",
      [](const std::vector<amsolcpp::Atom> &atoms,
         const amsolcpp::CalculationOptions &options) {
        py::gil_scoped_release release;
        return amsolcpp::calculate(atoms, options);
      },
      py::arg("atoms"), py::arg("options") = amsolcpp::CalculationOptions{});
  module.def(
      "calculate_water_and_hexadecane",
      [](const std::vector<amsolcpp::Atom> &atoms,
         amsolcpp::CalculationOptions options) {
        py::gil_scoped_release release;
        return amsolcpp::calculate_water_and_hexadecane(atoms, options);
      },
      py::arg("atoms"), py::arg("options") = amsolcpp::CalculationOptions{});

  // -----------------------------------------------------------------------
  // SolvDescriptors: structs are defined at file scope (above) so they can
  // also be used by the RDKit direct-mol compute path.
  // -----------------------------------------------------------------------

  py::class_<SolvAtomDescriptor>(module, "SolvAtomDescriptor")
      .def_readonly("charge", &SolvAtomDescriptor::charge)
      .def_readonly("polar_diff", &SolvAtomDescriptor::polar_diff)
      .def_readonly("surface", &SolvAtomDescriptor::surface)
      .def_readonly("apolar_diff", &SolvAtomDescriptor::apolar_diff)
      .def_readonly("solv_diff", &SolvAtomDescriptor::solv_diff);

  // SolvDescriptors: per-atom fields are also exposed as flat C++ vectors
  // (charges, polar_diffs, surfaces, apolar_diffs, solv_diffs) so Python
  // can assign them directly to Solv attributes without list comprehensions.
  py::class_<SolvDescriptors>(module, "SolvDescriptors")
      .def_readonly("atoms", &SolvDescriptors::atoms)
      .def_readonly("cs_coeff", &SolvDescriptors::cs_coeff)
      .def_readonly("converged_water", &SolvDescriptors::converged_water)
      .def_readonly("converged_hexadecane",
                    &SolvDescriptors::converged_hexadecane)
      .def_readonly("warnings_water", &SolvDescriptors::warnings_water)
      .def_readonly("warnings_hexadecane",
                    &SolvDescriptors::warnings_hexadecane)
      .def_readonly("total_surface", &SolvDescriptors::total_surface,
                    "Sum of per-atom rounded surface areas (Angstrom^2).\n"
                    "Computed in C++ as sum(rnd(area_i)) to match the legacy\n"
                    "file-read path where each value was already rounded before\n"
                    "being summed.")
      .def_readonly("total_diff_polar", &SolvDescriptors::total_diff_polar)
      .def_readonly("total_diff_apolar", &SolvDescriptors::total_diff_apolar)
      .def_readonly("total_solv_diff", &SolvDescriptors::total_solv_diff)
      // --- flat per-field array accessors (no Python-side list comprehension) ---
      .def_property_readonly("charges",
          [](const SolvDescriptors& d) {
            std::vector<double> v; v.reserve(d.atoms.size());
            for (const auto& a : d.atoms) v.push_back(a.charge);
            return v;
          }, "CM2 atomic charges from the hexadecane run (List[float]).")
      .def_property_readonly("polar_diffs",
          [](const SolvDescriptors& d) {
            std::vector<double> v; v.reserve(d.atoms.size());
            for (const auto& a : d.atoms) v.push_back(a.polar_diff);
            return v;
          }, "Per-atom polarization free-energy difference water-hex (List[float]).")
      .def_property_readonly("surfaces",
          [](const SolvDescriptors& d) {
            std::vector<double> v; v.reserve(d.atoms.size());
            for (const auto& a : d.atoms) v.push_back(a.surface);
            return v;
          }, "Per-atom solvent-accessible surface area in Angstrom^2 (List[float]).")
      .def_property_readonly("apolar_diffs",
          [](const SolvDescriptors& d) {
            std::vector<double> v; v.reserve(d.atoms.size());
            for (const auto& a : d.atoms) v.push_back(a.apolar_diff);
            return v;
          }, "Per-atom apolar solvation difference water-hex (List[float]).")
      .def_property_readonly("solv_diffs",
          [](const SolvDescriptors& d) {
            std::vector<double> v; v.reserve(d.atoms.size());
            for (const auto& a : d.atoms) v.push_back(a.solv_diff);
            return v;
          }, "Per-atom total solvation difference water-hex (List[float]).");

  module.def(
      "calculate_solv_descriptors",
      [](const std::vector<amsolcpp::Atom>& atoms,
         amsolcpp::CalculationOptions options) -> SolvDescriptors {
        py::gil_scoped_release release;
        return compute_solv_descriptors(atoms, options);
      },
      py::arg("atoms"), py::arg("options") = amsolcpp::CalculationOptions{},
      "Compute water/hexadecane dual-solvent descriptors.\n\n"
      "All per-atom source values are rounded to "
      "options.output_decimal_precision\n"
      "before any arithmetic, matching the rounding semantics of the "
      "file-based .solv reader. Returns a SolvDescriptors object.");
  module.def("calculate_from_input", [](const std::string &input) {
    py::gil_scoped_release release;
    return amsolcpp::calculate_from_input(input);
  });
  module.def("run_legacy_input", [](const std::string &input) {
    py::gil_scoped_release release;
    return amsolcpp::run_legacy_input(input);
  });
  module.def(
      "calculate_batch",
      [](const std::vector<std::vector<amsolcpp::Atom>> &molecules,
         const amsolcpp::CalculationOptions &options,
         const bool continue_on_error) {
        py::gil_scoped_release release;
        return amsolcpp::calculate_batch(molecules, options, continue_on_error);
      },
      py::arg("molecules"), py::arg("options") = amsolcpp::CalculationOptions{},
      py::arg("continue_on_error") = false);

  // -----------------------------------------------------------------------
  // RDKit mol object extraction
  // -----------------------------------------------------------------------
  // These three functions are always registered in the module so that
  // callers always get the same attribute-lookup behaviour regardless of
  // how the extension was built.  When AMSOLCPP_WITH_RDKIT is NOT defined
  // (i.e. RDKit was not found at build time) each function raises a clear
  // RuntimeError that explains what the user must do to enable support.

#ifdef AMSOLCPP_WITH_RDKIT
  // --- Full implementation (RDKit present at build time) ---

  module.def(
      "atoms_from_rdkit_mol",
      [](py::object mol_obj) -> std::vector<amsolcpp::Atom> {
        return atoms_from_rdkit_mol_impl(std::move(mol_obj));
      },
      py::arg("mol"),
      "Extract a list of Atom objects from a Python RDKit Mol object.\n\n"
      "The molecule must already have 3-D coordinates (at least one conformer).\n"
      "Atom positions are taken from the first conformer in angstroms.\n\n"
      "This replaces the Python-side mol_to_amsol_atoms() helper and avoids\n"
      "the overhead of iterating atoms in Python.");

  module.def(
      "calculate_from_rdkit_mol",
      [](py::object mol_obj, amsolcpp::CalculationOptions options)
          -> amsolcpp::CalculationResult {
        auto atoms = atoms_from_rdkit_mol_impl(mol_obj);
        py::gil_scoped_release release;
        return amsolcpp::calculate(atoms, options);
      },
      py::arg("mol"), py::arg("options") = amsolcpp::CalculationOptions{},
      "Run a single-solvent AM1/CM2/SM5.42R calculation directly from an\n"
      "RDKit Mol object.  The molecule must have at least one 3-D conformer.");

  module.def(
      "calculate_solv_descriptors_from_rdkit_mol",
      [](py::object mol_obj, amsolcpp::CalculationOptions options)
          -> SolvDescriptors {
        // Extract atoms from the Python RDKit mol object, then run the full
        // dual-solvent computation.  GIL must be held during mol extraction
        // (ToBinary call), then released for the numerical computation.
        auto atoms = atoms_from_rdkit_mol_impl(mol_obj);

        py::gil_scoped_release release;
        return compute_solv_descriptors(atoms, options);
      },
      py::arg("mol"), py::arg("options") = amsolcpp::CalculationOptions{},
      "Extract atoms from an RDKit Mol object and compute water/hexadecane\n"
      "dual-solvent descriptors in one call.\n\n"
      "The molecule must already have at least one 3-D conformer.\n"
      "Returns a SolvDescriptors object whose per-field array properties\n"
      "(charges, polar_diffs, surfaces, apolar_diffs, solv_diffs) can be\n"
      "assigned directly to Solv attributes without list comprehensions.");

  // Expose a flag so Python code can inspect RDKit support at runtime.
  module.attr("rdkit_support") = true;

#else  // AMSOLCPP_WITH_RDKIT not defined
  // --- Stub implementations (RDKit absent at build time) ---

  module.def(
      "atoms_from_rdkit_mol",
      [](py::object /*mol_obj*/) -> std::vector<amsolcpp::Atom> {
        throw std::runtime_error(RDKIT_NOT_COMPILED_MSG);
        return {};
      },
      py::arg("mol"),
      "[RDKit support not compiled in]\n\n"
      "Raises RuntimeError. Rebuild amsolcpp with RDKit headers available\n"
      "to enable this function.");

  module.def(
      "calculate_from_rdkit_mol",
      [](py::object /*mol_obj*/, amsolcpp::CalculationOptions /*options*/)
          -> amsolcpp::CalculationResult {
        throw std::runtime_error(RDKIT_NOT_COMPILED_MSG);
        return {};
      },
      py::arg("mol"), py::arg("options") = amsolcpp::CalculationOptions{},
      "[RDKit support not compiled in]\n\n"
      "Raises RuntimeError. Rebuild amsolcpp with RDKit headers available\n"
      "to enable this function.");

  module.def(
      "calculate_solv_descriptors_from_rdkit_mol",
      [](py::object /*mol_obj*/, amsolcpp::CalculationOptions /*options*/)
          -> SolvDescriptors {
        throw std::runtime_error(RDKIT_NOT_COMPILED_MSG);
        return {};
      },
      py::arg("mol"), py::arg("options") = amsolcpp::CalculationOptions{},
      "[RDKit support not compiled in]\n\n"
      "Raises RuntimeError. Rebuild amsolcpp with RDKit headers available\n"
      "to enable this function.");

  // Expose a flag so Python code can inspect RDKit support at runtime.
  module.attr("rdkit_support") = false;

#endif  // AMSOLCPP_WITH_RDKIT
}
