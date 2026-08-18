#include "science.hpp"
#include "amsol_eigensolver.hpp"
#include "performance_diagnostics.hpp"

#include <amsolcpp/exceptions.hpp>
#include <amsolcpp/packed_matrix.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

namespace amsolcpp::detail {
namespace {

constexpr double ev_to_kcal = 23.061;
constexpr std::array<std::pair<int, int>, 10> product_pairs{{
    {0, 0}, {1, 0}, {1, 1}, {2, 0}, {2, 1},
    {2, 2}, {3, 0}, {3, 1}, {3, 2}, {3, 3}}};
constexpr std::array<std::array<std::size_t, 4>, 4> orbital_pair_indices{{
    {{0, 1, 3, 6}},
    {{1, 2, 4, 7}},
    {{3, 4, 5, 8}},
    {{6, 7, 8, 9}},
}};
#ifdef _OPENMP
// Phase 5 measured whole-calculation benefit at N=97 and N=109.  Keep
// smaller matrices in a serialized region to avoid indiscriminate threading.
constexpr std::size_t openmp_commutator_min_dimension = 96;
#endif

struct Cm2Entry { int a; int b; double value; };
constexpr std::size_t cm2_table_extent = 54;
using Cm2Table = std::array<
    double, cm2_table_extent * cm2_table_extent>;

template <std::size_t Size>
constexpr Cm2Table make_cm2_table(
    const std::array<Cm2Entry, Size>& entries
) {
    Cm2Table result{};
    for (const auto& entry : entries) {
        result[static_cast<std::size_t>(entry.a) * cm2_table_extent
               + static_cast<std::size_t>(entry.b)] = entry.value;
        result[static_cast<std::size_t>(entry.b) * cm2_table_extent
               + static_cast<std::size_t>(entry.a)] = -entry.value;
    }
    return result;
}

constexpr std::array<Cm2Entry, 10> cm2_ck_entries{{
    {1,6,-.02},{1,7,.207},{1,8,.177},{6,7,.008},{7,8,-.197},
    {6,8,.026},{1,14,-.083},{1,16,.038},{6,14,.062},{6,16,-.059}}};
constexpr std::array<Cm2Entry, 14> cm2_dd_entries{{
    {6,7,.086},{6,8,.016},{6,9,.019},{6,17,.027},{6,35,.081},
    {6,53,.147},{6,16,.171},{1,15,.103},{6,15,-.019},{8,15,.088},
    {9,15,.252},{7,8,.134},{8,16,0.0},{15,16,-.080}}};
constexpr Cm2Table cm2_ck_table = make_cm2_table(cm2_ck_entries);
constexpr Cm2Table cm2_dd_table = make_cm2_table(cm2_dd_entries);

double cm2_parameter(const int zi, const int zj, const bool ck) {
    if (zi < 0 || zj < 0
        || static_cast<std::size_t>(zi) >= cm2_table_extent
        || static_cast<std::size_t>(zj) >= cm2_table_extent) {
        return 0.0;
    }
    const std::size_t index = static_cast<std::size_t>(zi)
        * cm2_table_extent + static_cast<std::size_t>(zj);
    return ck ? cm2_ck_table[index] : cm2_dd_table[index];
}

std::vector<double> initial_density(const ElectronicSystem& system, const int charge) {
    AMSOLCPP_PERF_SCOPE(PerformanceStage::InitialDensity);
    const std::size_t n = system.orbital_count;
    std::vector<double> result(n * n, 0.0);
    const double correction = static_cast<double>(charge) / static_cast<double>(n);
    for (std::size_t atom = 0; atom < system.atoms.size(); ++atom) {
        const auto& p = *system.parameters[atom];
        const double fraction = p.atomic_number == 1 ? 1.0 : 0.25;
        const double population = static_cast<double>(p.core_charge) * fraction - correction;
        for (int local = 0; local < p.basis_orbitals; ++local) {
            const std::size_t ao = system.first_orbital[atom] + static_cast<std::size_t>(local);
            result[ao * n + ao] = population;
        }
    }
    return result;
}

void gas_fock(const ElectronicSystem& system, std::span<const double> density,
              std::vector<double>& fock,
              std::vector<double>& weighted_pair_density) {
    AMSOLCPP_PERF_SCOPE(PerformanceStage::GasFock);
    AMSOLCPP_PERF_INCREMENT(gas_fock_constructions);
    const std::size_t n = system.orbital_count;
    fock.assign(system.hcore.begin(), system.hcore.end());
    for (std::size_t atom = 0; atom < system.atoms.size(); ++atom) {
        const auto& p = *system.parameters[atom];
        const std::size_t s = system.first_orbital[atom];
        fock[s * n + s] += 0.5 * density[s * n + s] * p.gss;
        if (p.basis_orbitals == 1) continue;
        double ptpop = 0.0;
        for (std::size_t q = 1; q < 4; ++q) ptpop += density[(s + q) * n + s + q];
        fock[s * n + s] += ptpop * p.gsp - 0.5 * ptpop * p.hsp;
        for (std::size_t q = 1; q < 4; ++q) {
            const std::size_t po = s + q;
            const double pss = density[s * n + s], ppp = density[po * n + po];
            fock[po * n + po] += pss * p.gsp - .5 * pss * p.hsp + .5 * ppp * p.gpp +
                (ptpop - ppp) * p.gp2 - .25 * (ptpop - ppp) * (p.gpp - p.gp2);
            const double value = 2.0 * density[po * n + s] * p.hsp -
                0.5 * density[po * n + s] * (p.hsp + p.gsp);
            fock[po * n + s] += value; fock[s * n + po] += value;
        }
        for (std::size_t high = 2; high < 4; ++high) for (std::size_t low = 1; low < high; ++low) {
            const std::size_t ph = s + high, pl = s + low;
            const double value = density[ph * n + pl] * (p.gpp - p.gp2) -
                .25 * density[ph * n + pl] * (p.gpp + p.gp2);
            fock[ph * n + pl] += value; fock[pl * n + ph] += value;
        }
    }
    constexpr std::size_t local_pair_count = product_pairs.size();
    weighted_pair_density.resize(system.atoms.size() * local_pair_count);
    for (std::size_t atom = 0; atom < system.atoms.size(); ++atom) {
        const std::size_t first = system.first_orbital[atom];
        const std::size_t orbital_count = static_cast<std::size_t>(
            system.parameters[atom]->basis_orbitals);
        const std::size_t pair_count =
            orbital_count * (orbital_count + 1) / 2;
        double* const packed = weighted_pair_density.data()
            + atom * local_pair_count;
        for (std::size_t pair = 0; pair < pair_count; ++pair) {
            const auto [high, low] = product_pairs[pair];
            const std::size_t row = first + static_cast<std::size_t>(high);
            const std::size_t column = first + static_cast<std::size_t>(low);
            const double multiplicity = high == low ? 1.0 : 2.0;
            packed[pair] = multiplicity * density[row * n + column];
        }
    }
    for (const auto& block : system.integral_blocks) {
        const std::size_t fi = system.first_orbital[block.atom_i];
        const std::size_t fj = system.first_orbital[block.atom_j];
        const std::size_t ni = block.orbital_count_i * (block.orbital_count_i + 1) / 2;
        const std::size_t nj = block.orbital_count_j * (block.orbital_count_j + 1) / 2;
        const double* const density_i = weighted_pair_density.data()
            + block.atom_i * local_pair_count;
        const double* const density_j = weighted_pair_density.data()
            + block.atom_j * local_pair_count;
        std::array<double, local_pair_count> row_sums{};
        std::array<double, local_pair_count> column_sums{};
        for (std::size_t row = 0; row < ni; ++row) {
            const double* const integrals = block.values.data() + row * nj;
            for (std::size_t column = 0; column < nj; ++column) {
                row_sums[row] += density_j[column] * integrals[column];
                column_sums[column] += density_i[row] * integrals[column];
            }
        }
        for (std::size_t row = 0; row < ni; ++row) {
            const auto [ml, nl] = product_pairs[row];
            const std::size_t mu = fi + static_cast<std::size_t>(ml);
            const std::size_t nu = fi + static_cast<std::size_t>(nl);
            fock[mu * n + nu] += row_sums[row];
            if (mu != nu) fock[nu * n + mu] += row_sums[row];
        }
        for (std::size_t column = 0; column < nj; ++column) {
            const auto [ll, sl] = product_pairs[column];
            const std::size_t lam = fj + static_cast<std::size_t>(ll);
            const std::size_t sig = fj + static_cast<std::size_t>(sl);
            fock[lam * n + sig] += column_sums[column];
            if (lam != sig) {
                fock[sig * n + lam] += column_sums[column];
            }
        }
        for (std::size_t ml = 0; ml < block.orbital_count_i; ++ml) {
            const std::size_t mu = fi + ml;
            for (std::size_t nl = 0; nl < block.orbital_count_j; ++nl) {
                const std::size_t nu = fj + nl;
                double value = 0.0;
                for (std::size_t ll = 0; ll < block.orbital_count_i; ++ll) {
                    const std::size_t lam = fi + ll;
                    const std::size_t row = orbital_pair_indices[ml][ll];
                    for (std::size_t sl = 0; sl < block.orbital_count_j; ++sl) {
                        const std::size_t sig = fj + sl;
                        const std::size_t column =
                            orbital_pair_indices[nl][sl];
                        value -= .5 * density[lam * n + sig]
                            * block.values[row * nj + column];
                    }
                }
                fock[mu * n + nu] += value;
                fock[nu * n + mu] += value;
            }
        }
    }
}

struct Charges {
    std::vector<double> populations, mulliken, bond, cm2, response;
};

void charges(const ElectronicSystem& system, std::span<const double> density,
             Charges& result) {
    AMSOLCPP_PERF_SCOPE(PerformanceStage::Cm2);
    AMSOLCPP_PERF_INCREMENT(cm2_evaluations);
    const std::size_t atoms = system.atoms.size(), n = system.orbital_count;
    result.populations.assign(atoms, 0.0);
    result.bond.assign(atoms * atoms, 0.0);
    for (std::size_t ai = 0; ai < atoms; ++ai) {
        const auto first_i = system.first_orbital[ai]; const int count_i = system.parameters[ai]->basis_orbitals;
        for (int u = 0; u < count_i; ++u) { const auto ao = first_i + static_cast<std::size_t>(u); result.populations[ai] += density[ao * n + ao]; }
        for (std::size_t aj = 0; aj < ai; ++aj) {
            const auto first_j = system.first_orbital[aj]; const int count_j = system.parameters[aj]->basis_orbitals; double value = 0.0;
            for (int u = 0; u < count_i; ++u) for (int v = 0; v < count_j; ++v) { const double p = density[(first_i + static_cast<std::size_t>(u)) * n + first_j + static_cast<std::size_t>(v)]; value += p * p; }
            result.bond[ai * atoms + aj] = value; result.bond[aj * atoms + ai] = value;
            result.bond[ai * atoms + ai] += value; result.bond[aj * atoms + aj] += value;
        }
    }
    result.mulliken.resize(atoms); result.cm2.resize(atoms);
    for (std::size_t i = 0; i < atoms; ++i) result.mulliken[i] = result.cm2[i] = static_cast<double>(system.parameters[i]->core_charge) - result.populations[i];
    for (std::size_t i = 0; i < atoms; ++i) for (std::size_t j = 0; j < i; ++j) {
        const int zi = system.atoms[i].atomic_number, zj = system.atoms[j].atomic_number;
        const double b = result.bond[i * atoms + j];
        const double delta = b * (cm2_parameter(zi, zj, false) + cm2_parameter(zi, zj, true) * b);
        result.cm2[i] += delta; result.cm2[j] -= delta;
    }
    result.response.resize(atoms * atoms);
    for (std::size_t i = 0; i < atoms; ++i) {
        for (std::size_t j = 0; j < atoms; ++j) {
            const int zi = system.atoms[i].atomic_number;
            const int zj = system.atoms[j].atomic_number;
            result.response[i * atoms + j] =
                cm2_parameter(zi, zj, false)
                + 2.0 * cm2_parameter(zi, zj, true)
                    * result.bond[i * atoms + j];
        }
    }
}

void add_polar_fock(const ElectronicSystem& system,
                    std::span<const double> density,
                    const Charges& charge,
                    std::span<const double> fgb,
                    std::vector<double>& fock,
                    std::vector<double>& vk) {
    AMSOLCPP_PERF_SCOPE(PerformanceStage::PolarFock);
    AMSOLCPP_PERF_INCREMENT(polar_fock_constructions);
    const std::size_t na = system.atoms.size(), n = system.orbital_count;
    vk.assign(na, 0.0);
    for (std::size_t i = 0; i < na; ++i) for (std::size_t j = 0; j < na; ++j) vk[i] += fgb[i * na + j] * charge.cm2[j];
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t ai = system.ao_to_atom[i];
        for (std::size_t j = 0; j <= i; ++j) {
            const std::size_t aj = system.ao_to_atom[j]; double value = i == j ? vk[ai] : 0.0;
            // The legacy delta terms below can contribute only at k == j and
            // k == i.  Evaluate those sites in the original order instead of
            // scanning every AO and performing unused CM2 table lookups.
            const double dji = charge.response[aj * na + ai];
            value -= .5 * (vk[aj] - vk[ai]) * dji * density[j * n + i];

            const double dij = charge.response[ai * na + aj];
            value -= .5 * (vk[ai] - vk[aj]) * dij * density[i * n + j];
            fock[i * n + j] += value;
            if (i != j) {
                fock[j * n + i] += value;
            }
        }
    }
}

void density_from_vectors(std::span<const double> vectors,
                          const std::size_t n,
                          const std::size_t occupied,
                          std::vector<double>& density) {
    AMSOLCPP_PERF_SCOPE(PerformanceStage::DensityConstruction);
    AMSOLCPP_PERF_INCREMENT(density_rebuilds);
    density.resize(n * n);
    for (std::size_t i = 0; i < n; ++i) {
        double* const density_row = density.data() + i * n;
        std::fill_n(density_row, i + 1, 0.0);
        for (std::size_t q = 0; q < occupied; ++q) {
            const double* const orbital = vectors.data() + q * n;
            const double doubled_coefficient = 2.0 * orbital[i];
            for (std::size_t j = 0; j <= i; ++j) {
                density_row[j] += doubled_coefficient * orbital[j];
            }
        }
        for (std::size_t j = 0; j < i; ++j) {
            density[j * n + i] = density_row[j];
        }
    }
}

double solution_energy(std::span<const double> density, std::span<const double> hcore,
                       std::span<const double> gas, std::span<const double> q,
                       std::span<const double> fgb) {
    AMSOLCPP_PERF_SCOPE(PerformanceStage::ElectronicEnergy);
    AMSOLCPP_PERF_INCREMENT(energy_evaluations);
    double value = 0.0; for (std::size_t i = 0; i < density.size(); ++i) value += .5 * density[i] * (hcore[i] + gas[i]);
    const std::size_t na = q.size(); for (std::size_t i = 0; i < na; ++i) for (std::size_t j = 0; j < na; ++j) value -= .5 * q[i] * fgb[i * na + j] * q[j]; return value;
}

double commutator(std::span<const double> fock, std::span<const double> density, const std::size_t n, std::vector<double>* output = nullptr) {
    AMSOLCPP_PERF_SCOPE(PerformanceStage::ConvergenceResidual);
    AMSOLCPP_PERF_INCREMENT(commutator_evaluations);
    if (output) {
        output->resize(n * n);
        for (std::size_t i = 0; i < n; ++i) {
            (*output)[i * n + i] = 0.0;
        }
    }
    double largest = 0.0;
    // F and P are symmetric, so FP-PF is skew-symmetric.  AMSOL 7.1
    // forms FP once and derives the packed commutator from its transpose;
    // compute one triangle here and mirror it instead of evaluating both.
#ifdef _OPENMP
// Row i has i*n contraction work, so cyclic static assignment balances the
// triangle without changing any individual k-summation order.
#pragma omp parallel for schedule(static, 1) reduction(max : largest) \
    if(n >= openmp_commutator_min_dimension)
#endif
    for (std::size_t i = 1; i < n; ++i) {
        for (std::size_t j = 0; j < i; ++j) {
            double value = 0.0;
            for (std::size_t k = 0; k < n; ++k) {
                value += fock[i * n + k] * density[j * n + k]
                    - density[i * n + k] * fock[j * n + k];
            }
            largest = std::max(largest, std::abs(value));
            if (output) {
                (*output)[i * n + j] = value;
                (*output)[j * n + i] = -value;
            }
        }
    }
    return largest;
}

void pseudo_diagonalize(std::span<const double> fock, std::vector<double>& vectors,
                        std::span<const double> eigenvalues,
                        const std::size_t n, const std::size_t occupied,
                        std::vector<double>& fock_times_virtual,
                        std::vector<double>& fmo) {
    AMSOLCPP_PERF_SCOPE(PerformanceStage::PseudoDiagonalization);
    AMSOLCPP_PERF_INCREMENT(pseudo_diagonalizations);
    const std::size_t virtual_count = n - occupied;
    fock_times_virtual.resize(virtual_count * n);
    fmo.resize(virtual_count * occupied);
    double largest = -std::numeric_limits<double>::infinity();

    // Match the factored AMSOL 7.1 DIAG path: form F*C_virtual once,
    // then contract C_occupied^T with that workspace.  The previous C++
    // expression expanded the same product independently for every
    // occupied/virtual pair, increasing the work from O(n^3) to O(n^4).
    for (std::size_t virtual_index = 0; virtual_index < virtual_count;
         ++virtual_index) {
        const std::size_t v = occupied + virtual_index;
        double* const transformed =
            fock_times_virtual.data() + virtual_index * n;
        for (std::size_t i = 0; i < n; ++i) {
            double value = 0.0;
            for (std::size_t j = 0; j < n; ++j) {
                value += fock[i * n + j] * vectors[v * n + j];
            }
            transformed[i] = value;
        }
        for (std::size_t o = 0; o < occupied; ++o) {
            double value = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                value += vectors[o * n + i] * transformed[i];
            }
            fmo[virtual_index * occupied + o] = value;
            largest = std::max(largest, value);
        }
    }
    // AMSOL 7.1 calls its local ISMAX here.  Despite the inherited IDAMAX
    // comment, new/ismax.f selects the largest signed FMO entry and DIAG only
    // then takes its absolute value.  Selecting the largest magnitude changes
    // which occupied/virtual rotations are retained and can select a different
    // solvent SCF root on difficult molecules.
    const double tiny = .04 * std::abs(largest);
    for(std::size_t v=occupied;v<n;++v)for(std::size_t o=0;o<occupied;++o){const double coupling=fmo[(v-occupied)*occupied+o];if(std::abs(coupling)<tiny)continue;const double difference=eigenvalues[v]-eigenvalues[o];const double fraction=std::min(1.0,.5*(1.0+difference/std::sqrt(4*coupling*coupling+difference*difference)));const double alpha=std::sqrt(fraction),beta=std::copysign(std::sqrt(1-fraction),coupling);for(std::size_t row=0;row<n;++row){const double ov=vectors[o*n+row],vv=vectors[v*n+row];vectors[o*n+row]=alpha*ov-beta*vv;vectors[v*n+row]=alpha*vv+beta*ov;}}
}

bool solve_linear(std::vector<double>& a, std::vector<double>& b, const std::size_t n) {
    for(std::size_t k=0;k<n;++k){std::size_t pivot=k;for(std::size_t i=k+1;i<n;++i)if(std::abs(a[i*n+k])>std::abs(a[pivot*n+k]))pivot=i;if(std::abs(a[pivot*n+k])<1e-14)return false;if(pivot!=k){for(std::size_t j=k;j<n;++j)std::swap(a[k*n+j],a[pivot*n+j]);std::swap(b[k],b[pivot]);}for(std::size_t i=k+1;i<n;++i){const double f=a[i*n+k]/a[k*n+k];for(std::size_t j=k;j<n;++j)a[i*n+j]-=f*a[k*n+j];b[i]-=f*b[k];}}
    for(std::size_t ii=n;ii-->0;){double x=b[ii];for(std::size_t j=ii+1;j<n;++j)x-=a[ii*n+j]*b[j];b[ii]=x/a[ii*n+ii];}return true;
}

std::vector<double> diis_fock(
    const std::vector<std::vector<double>>& focks,
    const std::vector<std::vector<double>>& errors
) {
    AMSOLCPP_PERF_SCOPE(PerformanceStage::Diis);
    AMSOLCPP_PERF_INCREMENT(diis_extrapolations);
    const std::size_t count = focks.size();
    const std::size_t m = count + 1;
    std::vector<double> a(m * m, 0.0);
    std::vector<double> b(m, 0.0);
    b[count] = -1.0;
    for (std::size_t i = 0; i < count; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            double dot = 0.0;
            for (std::size_t k = 0; k < errors[i].size(); ++k) {
                dot += errors[i][k] * errors[j][k];
            }
            a[i * m + j] = dot;
            a[j * m + i] = dot;
        }
        a[i * m + count] = -1.0;
        a[count * m + i] = -1.0;
    }
    if (!solve_linear(a, b, m)) {
        return focks.back();
    }
    std::vector<double> result(focks.back().size(), 0.0);
    for (std::size_t i = 0; i < count; ++i) {
        for (std::size_t k = 0; k < result.size(); ++k) {
            result[k] += b[i] * focks[i][k];
        }
    }
    return result;
}

ScfState make_scf_state(
    std::vector<double> density,
    Charges charge,
    SymmetricEigensystem eigensystem,
    const std::size_t iterations,
    const bool rescue,
    const bool damping,
    const double electronic_energy_ev,
    const double energy_residual,
    const double density_residual,
    const double commutator_residual,
    std::vector<ScfIterationRecord> records
) {
    ScfState state;
    state.converged = true;
    state.rescue_used = rescue;
    state.damping_used = damping;
    state.iterations = iterations;
    state.electronic_energy_ev = electronic_energy_ev;
    state.energy_residual_kcal = energy_residual;
    state.density_residual = density_residual;
    state.commutator_residual = commutator_residual;
    state.density = std::move(density);
    state.orbital_energies = std::move(eigensystem.eigenvalues);
    state.mulliken_charges = std::move(charge.mulliken);
    state.cm2_charges = std::move(charge.cm2);
    state.bond_orders = std::move(charge.bond);
    state.records = std::move(records);
    return state;
}

ScfState finalize_precomputed(
    const ElectronicSystem& system,
    std::vector<double> density,
    Charges charge,
    std::span<const double> final_fock,
    std::vector<double>& candidate_workspace,
    const std::size_t iterations,
    const bool rescue,
    const bool damping,
    const double electronic_energy_ev,
    const double energy_residual,
    std::vector<ScfIterationRecord> records
) {
    AMSOLCPP_PERF_SCOPE(PerformanceStage::FinalConsistency);
    AMSOLCPP_PERF_INCREMENT(final_consistency_attempts);
    auto eigensystem = diagonalize_amsol_symmetric_scf(
        final_fock, system.orbital_count);
    density_from_vectors(
        eigensystem.eigenvectors, system.orbital_count,
        system.occupied_orbitals, candidate_workspace);
    double density_residual = 0.0;
    for (std::size_t i = 0; i < density.size(); ++i) {
        density_residual = std::max(
            density_residual,
            std::abs(candidate_workspace[i] - density[i])
        );
    }

    const double commutator_residual = commutator(
        final_fock, density, system.orbital_count);
    return make_scf_state(
        std::move(density), std::move(charge), std::move(eigensystem),
        iterations, rescue, damping, electronic_energy_ev, energy_residual,
        density_residual, commutator_residual, std::move(records));
}

ScfState finalize_reused(
    std::vector<double> density,
    Charges charge,
    SymmetricEigensystem eigensystem,
    const std::size_t iterations,
    const double electronic_energy_ev,
    const double energy_residual,
    const double density_residual,
    const double commutator_residual,
    std::vector<ScfIterationRecord> records
) {
    AMSOLCPP_PERF_SCOPE(PerformanceStage::FinalConsistency);
    AMSOLCPP_PERF_INCREMENT(final_consistency_attempts);
    return make_scf_state(
        std::move(density), std::move(charge), std::move(eigensystem),
        iterations, true, true, electronic_energy_ev, energy_residual,
        density_residual, commutator_residual, std::move(records));
}

}  // namespace

ScfState run_scf(const ElectronicSystem& system, std::span<const double> fgb, const CalculationOptions& options) {
    AMSOLCPP_PERF_SCOPE(PerformanceStage::ScfOverall);
    const std::size_t na = system.atoms.size();
    const std::size_t n = system.orbital_count;
    if (fgb.size() != na * na) {
        throw InputError("FGB matrix dimensions do not match the molecule");
    }

    auto density = initial_density(system, options.molecular_charge);
    auto initial_reference = density;
    for (double& value : initial_reference) {
        value *= 0.5;
    }
    std::vector<double> old_diag(n, 0.0);
    std::vector<double> old_current(n);
    std::vector<double> diff(n);
    std::vector<double> candidate;
    std::vector<double> accepted;
    std::vector<double> fixed_density;
    std::vector<double> vectors;
    std::vector<double> eigenvalues;
    std::vector<double> pseudo_fock_times_virtual;
    std::vector<double> pseudo_fmo;
    std::vector<double> primary_fock;
    std::vector<double> weighted_pair_density;
    std::vector<double> polar_potential;
    Charges charge_workspace;
    std::vector<ScfIterationRecord> records;
    bool pseudo = false;
    bool damping = false;
    bool ready = false;
    std::size_t ready_count = 0;
    double pl = 1.0;
    double previous_energy = 100.0;
    double energy_delta = std::numeric_limits<double>::infinity();
    double last_comm = 0.0;
    double last_change = 0.0;
    const double selcon = 1.0e-6 * ev_to_kcal;
    const double pltest = std::pow(
        10.0, 0.605 * std::log10(selcon) - 1.6129);
    const std::size_t primary_limit = std::min(
        options.legacy_scf_iterations, options.max_scf_iterations);
    std::size_t primary_iterations = 0;
    bool legacy_root_available = false;

    for (std::size_t iteration = 1; iteration <= primary_limit; ++iteration) {
        AMSOLCPP_PERF_INCREMENT(primary_scf_iterations);
        gas_fock(system, density, primary_fock, weighted_pair_density);
        charges(system, density, charge_workspace);
        const double electronic_energy = solution_energy(
            density, system.hcore, primary_fock,
            charge_workspace.cm2, fgb);
        add_polar_fock(
            system, density, charge_workspace, fgb, primary_fock,
            polar_potential);
        if (options.collect_iteration_diagnostics) {
            last_comm = commutator(primary_fock, density, n);
        }
        const double energy = electronic_energy * ev_to_kcal;
        energy_delta = energy - previous_energy;

        if (pl < pltest && std::abs(energy_delta) < selcon && ready) {
            // The exact Fock, charges, and energy for this density were just
            // computed.  Reuse them for the strict final consistency check.
            auto legacy_state = finalize_precomputed(
                system, std::move(density), std::move(charge_workspace),
                primary_fock, candidate, iteration, false, damping,
                electronic_energy, std::abs(energy_delta),
                std::move(records));
            if (legacy_state.density_residual <= options.density_tolerance
                && legacy_state.commutator_residual
                    <= options.commutator_tolerance
                && legacy_state.energy_residual_kcal
                    <= options.energy_tolerance * ev_to_kcal) {
                return legacy_state;
            }
            // ITER/CNVG selected a source-valid root, but its final density is
            // not yet tight enough for the stronger C++ consistency contract.
            // Continue from that accepted root.  Restarting from PDIAG here
            // discards AMSOL's root selection and lets DIIS converge tightly
            // to a different stationary determinant on difficult polar cases.
            density = std::move(legacy_state.density);
            legacy_root_available = true;
            records = std::move(legacy_state.records);
            primary_iterations = iteration;
            break;
        }

        ready = ready_count > 0
            && std::abs(energy_delta) < selcon * 10.0;
        ++ready_count;
        previous_energy = energy;
        if (pseudo) {
            pseudo_diagonalize(
                primary_fock, vectors, eigenvalues, n,
                system.occupied_orbitals, pseudo_fock_times_virtual,
                pseudo_fmo);
        } else {
            auto eig = diagonalize_amsol_symmetric_scf(primary_fock, n);
            eigenvalues = std::move(eig.eigenvalues);
            vectors = std::move(eig.eigenvectors);
        }
        density_from_vectors(
            vectors, n, system.occupied_orbitals, candidate);
        {
            AMSOLCPP_PERF_SCOPE(PerformanceStage::ConvergenceEvaluation);
            const auto& reference = iteration == 1
                ? initial_reference : density;
            pl = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                old_current[i] = reference[i * n + i];
                diff[i] = std::abs(
                    candidate[i * n + i] - old_current[i]);
                pl = std::max(pl, diff[i]);
            }
            double fac = 0.0;
            if (iteration % 3 == 0) {
                double a = 0.0;
                double b = 0.0;
                for (std::size_t i = 0; i < n; ++i) {
                    a += diff[i] * diff[i];
                    const double q = candidate[i * n + i]
                        - 2.0 * old_current[i] + old_diag[i];
                    b += q * q;
                }
                if (b > 0.0 && a < 100.0 * b) {
                    fac = std::sqrt(a / b);
                }
            }
            for (std::size_t i = 0; i < candidate.size(); ++i) {
                const double value = candidate[i];
                candidate[i] += fac * (value - reference[i]);
            }
            const double damp = iteration <= 3 ? 1.0e10 : 0.05;
            if (iteration > 3) {
                damping = true;
            }
            for (std::size_t i = 0; i < n; ++i) {
                double value = candidate[i * n + i];
                const double step = value - old_current[i];
                if (std::abs(step) > damp) {
                    value = old_current[i] + std::copysign(damp, step);
                }
                candidate[i * n + i] = std::clamp(value, 0.0, 2.0);
            }
            last_change = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t j = 0; j <= i; ++j) {
                    const std::size_t index = i * n + j;
                    last_change = std::max(
                        last_change,
                        std::abs(candidate[index] - density[index]));
                }
            }
            if (options.collect_iteration_diagnostics) {
                records.push_back({iteration, electronic_energy, energy_delta,
                                   last_change, last_comm, fac});
            }
            old_diag.swap(old_current);
            density.swap(candidate);
            pseudo = pl < 0.1 || pseudo;
        }
    }

    if (!options.enable_scf_rescue) {
        throw ScfConvergenceError("legacy AM1 SCF did not converge");
    }
    if (primary_iterations == 0) {
        primary_iterations = primary_limit;
    }
    // Polish a source-converged root in place.  Only a genuinely exhausted
    // legacy trajectory starts the independent bounded rescue from PDIAG.
    if (!legacy_root_available) {
        density = initial_density(system, options.molecular_charge);
    }
    const std::size_t limit = options.max_scf_iterations - primary_iterations;
    if (limit == 0) {
        throw ScfConvergenceError(
            "legacy AM1 SCF reached the configured iteration limit");
    }

    std::vector<std::vector<double>> focks;
    std::vector<std::vector<double>> errors;
    focks.reserve(6);
    errors.reserve(6);
    std::vector<double> carried_fock;
    std::vector<double> carried_error;
    std::vector<double> reusable_candidate;
    bool have_carried_state = false;
    bool have_reusable_candidate = false;
    double previous_ev = std::numeric_limits<double>::quiet_NaN();
    for (std::size_t iteration = 1; iteration <= limit; ++iteration) {
        AMSOLCPP_PERF_INCREMENT(rescue_scf_iterations);
        std::vector<double> error;
        std::vector<double> fock;
        if (have_carried_state) {
            fock = std::move(carried_fock);
            error = std::move(carried_error);
            have_carried_state = false;
        } else {
            gas_fock(system, density, fock, weighted_pair_density);
            charges(system, density, charge_workspace);
            add_polar_fock(
                system, density, charge_workspace, fgb, fock,
                polar_potential);
            commutator(fock, density, n, &error);
        }
        if (focks.size() == 6) {
            focks.erase(focks.begin());
            errors.erase(errors.begin());
        }
        focks.push_back(std::move(fock));
        errors.push_back(std::move(error));

        if (iteration == 2 && have_reusable_candidate) {
            // Before DIIS begins, this is exactly the fixed density already
            // constructed from iteration one's accepted Fock.
            candidate.swap(reusable_candidate);
            have_reusable_candidate = false;
        } else if (iteration >= 3) {
            auto diagonal = diis_fock(focks, errors);
            auto eig = diagonalize_amsol_symmetric_scf(diagonal, n);
            density_from_vectors(
                eig.eigenvectors, n, system.occupied_orbitals, candidate);
        } else {
            auto eig = diagonalize_amsol_symmetric_scf(focks.back(), n);
            density_from_vectors(
                eig.eigenvectors, n, system.occupied_orbitals, candidate);
        }
        double mixing = 0.0;
        {
            AMSOLCPP_PERF_SCOPE(PerformanceStage::ConvergenceEvaluation);
            double raw = 0.0;
            for (std::size_t i = 0; i < density.size(); ++i) {
                raw = std::max(raw, std::abs(candidate[i] - density[i]));
            }
            mixing = raw > 0.25 ? 0.2 : 1.0;
            accepted.assign(density.begin(), density.end());
            for (std::size_t i = 0; i < density.size(); ++i) {
                accepted[i] += mixing * (candidate[i] - density[i]);
            }
        }

        std::vector<double> accepted_fock;
        gas_fock(system, accepted, accepted_fock, weighted_pair_density);
        charges(system, accepted, charge_workspace);
        const double ev = solution_energy(
            accepted, system.hcore, accepted_fock,
            charge_workspace.cm2, fgb);
        add_polar_fock(
            system, accepted, charge_workspace, fgb, accepted_fock,
            polar_potential);
        std::vector<double> accepted_error;
        last_comm = commutator(
            accepted_fock, accepted, n, &accepted_error);
        auto accepted_eig = diagonalize_amsol_symmetric_scf(accepted_fock, n);
        density_from_vectors(
            accepted_eig.eigenvectors, n, system.occupied_orbitals,
            fixed_density);
        {
            AMSOLCPP_PERF_SCOPE(PerformanceStage::ConvergenceEvaluation);
            last_change = 0.0;
            for (std::size_t i = 0; i < accepted.size(); ++i) {
                last_change = std::max(
                    last_change, std::abs(fixed_density[i] - accepted[i]));
            }
            energy_delta = std::isfinite(previous_ev)
                ? (ev - previous_ev) * ev_to_kcal
                : std::numeric_limits<double>::infinity();
            previous_ev = ev;
            const std::size_t total_iteration = primary_iterations + iteration;
            if (options.collect_iteration_diagnostics) {
                records.push_back({total_iteration, ev, energy_delta,
                                   last_change, last_comm, mixing});
            }
            if (last_change <= options.density_tolerance
                && last_comm <= options.commutator_tolerance
                && std::abs(energy_delta)
                    <= options.energy_tolerance * ev_to_kcal) {
                // The accepted state has already undergone every strict final
                // consistency calculation.  Return those exact artifacts.
                return finalize_reused(
                    std::move(accepted), std::move(charge_workspace),
                    std::move(accepted_eig), total_iteration, ev,
                    std::abs(energy_delta), last_change, last_comm,
                    std::move(records));
            }
            if (iteration == 1) {
                reusable_candidate.swap(fixed_density);
                have_reusable_candidate = true;
            }
            density.swap(accepted);
            carried_fock = std::move(accepted_fock);
            carried_error = std::move(accepted_error);
            have_carried_state = true;
        }
    }
    throw ScfConvergenceError("AM1 SCF failed both legacy CNVG and bounded DIIS convergence paths");
}

}  // namespace amsolcpp::detail
