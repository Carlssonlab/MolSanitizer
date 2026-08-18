#pragma once

#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>

namespace amsolcpp::detail {

enum class PerformanceStage : std::size_t {
    Parsing,
    GeometryPreprocessing,
    BasisAndParameters,
    OneElectronTerms,
    TwoElectronIntegrals,
    IntegralPreprocessing,
    InitialDensity,
    ScfOverall,
    GasFock,
    Cm2,
    PolarFock,
    ElectronicEnergy,
    Diagonalization,
    EigensolverReduction,
    EigensolverQl,
    PseudoDiagonalization,
    DensityConstruction,
    Diis,
    ConvergenceResidual,
    ConvergenceEvaluation,
    FinalConsistency,
    SolvationOverall,
    SurfaceGeometry,
    Cds,
    ResultAssembly,
    OutputFormatting,
    Unattributed,
    Count,
};

constexpr std::size_t performance_stage_count =
    static_cast<std::size_t>(PerformanceStage::Count);

struct PerformanceDiagnostics {
    std::array<std::uint64_t, performance_stage_count> nanoseconds{};
    std::array<std::uint64_t, performance_stage_count> calls{};
    std::array<std::uint64_t, performance_stage_count> allocations{};
    std::array<std::uint64_t, performance_stage_count> allocated_bytes{};

    std::uint64_t integral_pair_evaluations{};
    std::uint64_t gas_fock_constructions{};
    std::uint64_t polar_fock_constructions{};
    std::uint64_t full_diagonalizations{};
    std::uint64_t pseudo_diagonalizations{};
    std::uint64_t jacobi_sweeps{};
    std::uint64_t jacobi_rotations{};
    std::uint64_t density_rebuilds{};
    std::uint64_t cm2_evaluations{};
    std::uint64_t commutator_evaluations{};
    std::uint64_t energy_evaluations{};
    std::uint64_t diis_extrapolations{};
    std::uint64_t final_consistency_attempts{};
    std::uint64_t primary_scf_iterations{};
    std::uint64_t rescue_scf_iterations{};
};

#ifdef AMSOLCPP_PERFORMANCE_DIAGNOSTICS

PerformanceDiagnostics& active_performance_diagnostics() noexcept;
void reset_performance_diagnostics() noexcept;
PerformanceDiagnostics performance_diagnostics_snapshot() noexcept;
const char* performance_stage_name(PerformanceStage stage) noexcept;
void set_performance_allocation_tracking(bool enabled) noexcept;
void record_performance_allocation(std::size_t bytes) noexcept;

class PerformanceScope {
public:
    explicit PerformanceScope(PerformanceStage stage) noexcept;
    ~PerformanceScope() noexcept;

    PerformanceScope(const PerformanceScope&) = delete;
    PerformanceScope& operator=(const PerformanceScope&) = delete;

private:
    PerformanceStage stage_;
    PerformanceStage previous_stage_;
    std::chrono::steady_clock::time_point start_;
};

#define AMSOLCPP_PERF_JOIN_INNER(left, right) left##right
#define AMSOLCPP_PERF_JOIN(left, right) AMSOLCPP_PERF_JOIN_INNER(left, right)
#define AMSOLCPP_PERF_SCOPE(stage) \
    ::amsolcpp::detail::PerformanceScope \
        AMSOLCPP_PERF_JOIN(amsolcpp_performance_scope_, __LINE__){stage}
#define AMSOLCPP_PERF_INCREMENT(field) \
    ++::amsolcpp::detail::active_performance_diagnostics().field

#else

#define AMSOLCPP_PERF_SCOPE(stage) static_cast<void>(stage)
#define AMSOLCPP_PERF_INCREMENT(field) static_cast<void>(0)

#endif

}  // namespace amsolcpp::detail
