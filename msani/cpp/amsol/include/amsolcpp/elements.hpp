#pragma once

#include <array>
#include <cstddef>
#include <span>
#include <string_view>

namespace amsolcpp {

struct GaussianCoreTerm {
    double amplitude{};
    double exponent{};
    double displacement{};
};

struct ElementParameters {
    int atomic_number{};
    std::string_view symbol;
    int core_charge{};
    int valence_electrons{};
    int basis_orbitals{};
    double uss{};
    double upp{};
    double betas{};
    double betap{};
    double zs{};
    double zp{};
    double alpha{};
    double isolated_atom_energy{};
    double gss{};
    double gsp{};
    double gpp{};
    double gp2{};
    double hsp{};
    double dd{};
    double qq{};
    double am{};
    double ad{};
    double aq{};
    std::array<GaussianCoreTerm, 4> gaussian_terms{};
    std::size_t gaussian_term_count{};
};

[[nodiscard]] const ElementParameters& element_parameters(int atomic_number);
[[nodiscard]] const ElementParameters& element_parameters(std::string_view symbol);
[[nodiscard]] bool is_supported_element(int atomic_number) noexcept;
[[nodiscard]] int atomic_number(std::string_view symbol);
[[nodiscard]] std::string_view element_symbol(int atomic_number);
[[nodiscard]] std::span<const int> supported_atomic_numbers() noexcept;

}  // namespace amsolcpp
