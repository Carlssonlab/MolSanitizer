#pragma once

#include <span>
#include <string>
#include <string_view>
#include <vector>

#include <amsolcpp/atom.hpp>
#include <amsolcpp/exceptions.hpp>
#include <amsolcpp/options.hpp>
#include <amsolcpp/result.hpp>

namespace amsolcpp {

[[nodiscard]] CalculationResult calculate(
    std::span<const Atom> atoms,
    const CalculationOptions& options = {}
);

[[nodiscard]] DualSolventResult calculate_water_and_hexadecane(
    std::span<const Atom> atoms,
    CalculationOptions options = {}
);

[[nodiscard]] CalculationResult calculate_from_input(std::string_view input);
[[nodiscard]] std::string run_legacy_input(std::string_view input);

[[nodiscard]] std::vector<CalculationResult> calculate_batch(
    const std::vector<std::vector<Atom>>& molecules,
    const CalculationOptions& options = {},
    bool continue_on_error = false
);

}  // namespace amsolcpp
