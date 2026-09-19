#pragma once

#include <span>
#include <string>
#include <string_view>

#include <amsolcpp/atom.hpp>
#include <amsolcpp/options.hpp>
#include <amsolcpp/result.hpp>

namespace amsolcpp {

[[nodiscard]] std::string result_to_json(
    std::span<const Atom> atoms,
    const CalculationOptions& options,
    const CalculationResult& result
);

[[nodiscard]] std::string format_legacy_output(
    std::string_view molecule_name,
    std::span<const Atom> atoms,
    const CalculationResult& result
);

}  // namespace amsolcpp
