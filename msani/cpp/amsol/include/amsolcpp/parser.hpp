#pragma once

#include <string>
#include <string_view>
#include <vector>

#include <amsolcpp/atom.hpp>
#include <amsolcpp/options.hpp>

namespace amsolcpp {

struct ParsedInput {
    std::string molecule_name;
    std::vector<Atom> atoms;
    CalculationOptions options;
};

[[nodiscard]] ParsedInput parse_legacy_input(std::string_view input);

}  // namespace amsolcpp
