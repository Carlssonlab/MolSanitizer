#include <amsolcpp/parser.hpp>

#include "performance_diagnostics.hpp"

#include <amsolcpp/coordinates.hpp>
#include <amsolcpp/elements.hpp>
#include <amsolcpp/exceptions.hpp>

#include <algorithm>
#include <array>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace amsolcpp {
namespace {

std::string trim(const std::string_view value) {
    std::size_t first = 0;
    while (first < value.size() && std::isspace(static_cast<unsigned char>(value[first]))) {
        ++first;
    }
    std::size_t last = value.size();
    while (last > first && std::isspace(static_cast<unsigned char>(value[last - 1]))) {
        --last;
    }
    return std::string(value.substr(first, last - first));
}

std::string upper(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](const char c) {
        return static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    });
    return value;
}

std::vector<std::string> split_words(const std::string_view value) {
    std::istringstream stream{std::string(value)};
    std::vector<std::string> result;
    for (std::string word; stream >> word;) {
        result.push_back(std::move(word));
    }
    return result;
}

int parse_integer(const std::string_view value, const std::string_view label) {
    int result = 0;
    const auto conversion = std::from_chars(value.data(), value.data() + value.size(), result);
    if (conversion.ec != std::errc{} || conversion.ptr != value.data() + value.size()) {
        throw InputError(std::string(label) + " must be an integer");
    }
    return result;
}

double parse_double(const std::string& value, const std::string_view label) {
    char* end = nullptr;
    const double result = std::strtod(value.c_str(), &end);
    if (end != value.c_str() + value.size() || !std::isfinite(result)) {
        throw InputError(std::string(label) + " must be a finite number");
    }
    return result;
}

bool molecule_header(const std::vector<std::string>& words) {
    if (words.size() != 2) {
        return false;
    }
    try {
        static_cast<void>(parse_integer(words[1], "atom count"));
        return true;
    } catch (const InputError&) {
        return false;
    }
}

void require_close(
    const std::unordered_map<std::string, std::string>& values,
    const std::string& key,
    const double expected
) {
    const auto found = values.find(key);
    if (found == values.end()) {
        throw UnsupportedCalculationError("GENORG hexadecane input requires " + key);
    }
    const double actual = parse_double(found->second, key);
    if (std::abs(actual - expected) > 1.0e-12) {
        throw UnsupportedCalculationError("GENORG " + key + " differs from the audited hexadecane value");
    }
}

}  // namespace

ParsedInput parse_legacy_input(const std::string_view input) {
    AMSOLCPP_PERF_SCOPE(detail::PerformanceStage::Parsing);
    constexpr std::size_t maximum_input_bytes = 16U * 1024U * 1024U;
    if (input.size() > maximum_input_bytes) {
        throw InputError("AMSOL input exceeds the 16 MiB parser limit");
    }
    std::vector<std::string> lines;
    std::istringstream input_stream{std::string(input)};
    for (std::string line; std::getline(input_stream, line);) {
        if (line.size() > 1024U * 1024U) {
            throw InputError("AMSOL input line exceeds the 1 MiB parser limit");
        }
        lines.push_back(std::move(line));
    }
    if (lines.empty()) {
        throw InputError("empty AMSOL input");
    }

    std::unordered_set<std::string> bare;
    std::unordered_map<std::string, std::string> values;
    std::size_t header_index = lines.size();
    std::vector<std::string> header_words;
    const std::unordered_set<std::string> supported_bare{"AM1", "1SCF", "GEO-OK", "SM5.42R", "DEV"};
    const std::unordered_set<std::string> supported_values{
        "CHARGE", "TLIMIT", "SOLVNT", "IOFR", "ALPHA", "BETA", "GAMMA", "DIELEC", "FACARB", "FEHALO"};

    for (std::size_t index = 0; index < lines.size(); ++index) {
        auto normalized = trim(lines[index]);
        if (normalized.empty()) {
            continue;
        }
        const auto words = split_words(normalized);
        if (molecule_header(words)) {
            header_index = index;
            header_words = words;
            break;
        }
        if (std::isspace(static_cast<unsigned char>(lines[index].front()))) {
            throw InputError("unexpected indented line before molecule header");
        }
        if (!normalized.empty() && normalized.front() == '&') {
            normalized = trim(std::string_view(normalized).substr(1));
        }
        for (auto token : split_words(normalized)) {
            token = upper(std::move(token));
            const auto separator = token.find('=');
            if (separator == std::string::npos) {
                if (!supported_bare.contains(token)) {
                    throw UnsupportedCalculationError("unsupported AMSOL keyword '" + token + "'");
                }
                if (!bare.insert(token).second) {
                    throw InputError("duplicate keyword " + token);
                }
                continue;
            }
            const std::string key = token.substr(0, separator);
            const std::string value = token.substr(separator + 1);
            if (!supported_values.contains(key)) {
                throw UnsupportedCalculationError("unsupported AMSOL option '" + key + "'");
            }
            if (value.empty() || values.contains(key)) {
                throw InputError(value.empty() ? "missing value for " + key : "duplicate option " + key);
            }
            values.emplace(key, value);
        }
    }
    if (header_index == lines.size()) {
        throw InputError("missing molecule name and atom-count line");
    }
    for (const auto* required : {"AM1", "1SCF", "GEO-OK", "SM5.42R"}) {
        if (!bare.contains(required)) {
            throw UnsupportedCalculationError("missing required audited keyword " + std::string(required));
        }
    }
    if (!values.contains("CHARGE")) {
        throw UnsupportedCalculationError("missing required CHARGE=<integer> option");
    }
    if (!values.contains("TLIMIT") || values.at("TLIMIT") != "15") {
        throw UnsupportedCalculationError("the audited legacy input accepts TLIMIT=15 only");
    }
    if (!values.contains("SOLVNT")) {
        throw UnsupportedCalculationError("missing required SOLVNT option");
    }

    CalculationOptions options;
    options.molecular_charge = parse_integer(values.at("CHARGE"), "CHARGE");
    if (values.at("SOLVNT") == "WATER") {
        options.solvent = Solvent::Water;
        const std::unordered_set<std::string> allowed{"CHARGE", "TLIMIT", "SOLVNT"};
        for (const auto& [key, unused] : values) {
            static_cast<void>(unused);
            if (!allowed.contains(key)) {
                throw UnsupportedCalculationError("unsupported WATER solvent option " + key);
            }
        }
        if (bare.contains("DEV")) {
            throw UnsupportedCalculationError("DEV is not accepted for the audited WATER branch");
        }
    } else if (values.at("SOLVNT") == "GENORG") {
        options.solvent = Solvent::Hexadecane;
        require_close(values, "IOFR", 1.4345);
        require_close(values, "ALPHA", 0.0);
        require_close(values, "BETA", 0.0);
        require_close(values, "GAMMA", 38.93);
        require_close(values, "DIELEC", 2.06);
        require_close(values, "FACARB", 0.0);
        require_close(values, "FEHALO", 0.0);
        if (!bare.contains("DEV")) {
            throw UnsupportedCalculationError("audited GENORG hexadecane input requires DEV");
        }
    } else {
        throw UnsupportedCalculationError("unsupported solvent SOLVNT=" + values.at("SOLVNT"));
    }

    const int declared_count = parse_integer(header_words[1], "atom count");
    if (declared_count <= 0) {
        throw InputError("molecule atom count must be positive");
    }
    if (declared_count > 4096) {
        throw InputError("molecule atom count exceeds the legacy parser limit of 4096");
    }
    std::vector<ZMatrixAtom> zmat;
    zmat.reserve(static_cast<std::size_t>(declared_count));
    for (std::size_t line_index = header_index + 1; line_index < lines.size(); ++line_index) {
        const auto normalized = trim(lines[line_index]);
        if (normalized.empty()) {
            continue;
        }
        const auto words = split_words(normalized);
        if (words.size() != 10) {
            throw InputError("expected 10 fields in MOPAC Z-matrix atom row");
        }
        const int number = atomic_number(words[0]);
        const double bond = parse_double(words[1], "bond length");
        const int bond_opt = parse_integer(words[2], "bond optimization flag");
        const double angle = parse_double(words[3], "bond angle");
        const int angle_opt = parse_integer(words[4], "angle optimization flag");
        const double dihedral = parse_double(words[5], "dihedral angle");
        const int dihedral_opt = parse_integer(words[6], "dihedral optimization flag");
        if ((bond_opt != 0 && bond_opt != 1) || (angle_opt != 0 && angle_opt != 1) ||
            (dihedral_opt != 0 && dihedral_opt != 1)) {
            throw UnsupportedCalculationError("unsupported Z-matrix optimization flag");
        }
        zmat.push_back({number, bond, angle, dihedral,
                        parse_integer(words[7], "bond reference"),
                        parse_integer(words[8], "angle reference"),
                        parse_integer(words[9], "dihedral reference")});
    }
    if (zmat.size() != static_cast<std::size_t>(declared_count)) {
        throw InputError("molecule header atom count does not match parsed atom rows");
    }

    ParsedInput result;
    result.molecule_name = header_words[0];
    result.atoms = zmat_to_cartesian(zmat);
    result.options = options;
    return result;
}

}  // namespace amsolcpp
