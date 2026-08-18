#include <amsolcpp/coordinates.hpp>

#include <amsolcpp/elements.hpp>
#include <amsolcpp/exceptions.hpp>

#include "performance_diagnostics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <sstream>

namespace amsolcpp {
namespace {

using Point = std::array<double, 3>;

double radians(const double degrees) {
    return degrees * std::numbers::pi_v<double> / 180.0;
}

Point place_general(const ZMatrixAtom& atom, const std::span<const Point> coordinates, const std::size_t index) {
    const double distance = atom.bond_length_angstrom;
    const double theta = radians(atom.bond_angle_degrees);
    const double phi = radians(atom.dihedral_degrees);
    const double cosine_angle = std::cos(theta);

    const std::size_t mb = static_cast<std::size_t>(atom.angle_reference - 1);
    const std::size_t mc = static_cast<std::size_t>(atom.bond_reference - 1);
    const std::size_t ma = static_cast<std::size_t>(atom.dihedral_reference - 1);

    double xb = coordinates[mb][0] - coordinates[mc][0];
    const double yb = coordinates[mb][1] - coordinates[mc][1];
    double zb = coordinates[mb][2] - coordinates[mc][2];
    const double rb = std::sqrt(xb * xb + yb * yb + zb * zb);
    if (!(rb > 0.0)) {
        throw InputError("coincident bond and angle references for Z-matrix atom " + std::to_string(index + 1));
    }
    const double inverse_rb = 1.0 / rb;
    if (std::abs(cosine_angle) >= 0.99999999991) {
        const double scale = distance * inverse_rb * cosine_angle;
        return {coordinates[mc][0] + xb * scale,
                coordinates[mc][1] + yb * scale,
                coordinates[mc][2] + zb * scale};
    }

    double xa = coordinates[ma][0] - coordinates[mc][0];
    const double ya = coordinates[ma][1] - coordinates[mc][1];
    double za = coordinates[ma][2] - coordinates[mc][2];
    double xyb = std::sqrt(xb * xb + yb * yb);
    bool rotate_y_first = false;
    if (xyb <= 0.1) {
        std::swap(xa, za);
        za = -za;
        std::swap(xb, zb);
        zb = -zb;
        xyb = std::sqrt(xb * xb + yb * yb);
        rotate_y_first = true;
    }
    if (!(xyb > 0.0)) {
        throw InputError("degenerate reference axis for Z-matrix atom " + std::to_string(index + 1));
    }

    const double costh = xb / xyb;
    const double sinth = yb / xyb;
    const double xpa = xa * costh + ya * sinth;
    const double ypa = ya * costh - xa * sinth;
    const double sinph = zb * inverse_rb;
    const double cosph = std::sqrt(std::abs(1.0 - sinph * sinph));
    const double xqa = xpa * cosph + za * sinph;
    const double zqa = za * cosph - xpa * sinph;
    static_cast<void>(xqa);
    const double yza = std::sqrt(ypa * ypa + zqa * zqa);
    const double coskh = yza < 2.0e-2 ? 1.0 : ypa / yza;
    const double sinkh = yza < 2.0e-2 ? 0.0 : zqa / yza;

    const double sine_angle = std::sin(theta);
    const double xd = distance * cosine_angle;
    const double yd = distance * sine_angle * std::cos(phi);
    const double zd = -distance * sine_angle * std::sin(phi);
    const double ypd = yd * coskh - zd * sinkh;
    const double zpd = zd * coskh + yd * sinkh;
    const double xpd = xd * cosph - zpd * sinph;
    double zqd = zpd * cosph + xd * sinph;
    double xqd = xpd * costh - ypd * sinth;
    const double yqd = ypd * costh + xpd * sinth;
    if (rotate_y_first) {
        const double old_x = xqd;
        xqd = -zqd;
        zqd = old_x;
    }
    return {xqd + coordinates[mc][0], yqd + coordinates[mc][1], zqd + coordinates[mc][2]};
}

}  // namespace

void validate_cartesian_atoms(const std::span<const Atom> atoms) {
    if (atoms.empty()) {
        throw InputError("a calculation requires at least one atom");
    }
    for (std::size_t index = 0; index < atoms.size(); ++index) {
        const auto& atom = atoms[index];
        if (!is_supported_element(atom.atomic_number)) {
            throw UnsupportedElementError(atom.atomic_number, "AM1 and SM5.42R");
        }
        if (!std::isfinite(atom.x_angstrom) || !std::isfinite(atom.y_angstrom) || !std::isfinite(atom.z_angstrom)) {
            throw InputError("atom " + std::to_string(index + 1) + " has a non-finite coordinate");
        }
        constexpr double coordinate_limit_angstrom = 1.0e6;
        if (std::abs(atom.x_angstrom) > coordinate_limit_angstrom
            || std::abs(atom.y_angstrom) > coordinate_limit_angstrom
            || std::abs(atom.z_angstrom) > coordinate_limit_angstrom) {
            throw InputError("atom " + std::to_string(index + 1)
                + " coordinate exceeds the supported magnitude of 1e6 angstrom");
        }
        for (std::size_t previous = 0; previous < index; ++previous) {
            const double dx = atom.x_angstrom - atoms[previous].x_angstrom;
            const double dy = atom.y_angstrom - atoms[previous].y_angstrom;
            const double dz = atom.z_angstrom - atoms[previous].z_angstrom;
            if (dx * dx + dy * dy + dz * dz < 1.0e-2) {
                throw InputError("atoms " + std::to_string(previous + 1) + " and "
                    + std::to_string(index + 1) + " are separated by less than 0.1 angstrom");
            }
        }
    }
}

std::vector<Atom> zmat_to_cartesian(const std::span<const ZMatrixAtom> atoms) {
    AMSOLCPP_PERF_SCOPE(detail::PerformanceStage::GeometryPreprocessing);
    if (atoms.empty()) {
        throw InputError("cannot convert an empty Z-matrix");
    }
    std::vector<Point> coordinates;
    coordinates.reserve(atoms.size());
    for (std::size_t index = 0; index < atoms.size(); ++index) {
        const auto& atom = atoms[index];
        if (!is_supported_element(atom.atomic_number)) {
            throw UnsupportedElementError(atom.atomic_number, "Z-matrix conversion");
        }
        if (index == 0) {
            if (atom.bond_length_angstrom != 0.0 || atom.bond_angle_degrees != 0.0
                || atom.dihedral_degrees != 0.0 || atom.bond_reference != 0
                || atom.angle_reference != 0 || atom.dihedral_reference != 0) {
                throw InputError("Z-matrix atom 1 must have zero geometry fields and references");
            }
            coordinates.push_back({0.0, 0.0, 0.0});
            continue;
        }
        if (!(atom.bond_length_angstrom > 0.0) || !std::isfinite(atom.bond_length_angstrom)) {
            throw InputError("invalid bond length for Z-matrix atom " + std::to_string(index + 1));
        }
        if (index == 1) {
            if (atom.bond_reference != 0 && atom.bond_reference != 1) {
                throw InputError("Z-matrix atom 2 must reference atom 1");
            }
            if (atom.bond_angle_degrees != 0.0 || atom.dihedral_degrees != 0.0
                || atom.angle_reference != 0 || atom.dihedral_reference != 0) {
                throw InputError("Z-matrix atom 2 must have zero angle, dihedral, and unused references");
            }
            coordinates.push_back({atom.bond_length_angstrom, 0.0, 0.0});
            continue;
        }
        if (!std::isfinite(atom.bond_angle_degrees) || atom.bond_angle_degrees < 0.0 ||
            atom.bond_angle_degrees > 180.0) {
            throw InputError("invalid bond angle for Z-matrix atom " + std::to_string(index + 1));
        }
        if (index == 2) {
            if (atom.bond_reference < 1 || atom.bond_reference > 2
                || atom.angle_reference < 1 || atom.angle_reference > 2
                || atom.bond_reference == atom.angle_reference
                || atom.dihedral_reference != 0
                || atom.dihedral_degrees != 0.0) {
                throw InputError("Z-matrix atom 3 must reference distinct atoms 1 and 2 and use zero dihedral data");
            }
            const double theta = radians(atom.bond_angle_degrees);
            const int na = atom.bond_reference == 0 ? 2 : atom.bond_reference;
            const double x = na == 1 ? coordinates[0][0] + atom.bond_length_angstrom * std::cos(theta)
                                      : coordinates[1][0] - atom.bond_length_angstrom * std::cos(theta);
            coordinates.push_back({x, atom.bond_length_angstrom * std::sin(theta), 0.0});
            continue;
        }
        const std::array<int, 3> references{
            atom.bond_reference, atom.angle_reference, atom.dihedral_reference};
        for (const int reference : references) {
            if (reference <= 0 || reference > static_cast<int>(index)) {
                throw InputError("invalid reference for Z-matrix atom " + std::to_string(index + 1));
            }
        }
        if (references[0] == references[1] || references[0] == references[2] || references[1] == references[2]) {
            throw InputError("Z-matrix atom must use three distinct references");
        }
        if (!std::isfinite(atom.dihedral_degrees)) {
            throw InputError("invalid dihedral for Z-matrix atom " + std::to_string(index + 1));
        }
        coordinates.push_back(place_general(atom, coordinates, index));
    }

    std::vector<Atom> result;
    result.reserve(atoms.size());
    for (std::size_t index = 0; index < atoms.size(); ++index) {
        result.push_back({atoms[index].atomic_number, coordinates[index][0], coordinates[index][1], coordinates[index][2]});
    }
    validate_cartesian_atoms(result);
    return result;
}

}  // namespace amsolcpp
