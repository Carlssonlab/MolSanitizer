#pragma once

#include <span>
#include <vector>

#include <amsolcpp/atom.hpp>

namespace amsolcpp {

struct ZMatrixAtom {
    int atomic_number{};
    double bond_length_angstrom{};
    double bond_angle_degrees{};
    double dihedral_degrees{};
    int bond_reference{};
    int angle_reference{};
    int dihedral_reference{};
};

[[nodiscard]] std::vector<Atom> zmat_to_cartesian(std::span<const ZMatrixAtom> atoms);
void validate_cartesian_atoms(std::span<const Atom> atoms);

}  // namespace amsolcpp
