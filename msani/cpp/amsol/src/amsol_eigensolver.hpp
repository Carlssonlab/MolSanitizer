#pragma once

#include <amsolcpp/packed_matrix.hpp>

#include <cstddef>
#include <span>

namespace amsolcpp::detail {

// Internal port of the AMSOL 7.1 HQRII -> RS -> TRED2/TQL2 path.
// The input is dense row-major; returned eigenvectors use the public
// column-major layout (column * dimension + row).
[[nodiscard]] SymmetricEigensystem diagonalize_amsol_symmetric(
    std::span<const double> dense,
    std::size_t dimension
);

// SCF constructs and updates both halves of every Fock matrix together.  This
// entry point relies on that internal invariant and avoids validating the
// unused upper triangle on every iteration.  Other callers use the checked
// function above.
[[nodiscard]] SymmetricEigensystem diagonalize_amsol_symmetric_scf(
    std::span<const double> dense,
    std::size_t dimension
);

}  // namespace amsolcpp::detail
