#include <amsolcpp/packed_matrix.hpp>

#include <amsolcpp/exceptions.hpp>

#include "amsol_eigensolver.hpp"
#include "performance_diagnostics.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace amsolcpp {
namespace {

std::size_t checked_square(const std::size_t dimension) {
    if (dimension != 0
        && dimension > std::numeric_limits<std::size_t>::max() / dimension) {
        throw InputError("matrix dimension overflows addressable storage");
    }
    return dimension * dimension;
}

std::size_t checked_triangle(const std::size_t dimension) {
    if (dimension == std::numeric_limits<std::size_t>::max()) {
        throw InputError("packed matrix dimension overflows addressable storage");
    }
    std::size_t left = dimension;
    std::size_t right = dimension + 1;
    if (left % 2 == 0) {
        left /= 2;
    } else {
        right /= 2;
    }
    if (right != 0 && left > std::numeric_limits<std::size_t>::max() / right) {
        throw InputError("packed matrix dimension overflows addressable storage");
    }
    return left * right;
}

constexpr std::size_t scaled_pythag_min_dimension = 48;

double amsol_pythag_iterative(const double a, const double b) noexcept {
    double p = std::max(std::abs(a), std::abs(b));
    if (!std::isfinite(p) || p == 0.0) {
        return p;
    }

    const double ratio = std::min(std::abs(a), std::abs(b)) / p;
    double r = ratio * ratio;
    for (;;) {
        const double t = 4.0 + r;
        if (t == 4.0) {
            return p;
        }
        const double s = r / t;
        const double u = 1.0 + 2.0 * s;
        p *= u;
        const double correction = s / u;
        r *= correction * correction;
    }
}

double amsol_pythag_scaled(const double a, const double b) noexcept {
    const double p = std::max(std::abs(a), std::abs(b));
    if (!std::isfinite(p) || p == 0.0) {
        return p;
    }
    const double ratio = std::min(std::abs(a), std::abs(b)) / p;
    return p * std::sqrt(1.0 + ratio * ratio);
}

// EISPACK TRED2 as shipped with AMSOL 7.1. Z is stored column-major so the
// Fortran Z(row,column) access and the public eigenvector layout coincide.
void amsol_tred2(
    const std::size_t dimension,
    std::vector<double>& diagonal,
    std::vector<double>& subdiagonal,
    std::vector<double>& vectors
) {
    const auto z = [&](const std::size_t row,
                       const std::size_t column) -> double& {
        return vectors[column * dimension + row];
    };

    for (std::size_t column = 0; column < dimension; ++column) {
        diagonal[column] = z(dimension - 1, column);
    }

    if (dimension != 1) {
        // Fortran: DO 300 II=2,N; I=N+2-II (I descends from N to 2).
        for (std::size_t i = dimension - 1; i > 0; --i) {
            double h = 0.0;
            double scale = 0.0;
            if (i >= 2) {
                for (std::size_t k = 0; k < i; ++k) {
                    scale += std::abs(diagonal[k]);
                }
            }

            if (i < 2 || scale == 0.0) {
                subdiagonal[i] = diagonal[i - 1];
                for (std::size_t j = 0; j < i; ++j) {
                    diagonal[j] = z(i - 1, j);
                    z(i, j) = 0.0;
                    z(j, i) = 0.0;
                }
            } else {
                for (std::size_t k = 0; k < i; ++k) {
                    diagonal[k] /= scale;
                    h += diagonal[k] * diagonal[k];
                }

                const double f = diagonal[i - 1];
                const double g = -std::copysign(std::sqrt(h), f);
                subdiagonal[i] = scale * g;
                h -= f * g;
                diagonal[i - 1] = f - g;

                std::fill_n(subdiagonal.begin(), i, 0.0);
                for (std::size_t j = 0; j < i; ++j) {
                    const double dj = diagonal[j];
                    z(j, i) = dj;
                    double product = subdiagonal[j] + z(j, j) * dj;
                    for (std::size_t k = j + 1; k < i; ++k) {
                        product += z(k, j) * diagonal[k];
                        subdiagonal[k] += z(k, j) * dj;
                    }
                    subdiagonal[j] = product;
                }

                double projection = 0.0;
                for (std::size_t j = 0; j < i; ++j) {
                    subdiagonal[j] /= h;
                    projection += subdiagonal[j] * diagonal[j];
                }
                const double half_projection = projection / (h + h);
                for (std::size_t j = 0; j < i; ++j) {
                    subdiagonal[j] -= half_projection * diagonal[j];
                }

                for (std::size_t j = 0; j < i; ++j) {
                    const double dj = diagonal[j];
                    const double ej = subdiagonal[j];
                    for (std::size_t k = j; k < i; ++k) {
                        z(k, j) -= dj * subdiagonal[k]
                            + ej * diagonal[k];
                    }
                    diagonal[j] = z(i - 1, j);
                    z(i, j) = 0.0;
                }
            }
            diagonal[i] = h;
        }
    }

    // Accumulate the Householder transformations.
    for (std::size_t i = 1; i < dimension; ++i) {
        z(dimension - 1, i - 1) = z(i - 1, i - 1);
        z(i - 1, i - 1) = 1.0;
        const double h = diagonal[i];
        if (h != 0.0) {
            for (std::size_t k = 0; k < i; ++k) {
                diagonal[k] = z(k, i) / h;
            }
            for (std::size_t j = 0; j < i; ++j) {
                double product = 0.0;
                for (std::size_t k = 0; k < i; ++k) {
                    product += z(k, i) * z(k, j);
                }
                for (std::size_t k = 0; k < i; ++k) {
                    z(k, j) -= product * diagonal[k];
                }
            }
        }
        for (std::size_t k = 0; k < i; ++k) {
            z(k, i) = 0.0;
        }
    }

    for (std::size_t i = 0; i < dimension; ++i) {
        diagonal[i] = z(dimension - 1, i);
        z(dimension - 1, i) = 0.0;
    }
    z(dimension - 1, dimension - 1) = 1.0;
    subdiagonal[0] = 0.0;
}

// EISPACK TQL2 as shipped with AMSOL 7.1. Returns the one-based index of an
// eigenvalue that did not converge, or zero on success.
template<bool use_scaled_pythag>
std::size_t amsol_tql2(
    const std::size_t dimension,
    std::vector<double>& diagonal,
    std::vector<double>& subdiagonal,
    std::vector<double>& vectors
) noexcept {
    const auto z = [&](const std::size_t row,
                       const std::size_t column) -> double& {
        return vectors[column * dimension + row];
    };

    if (dimension == 1) {
        return 0;
    }
    for (std::size_t i = 1; i < dimension; ++i) {
        subdiagonal[i - 1] = subdiagonal[i];
    }

    double accumulated_shift = 0.0;
    double test_scale = 0.0;
    subdiagonal[dimension - 1] = 0.0;
    for (std::size_t l = 0; l < dimension; ++l) {
        std::size_t iterations = 0;
        test_scale = std::max(
            test_scale,
            std::abs(diagonal[l]) + std::abs(subdiagonal[l]));

        std::size_t m = l;
        for (; m < dimension; ++m) {
            if (test_scale + std::abs(subdiagonal[m]) == test_scale) {
                break;
            }
        }

        while (m != l) {
            if (iterations == 30) {
                return l + 1;
            }
            ++iterations;

            const std::size_t l1 = l + 1;
            const double g = diagonal[l];
            double p = (diagonal[l1] - g) / (2.0 * subdiagonal[l]);
            const double r = use_scaled_pythag
                ? amsol_pythag_scaled(p, 1.0)
                : amsol_pythag_iterative(p, 1.0);
            const double signed_r = std::copysign(r, p);
            diagonal[l] = subdiagonal[l] / (p + signed_r);
            diagonal[l1] = subdiagonal[l] * (p + signed_r);
            const double next_diagonal = diagonal[l1];
            const double shift = g - diagonal[l];
            for (std::size_t i = l + 2; i < dimension; ++i) {
                diagonal[i] -= shift;
            }
            accumulated_shift += shift;

            p = diagonal[m];
            double cosine = 1.0;
            double previous_cosine = cosine;
            const double next_subdiagonal = subdiagonal[l1];
            double sine = 0.0;
            double older_cosine = 0.0;
            double previous_sine = 0.0;
            for (std::size_t i = m; i-- > l;) {
                older_cosine = previous_cosine;
                previous_cosine = cosine;
                previous_sine = sine;
                const double scaled_subdiagonal = cosine * subdiagonal[i];
                const double scaled_p = cosine * p;
                const double rotation_norm = use_scaled_pythag
                    ? amsol_pythag_scaled(p, subdiagonal[i])
                    : amsol_pythag_iterative(p, subdiagonal[i]);
                subdiagonal[i + 1] = sine * rotation_norm;
                sine = subdiagonal[i] / rotation_norm;
                cosine = p / rotation_norm;
                p = cosine * diagonal[i]
                    - sine * scaled_subdiagonal;
                diagonal[i + 1] = scaled_p
                    + sine * (cosine * scaled_subdiagonal
                              + sine * diagonal[i]);

                for (std::size_t row = 0; row < dimension; ++row) {
                    const double next = z(row, i + 1);
                    z(row, i + 1) = sine * z(row, i) + cosine * next;
                    z(row, i) = cosine * z(row, i) - sine * next;
                }
            }

            p = -sine * previous_sine * older_cosine
                * next_subdiagonal * subdiagonal[l] / next_diagonal;
            subdiagonal[l] = sine * p;
            diagonal[l] = cosine * p;
            if (test_scale + std::abs(subdiagonal[l]) == test_scale) {
                m = l;
            }
        }
        diagonal[l] += accumulated_shift;
    }

    // Selection-sort eigenvalues and swap the corresponding columns of Z.
    for (std::size_t i = 0; i + 1 < dimension; ++i) {
        std::size_t lowest = i;
        double value = diagonal[i];
        for (std::size_t j = i + 1; j < dimension; ++j) {
            if (diagonal[j] < value) {
                lowest = j;
                value = diagonal[j];
            }
        }
        if (lowest != i) {
            diagonal[lowest] = diagonal[i];
            diagonal[i] = value;
            for (std::size_t row = 0; row < dimension; ++row) {
                std::swap(z(row, i), z(row, lowest));
            }
        }
    }
    return 0;
}

}  // namespace

PackedSymmetricMatrix::PackedSymmetricMatrix(
    const std::size_t dimension,
    const double initial_value
) : dimension_(dimension), values_(checked_triangle(dimension), initial_value) {}

std::size_t PackedSymmetricMatrix::packed_index(
    const std::size_t row,
    const std::size_t column
) {
    const auto major = std::max(row, column);
    const auto minor = std::min(row, column);
    return major * (major + 1) / 2 + minor;
}

double& PackedSymmetricMatrix::operator()(
    const std::size_t row,
    const std::size_t column
) {
    if (row >= dimension_ || column >= dimension_) {
        throw std::out_of_range("packed symmetric matrix index out of range");
    }
    return values_[packed_index(row, column)];
}

double PackedSymmetricMatrix::operator()(
    const std::size_t row,
    const std::size_t column
) const {
    if (row >= dimension_ || column >= dimension_) {
        throw std::out_of_range("packed symmetric matrix index out of range");
    }
    return values_[packed_index(row, column)];
}

std::vector<double> PackedSymmetricMatrix::dense() const {
    std::vector<double> result(checked_square(dimension_));
    for (std::size_t row = 0; row < dimension_; ++row) {
        for (std::size_t column = 0; column < dimension_; ++column) {
            result[row * dimension_ + column] = (*this)(row, column);
        }
    }
    return result;
}

PackedSymmetricMatrix PackedSymmetricMatrix::from_dense(
    const std::span<const double> dense_values,
    const std::size_t dimension
) {
    if (dense_values.size() != checked_square(dimension)) {
        throw InputError("dense symmetric matrix size does not match its dimension");
    }
    PackedSymmetricMatrix result(dimension);
    for (std::size_t row = 0; row < dimension; ++row) {
        for (std::size_t column = 0; column <= row; ++column) {
            const double lower = dense_values[row * dimension + column];
            const double upper = dense_values[column * dimension + row];
            const double scale = std::max({1.0, std::abs(lower), std::abs(upper)});
            if (!std::isfinite(lower) || !std::isfinite(upper) || std::abs(lower - upper) > 1.0e-12 * scale) {
                throw InputError("matrix is non-finite or not symmetric");
            }
            result(row, column) = 0.5 * (lower + upper);
        }
    }
    return result;
}

SymmetricEigensystem diagonalize_symmetric(
    const std::span<const double> dense_values,
    const std::size_t dimension,
    const double tolerance,
    const std::size_t max_sweeps
) {
    AMSOLCPP_PERF_SCOPE(detail::PerformanceStage::Diagonalization);
    AMSOLCPP_PERF_INCREMENT(full_diagonalizations);
    std::size_t square = 0;
    try {
        square = checked_square(dimension);
    } catch (const InputError& error) {
        throw DiagonalizationError(error.what());
    }
    if (dimension == 0 || dense_values.size() != square) {
        throw DiagonalizationError("invalid symmetric eigensystem dimensions");
    }
    if (!std::isfinite(tolerance) || !(tolerance > 0.0) || max_sweeps == 0) {
        throw DiagonalizationError("invalid symmetric eigensolver controls");
    }

    std::vector<double> packed_matrix(checked_triangle(dimension));
    std::vector<std::size_t> row_offsets(dimension);
    for (std::size_t row = 0; row < dimension; ++row) {
        row_offsets[row] = row * (row + 1) / 2;
        for (std::size_t column = 0; column <= row; ++column) {
            const double lower = dense_values[row * dimension + column];
            const double upper = dense_values[column * dimension + row];
            const double scale = std::max({
                1.0, std::abs(lower), std::abs(upper)});
            if (!std::isfinite(lower) || !std::isfinite(upper)
                || std::abs(lower - upper) > 1.0e-12 * scale) {
                throw DiagonalizationError(
                    "matrix is non-finite or not symmetric");
            }
            packed_matrix[row_offsets[row] + column] =
                0.5 * (lower + upper);
        }
    }
    std::span<double> matrix = packed_matrix;
    // Store the Jacobi scratch transposed so each eigenvector is contiguous.
    // The public eigensystem layout is restored when assembling the result.
    std::vector<double> vectors(square, 0.0);
    for (std::size_t i = 0; i < dimension; ++i) {
        vectors[i * dimension + i] = 1.0;
    }

    bool converged = dimension == 1;
    for (std::size_t sweep = 0; sweep < max_sweeps && !converged; ++sweep) {
        AMSOLCPP_PERF_INCREMENT(jacobi_sweeps);
        double largest = 0.0;
        double diagonal_scale = 1.0;
        for (std::size_t i = 0; i < dimension; ++i) {
            diagonal_scale = std::max(
                diagonal_scale, std::abs(matrix[row_offsets[i] + i]));
            for (std::size_t j = i + 1; j < dimension; ++j) {
                largest = std::max(
                    largest, std::abs(matrix[row_offsets[j] + i]));
            }
        }
        if (largest <= tolerance * diagonal_scale) {
            converged = true;
            break;
        }

        for (std::size_t p = 0; p + 1 < dimension; ++p) {
            for (std::size_t q = p + 1; q < dimension; ++q) {
                const std::size_t pq = row_offsets[q] + p;
                const double apq = matrix[pq];
                if (std::abs(apq) <= tolerance * diagonal_scale) {
                    continue;
                }
                AMSOLCPP_PERF_INCREMENT(jacobi_rotations);
                const std::size_t pp = row_offsets[p] + p;
                const std::size_t qq = row_offsets[q] + q;
                const double app = matrix[pp];
                const double aqq = matrix[qq];
                const double rotation_scale = std::max({
                    std::abs(app), std::abs(aqq), std::abs(apq)});
                const double scaled_apq = apq / rotation_scale;
                const double tau =
                    (aqq / rotation_scale - app / rotation_scale)
                    / (2.0 * scaled_apq);
                const double absolute_tau = std::abs(tau);
                const double tau_norm = absolute_tau < 1.0
                    ? std::sqrt(1.0 + absolute_tau * absolute_tau)
                    : absolute_tau * std::sqrt(
                        1.0 + (1.0 / absolute_tau)
                            * (1.0 / absolute_tau));
                const double tangent = std::copysign(1.0, tau) /
                    (absolute_tau + tau_norm);
                const double cosine = 1.0 / std::sqrt(1.0 + tangent * tangent);
                const double sine = tangent * cosine;

                const auto rotate_elements = [&](const std::size_t kp,
                                                 const std::size_t kq) {
                    const double akp = matrix[kp];
                    const double akq = matrix[kq];
                    const double new_kp = cosine * akp - sine * akq;
                    const double new_kq = sine * akp + cosine * akq;
                    matrix[kp] = new_kp;
                    matrix[kq] = new_kq;
                };
                for (std::size_t k = 0; k < p; ++k) {
                    rotate_elements(row_offsets[p] + k,
                                    row_offsets[q] + k);
                }
                for (std::size_t k = p + 1; k < q; ++k) {
                    rotate_elements(row_offsets[k] + p,
                                    row_offsets[q] + k);
                }
                for (std::size_t k = q + 1; k < dimension; ++k) {
                    rotate_elements(row_offsets[k] + p,
                                    row_offsets[k] + q);
                }
                const double new_app = cosine * cosine * app -
                    2.0 * sine * cosine * apq + sine * sine * aqq;
                const double new_aqq = sine * sine * app +
                    2.0 * sine * cosine * apq + cosine * cosine * aqq;
                if (!std::isfinite(new_app) || !std::isfinite(new_aqq)) {
                    throw DiagonalizationError(
                        "Jacobi rotation produced a non-finite eigenvalue");
                }
                matrix[pp] = new_app;
                matrix[qq] = new_aqq;
                matrix[pq] = 0.0;

                for (std::size_t row = 0; row < dimension; ++row) {
                    const double vkp = vectors[p * dimension + row];
                    const double vkq = vectors[q * dimension + row];
                    vectors[p * dimension + row] = cosine * vkp - sine * vkq;
                    vectors[q * dimension + row] = sine * vkp + cosine * vkq;
                }
            }
        }
    }
    if (!converged) {
        double largest = 0.0;
        double diagonal_scale = 1.0;
        for (std::size_t i = 0; i < dimension; ++i) {
            diagonal_scale = std::max(
                diagonal_scale, std::abs(matrix[row_offsets[i] + i]));
            for (std::size_t j = i + 1; j < dimension; ++j) {
                largest = std::max(
                    largest, std::abs(matrix[row_offsets[j] + i]));
            }
        }
        converged = std::isfinite(largest)
            && std::isfinite(diagonal_scale)
            && largest <= tolerance * diagonal_scale;
        if (!converged) {
            throw DiagonalizationError("deterministic Jacobi eigensolver did not converge");
        }
    }

    std::vector<std::size_t> order(dimension);
    std::iota(order.begin(), order.end(), std::size_t{0});
    // Stable insertion sort avoids std::stable_sort's temporary allocation.
    // The dimension is small and the input is already nearly ordered.
    for (std::size_t i = 1; i < dimension; ++i) {
        const std::size_t item = order[i];
        const double item_value = matrix[row_offsets[item] + item];
        std::size_t position = i;
        while (position > 0
               && item_value
                   < matrix[row_offsets[order[position - 1]]
                            + order[position - 1]]) {
            order[position] = order[position - 1];
            --position;
        }
        order[position] = item;
    }

    SymmetricEigensystem result;
    result.eigenvalues.resize(dimension);
    result.eigenvectors.resize(square);
    for (std::size_t column = 0; column < dimension; ++column) {
        const std::size_t old_column = order[column];
        result.eigenvalues[column] =
            matrix[row_offsets[old_column] + old_column];
        for (std::size_t row = 0; row < dimension; ++row) {
            result.eigenvectors[column * dimension + row] =
                vectors[old_column * dimension + row];
        }
    }
    return result;
}

namespace detail {

static SymmetricEigensystem diagonalize_amsol_symmetric_impl(
    const std::span<const double> dense_values,
    const std::size_t dimension,
    const bool validate_input
) {
    AMSOLCPP_PERF_SCOPE(PerformanceStage::Diagonalization);
    AMSOLCPP_PERF_INCREMENT(full_diagonalizations);

    std::size_t square = 0;
    try {
        square = checked_square(dimension);
    } catch (const InputError& error) {
        throw DiagonalizationError(error.what());
    }
    if (dimension == 0 || dense_values.size() != square) {
        throw DiagonalizationError("invalid symmetric eigensystem dimensions");
    }

    SymmetricEigensystem result;
    result.eigenvalues.resize(dimension);
    result.eigenvectors.assign(square, 0.0);
    const auto z = [&](const std::size_t row,
                       const std::size_t column) -> double& {
        return result.eigenvectors[column * dimension + row];
    };

    // HQRII expands a canonical packed matrix into a full symmetric B before
    // RS. TRED2 only consumes B's lower triangle, so copy that triangle
    // directly.  The checked entry validates both dense halves; the SCF-only
    // entry relies on the symmetry invariant established by its Fock builders.
    for (std::size_t row = 0; row < dimension; ++row) {
        for (std::size_t column = 0; column <= row; ++column) {
            const double lower = dense_values[row * dimension + column];
            if (validate_input) {
                const double upper = dense_values[column * dimension + row];
                const double scale = std::max({
                    1.0, std::abs(lower), std::abs(upper)});
                if (!std::isfinite(lower) || !std::isfinite(upper)
                    || std::abs(lower - upper) > 1.0e-12 * scale) {
                    throw DiagonalizationError(
                        "matrix is non-finite or not symmetric");
                }
            }
            z(row, column) = lower;
        }
    }

    std::vector<double> subdiagonal(dimension);
    {
        AMSOLCPP_PERF_SCOPE(PerformanceStage::EigensolverReduction);
        amsol_tred2(
            dimension, result.eigenvalues, subdiagonal,
            result.eigenvectors);
    }
    if (!std::all_of(
            result.eigenvalues.begin(), result.eigenvalues.end(),
            [](const double value) { return std::isfinite(value); })
        || !std::all_of(
            subdiagonal.begin(), subdiagonal.end(),
            [](const double value) { return std::isfinite(value); })
        || !std::all_of(
            result.eigenvectors.begin(), result.eigenvectors.end(),
            [](const double value) { return std::isfinite(value); })) {
        throw DiagonalizationError(
            "AMSOL EISPACK reduction produced a non-finite result");
    }
    std::size_t error_index = 0;
    {
        AMSOLCPP_PERF_SCOPE(PerformanceStage::EigensolverQl);
        // Preserve AMSOL's exact iterative PYTHAG trajectory for small
        // systems.  At the measured hot dimensions, select the faster
        // overflow-safe scaled norm once outside the QL recurrence so the
        // deepest loop contains no run-time mode branch.
        if (dimension < scaled_pythag_min_dimension) {
            error_index = amsol_tql2<false>(
                dimension, result.eigenvalues, subdiagonal,
                result.eigenvectors);
        } else {
            error_index = amsol_tql2<true>(
                dimension, result.eigenvalues, subdiagonal,
                result.eigenvectors);
        }
    }
    if (error_index != 0) {
        throw DiagonalizationError(
            "AMSOL EISPACK QL eigensolver did not converge");
    }
    if (!std::all_of(
            result.eigenvalues.begin(), result.eigenvalues.end(),
            [](const double value) { return std::isfinite(value); })
        || !std::all_of(
            result.eigenvectors.begin(), result.eigenvectors.end(),
            [](const double value) { return std::isfinite(value); })) {
        throw DiagonalizationError(
            "AMSOL EISPACK eigensolver produced a non-finite result");
    }
    return result;
}

SymmetricEigensystem diagonalize_amsol_symmetric(
    const std::span<const double> dense_values,
    const std::size_t dimension
) {
    return diagonalize_amsol_symmetric_impl(
        dense_values, dimension, true);
}

SymmetricEigensystem diagonalize_amsol_symmetric_scf(
    const std::span<const double> dense_values,
    const std::size_t dimension
) {
#ifdef NDEBUG
    return diagonalize_amsol_symmetric_impl(
        dense_values, dimension, false);
#else
    // Keep internal symmetry mistakes immediately visible in developer builds.
    return diagonalize_amsol_symmetric_impl(
        dense_values, dimension, true);
#endif
}

}  // namespace detail

}  // namespace amsolcpp
