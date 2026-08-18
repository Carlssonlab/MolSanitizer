#pragma once

#include <cstddef>
#include <span>
#include <vector>

namespace amsolcpp {

class PackedSymmetricMatrix {
public:
    explicit PackedSymmetricMatrix(std::size_t dimension, double initial_value = 0.0);

    [[nodiscard]] std::size_t dimension() const noexcept { return dimension_; }
    [[nodiscard]] std::size_t size() const noexcept { return values_.size(); }
    [[nodiscard]] std::span<double> values() noexcept { return values_; }
    [[nodiscard]] std::span<const double> values() const noexcept { return values_; }

    double& operator()(std::size_t row, std::size_t column);
    [[nodiscard]] double operator()(std::size_t row, std::size_t column) const;

    [[nodiscard]] static std::size_t packed_index(std::size_t row, std::size_t column);
    [[nodiscard]] std::vector<double> dense() const;
    [[nodiscard]] static PackedSymmetricMatrix from_dense(
        std::span<const double> dense,
        std::size_t dimension
    );

private:
    std::size_t dimension_{};
    std::vector<double> values_;
};

struct SymmetricEigensystem {
    std::vector<double> eigenvalues;
    // Column-major eigenvectors. Column k corresponds to eigenvalues[k].
    std::vector<double> eigenvectors;
};

[[nodiscard]] SymmetricEigensystem diagonalize_symmetric(
    std::span<const double> dense,
    std::size_t dimension,
    double tolerance = 1.0e-13,
    std::size_t max_sweeps = 100
);

}  // namespace amsolcpp
