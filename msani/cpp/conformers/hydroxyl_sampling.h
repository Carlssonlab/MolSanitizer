#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace StochasticSampling {

// Stores per-bond choices, never their Cartesian product. Small spaces retain
// the historical std::shuffle permutation and RNG consumption exactly.
class HydroxylCombinationSampler {
 public:
  static constexpr std::uint64_t legacyShuffleLimit = 4096;
  static constexpr std::uint64_t countLimit = std::numeric_limits<long long>::max();

  explicit HydroxylCombinationSampler(std::vector<std::vector<double>> choices)
      : choices_(std::move(choices)) {
    for (const auto& values : choices_) {
      if (values.empty() || values.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument("Hydroxyl choices must be nonempty and fit an integer index");
      }
      count_ = cappedProduct(count_, values.size());
    }
  }

  static std::uint64_t cappedProduct(std::uint64_t left, std::uint64_t right) {
    return right && left > countLimit / right ? countLimit : left * right;
  }

  std::uint64_t count() const { return count_; }

  template <typename Random>
  std::vector<std::vector<double>> sample(std::size_t maximum, Random& random) const {
    const auto wanted = static_cast<std::size_t>(std::min<std::uint64_t>(maximum, count_));
    std::vector<std::vector<double>> result;
    result.reserve(wanted);
    if (!wanted) return result;
    if (count_ <= legacyShuffleLimit) {
      std::vector<std::size_t> indices(static_cast<std::size_t>(count_));
      std::iota(indices.begin(), indices.end(), 0);
      random.shuffle(indices);
      for (std::size_t i = 0; i < wanted; ++i) {
        auto index = indices[i];
        std::vector<double> angles(choices_.size());
        for (std::size_t depth = choices_.size(); depth > 0; --depth) {
          const auto& values = choices_[depth - 1];
          angles[depth - 1] = values[index % values.size()];
          index /= values.size();
        }
        result.push_back(std::move(angles));
      }
    } else {
      // Independent uniform digits plus rejection give uniform distinct
      // Cartesian-product indices, including products too large for uint64_t.
      std::set<std::vector<int>> seen;
      while (result.size() < wanted) {
        std::vector<int> indices;
        indices.reserve(choices_.size());
        for (const auto& values : choices_) {
          indices.push_back(random.randint(0, static_cast<int>(values.size() - 1)));
        }
        if (!seen.insert(indices).second) continue;
        std::vector<double> angles;
        angles.reserve(choices_.size());
        for (std::size_t depth = 0; depth < choices_.size(); ++depth) {
          angles.push_back(choices_[depth][indices[depth]]);
        }
        result.push_back(std::move(angles));
      }
    }
    return result;
  }

 private:
  std::vector<std::vector<double>> choices_;
  std::uint64_t count_ = 1;
};

}  // namespace StochasticSampling
