#ifndef ADHOC3_PARTITION_COMBINATIONS_HPP
#define ADHOC3_PARTITION_COMBINATIONS_HPP

#include "binomial_coefficient.hpp"

#include <cstddef>

namespace adhoc3 {

constexpr inline auto combinations(std::size_t bins, std::size_t balls) noexcept
    -> std::size_t {
    return binomial_coefficient(bins + balls - 1, bins - 1);
}

} // namespace adhoc3

#endif // ADHOC3_PARTITION_COMBINATIONS_HPP
