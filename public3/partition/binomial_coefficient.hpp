#ifndef ADHOC3_PARTITION_BINOMIAL_COEFFICIENT_HPP
#define ADHOC3_PARTITION_BINOMIAL_COEFFICIENT_HPP

#include <cstddef>
#include <utility>

namespace adhoc3 {

constexpr inline auto binomial_coefficient(std::size_t n,
                                           std::size_t k) noexcept
    -> std::size_t {
    return (k > n) ? 0 : // out of range
               (k == 0 || k == n) ? 1
                                  : // edge
               (k == 1 || k == n - 1) ? n
                                      : // first
               (k + k < n) ?            // recursive:
               (binomial_coefficient(n - 1, k - 1) * n) / k
                           : //  path to k=1   is faster
               (binomial_coefficient(n - 1, k) * n) /
                   (n - k); //  path to k=n-1 is faster
}

namespace detail {

template <std::size_t Balls, std::size_t... I>
constexpr auto BinomialCoefficients_aux(std::index_sequence<I...> /* in */) {
    return std::index_sequence<binomial_coefficient(Balls, I)...>{};
}

} // namespace detail

template <std::size_t Balls> constexpr auto BinomialCoefficients() {
    return detail::BinomialCoefficients_aux<Balls>(
        std::make_index_sequence<Balls + 1>{});
}

} // namespace adhoc3

#endif // ADHOC3_PARTITION_BINOMIAL_COEFFICIENT_HPP
