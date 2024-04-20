/*
 * AD-HOC, Automatic Differentiation for High Order Calculations
 *
 * This file is part of the AD-HOC distribution
 * (https://github.com/juanlucasrey/AD-HOC).
 * Copyright (c) 2024 Juan Lucas Rey
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

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

namespace detail {

template <std::size_t Balls, bool IsEven, std::size_t... I>
constexpr auto
BinomialCoefficientsSquare_aux(std::index_sequence<I...> /* in */) {
    if constexpr (IsEven) {
        constexpr auto size = sizeof...(I);
        return std::index_sequence<2 * binomial_coefficient(Balls, I)...,
                                   binomial_coefficient(Balls, size)>{};
    } else {
        return std::index_sequence<2 * binomial_coefficient(Balls, I)...>{};
    }
}

} // namespace detail

template <std::size_t Balls> constexpr auto BinomialCoefficientsSquare() {
    constexpr bool is_even = (Balls / 2) * 2 == Balls;
    constexpr auto seq = std::make_index_sequence<(Balls + 1) / 2>{};
    return detail::BinomialCoefficientsSquare_aux<Balls, is_even>(seq);
}

} // namespace adhoc3

#endif // ADHOC3_PARTITION_BINOMIAL_COEFFICIENT_HPP
