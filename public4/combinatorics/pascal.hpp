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

#ifndef ADHOC4_COMBINATORICS_PASCAL_HPP
#define ADHOC4_COMBINATORICS_PASCAL_HPP

#include <array>
#include <cmath>

namespace adhoc4 {

namespace detail {

template <std::size_t N>
constexpr auto NextPascalArray(std::array<std::size_t, N> const &prev) {
    std::array<std::size_t, N + 1> arr;
    arr.fill(0);

    arr[0] = prev[1];

    for (std::size_t i = 1; i < (prev.size() - 1); i++) {
        arr[i] = prev[i - 1] * (i - 1) + prev[i + 1] * (i + 1);
    }

    arr.back() = prev.back() * (prev.size() - 1);
    return arr;
}

template <std::size_t... I, std::size_t... I2>
constexpr auto NextPascal_aux(std::index_sequence<I...> /* i */,
                              std::index_sequence<I2...> /* i2 */) {
    constexpr std::array<std::size_t, sizeof...(I)> temp{I...};
    constexpr auto next = NextPascalArray(temp);
    return std::index_sequence<next[I2]...>{};
}

} // namespace detail

template <class IndexSequence> constexpr auto NextPascal(IndexSequence i) {
    return detail::NextPascal_aux(i, std::make_index_sequence<i.size() + 1>{});
}

} // namespace adhoc4

#endif // ADHOC4_COMBINATORICS_PASCAL_HPP
