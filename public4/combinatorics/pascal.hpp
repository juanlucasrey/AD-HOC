/*
 * AD-HOC, Automatic Differentiation for High Order Calculations
 *
 * This file is part of the AD-HOC distribution
 * (https://github.com/juanlucasrey/AD-HOC).
 * Copyright (c) 2024 Juan Lucas Rey
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef ADHOC4_COMBINATORICS_PASCAL_HPP
#define ADHOC4_COMBINATORICS_PASCAL_HPP

#include <array>
#include <cstddef>
#include <utility>

namespace adhoc4 {

namespace detail {

template <bool PrevPlus, bool NextPlus, std::size_t N>
constexpr auto NextPascalArray(std::array<int, N> const &prev) {
    constexpr int next_c = NextPlus ? 1 : -1;
    constexpr int prev_c = PrevPlus ? 1 : -1;

    std::array<int, N + 1> arr{};
    arr[0] = next_c * prev[1];
    for (std::size_t i = 1; i < (prev.size() - 1); i++) {
        arr[i] = next_c * prev[i + 1] * static_cast<int>(i + 1) +
                 prev_c * prev[i - 1] * static_cast<int>(i - 1);
    }
    arr.back() = prev_c * prev.back() * (static_cast<int>(prev.size()) - 1);
    return arr;
}

template <bool PrevPlus, bool NextPlus, int... I, std::size_t... I2>
constexpr auto NextPascal_aux(std::integer_sequence<int, I...> /* i */,
                              std::index_sequence<I2...> /* i2 */) {
    constexpr std::array<int, sizeof...(I)> temp{I...};
    constexpr auto next = NextPascalArray<PrevPlus, NextPlus>(temp);
    return std::integer_sequence<int, next[I2]...>{};
}

} // namespace detail

template <bool PrevPlus, bool NextPlus, int... I>
constexpr auto NextPascal(std::integer_sequence<int, I...> i) {
    return detail::NextPascal_aux<PrevPlus, NextPlus>(
        i, std::make_index_sequence<i.size() + 1>{});
}

} // namespace adhoc4

#endif // ADHOC4_COMBINATORICS_PASCAL_HPP
