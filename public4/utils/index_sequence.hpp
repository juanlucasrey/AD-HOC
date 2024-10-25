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

#ifndef ADHOC4_UTILS_INDEX_SEQUENCE_HPP
#define ADHOC4_UTILS_INDEX_SEQUENCE_HPP

#include "../combinatorics/pow.hpp"

#include <array>
#include <cstddef>
#include <utility>

namespace adhoc4 {

template <std::size_t Idx, class T, T... I>
constexpr auto get(std::integer_sequence<T, I...> /* i */) {
    return std::array<T, sizeof...(I)>{I...}[Idx];
}

template <std::size_t... I>
constexpr auto sum(std::index_sequence<I...> /* idx_seq */) -> std::size_t {
    return (I + ...);
}

template <std::size_t Val, std::size_t... I>
constexpr auto all_equal_to(std::index_sequence<I...> /* idx_seq */) -> bool {
    return ((I == Val) && ...);
}

namespace detail {

template <std::size_t Idx = 0, std::size_t Size, class IS>
auto inner_product_power_aux(double &res, IS current,
                             std::array<double, Size> const &in) {
    constexpr auto current_pow = get<Idx>(current);
    if constexpr (current_pow > 0) {
        res *= pow<current_pow>(in[Idx]);
    }

    if constexpr ((Idx + 1) < current.size()) {
        inner_product_power_aux<Idx + 1>(res, current, in);
    }
}

} // namespace detail

template <std::size_t Size, class IS>
auto inner_product_power(IS current, std::array<double, Size> const &in) {
    static_assert(current.size() <= Size);
    double res = 1.;
    detail::inner_product_power_aux(res, current, in);
    return res;
}

namespace detail {

template <std::size_t Add, std::size_t... I>
constexpr auto make_index_sequence_aux(std::index_sequence<I...> /* i */) {
    return std::index_sequence<(Add + I)...>{};
}

template <std::size_t Add, std::size_t... I>
constexpr auto
make_index_sequence_inverse_aux(std::index_sequence<I...> /* i */) {
    return std::index_sequence<(Add - I)...>{};
}

} // namespace detail

template <std::size_t From, std::size_t To>
constexpr auto make_index_sequence() {
    if constexpr (From <= To) {
        return detail::make_index_sequence_aux<From>(
            std::make_index_sequence<1 + To - From>{});
    } else {
        return detail::make_index_sequence_inverse_aux<From>(
            std::make_index_sequence<1 + From - To>{});
    }
}

namespace detail {

template <bool Value, std::size_t... I>
constexpr auto make_bool_sequence_aux(std::index_sequence<I...>) {
    constexpr auto return_bool = [](std::size_t /* I */, bool b) { return b; };
    return std::integer_sequence<bool, return_bool(I, Value)...>{};
}

} // namespace detail

template <bool Value, std::size_t Size> constexpr auto make_bool_sequence() {
    return detail::make_bool_sequence_aux<Value>(
        std::make_index_sequence<Size>());
}

namespace detail {

template <std::size_t Idx, bool New, bool... I, std::size_t... II>
constexpr auto replace_aux2(std::integer_sequence<bool, I...> /* i */,
                            std::index_sequence<II...> /* is_last */) {
    std::array<bool, sizeof...(I)> temp{I...};
    temp[Idx] = New;
    return temp;
}

template <std::size_t Idx, bool New, bool... I, std::size_t... II>
constexpr auto replace_aux(std::integer_sequence<bool, I...> i,
                           std::index_sequence<II...> idx) {
    constexpr auto result_array = replace_aux2<Idx, New>(i, idx);
    return std::integer_sequence<bool, result_array[II]...>{};
}

} // namespace detail

template <std::size_t Idx, bool New, bool... I>
constexpr auto replace(std::integer_sequence<bool, I...> i) {
    constexpr auto is_size = sizeof...(I);
    static_assert(Idx < is_size);

    return detail::replace_aux<Idx, New>(i,
                                         std::make_index_sequence<is_size>());
}

template <bool T, std::size_t Offset = 0, bool... I>
constexpr auto find(std::integer_sequence<bool, I...> /* i */) -> std::size_t {
    constexpr auto size_is = sizeof...(I);
    std::array<bool, size_is> temp{I...};

    std::size_t loc = Offset;
    while (loc < size_is && temp[loc] != T) {
        loc++;
    }
    return loc;
}

} // namespace adhoc4

#endif // ADHOC4_UTILS_INDEX_SEQUENCE_HPP
