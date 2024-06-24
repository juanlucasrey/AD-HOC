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

#ifndef ADHOC4_UTILS_TUPLE_HPP
#define ADHOC4_UTILS_TUPLE_HPP

#include <array>
#include <cmath>
#include <concepts>

namespace adhoc4 {

namespace detail {

template <std::size_t Idx = 0, class T>
auto constexpr get_first_type_idx_aux(std::tuple<> /* tuple */, T /* value */)
    -> std::size_t {
    return Idx;
}

template <std::size_t Idx = 0, class FirstT, class... Ts, class T>
auto constexpr get_first_type_idx_aux(std::tuple<FirstT, Ts...> /* tuple */,
                                      T /* value */) -> std::size_t {
    if constexpr (std::same_as<FirstT, T>) {
        return Idx;
    } else {
        return get_first_type_idx_aux<Idx + 1>(std::tuple<Ts...>{}, T{});
    }
}

} // namespace detail

template <class... Ts, class T>
auto constexpr get_first_type_idx(std::tuple<Ts...> /* tuple */, T /* value */)
    -> std::size_t {
    return detail::get_first_type_idx_aux(std::tuple<Ts...>{}, T{});
}

} // namespace adhoc4

#endif // ADHOC4_UTILS_TUPLE_HPP
