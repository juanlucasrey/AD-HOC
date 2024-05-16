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

#ifndef ADHOC_UTILS_HPP
#define ADHOC_UTILS_HPP

#include <type_traits>

namespace adhoc {

namespace detail {
template <typename T, typename... Ts> struct get_index;

template <typename T, typename... Ts>
struct get_index<T, T, Ts...> : std::integral_constant<std::size_t, 0> {};

template <typename T, typename... Ts>
struct get_index<T, const T, Ts...> : std::integral_constant<std::size_t, 0> {};

template <typename T, typename... Ts>
struct get_index<const T, T, Ts...> : std::integral_constant<std::size_t, 0> {};

template <typename T, typename Tail, typename... Ts>
struct get_index<T, Tail, Ts...>
    : std::integral_constant<std::size_t, 1 + get_index<T, Ts...>::value> {};
} // namespace detail

template <typename T, typename... Ts> constexpr auto has_type() -> bool {
    return (std::is_same_v<T, Ts> || ...) ||
           (std::is_same_v<T, Ts const> || ...);
};

template <typename T, typename... Ts> constexpr auto idx_type() -> std::size_t {
    static_assert(has_type<T, Ts...>(), "variable not on this tape");
    return detail::get_index<T, Ts...>::value;
};

template <typename T, typename... Ts>
constexpr auto idx_type2() -> std::size_t {
    if constexpr (has_type<T, Ts...>()) {
        return detail::get_index<T, Ts...>::value;
    } else {
        return sizeof...(Ts);
    }
};

} // namespace adhoc

#endif // ADHOC_UTILS_HPP
