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

#ifndef ADHOC3_UTILS_HPP
#define ADHOC3_UTILS_HPP

#include <type_traits>

namespace adhoc3 {

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

} // namespace adhoc3

#endif // ADHOC3_UTILS_HPP
