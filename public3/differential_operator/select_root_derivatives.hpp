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

#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_SELECT_ROOT_DERIVATIVES_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_SELECT_ROOT_DERIVATIVES_HPP

#include "../adhoc.hpp"
#include "differential_operator.hpp"

namespace adhoc3 {

template <typename T> struct is_root_class : std::false_type {};

template <std::size_t... Ns, std::size_t... Orders, std::size_t... Powers>
struct is_root_class<
    std::tuple<der2::p<Powers, der2::d<Orders, double_t<Ns>>>...>>
    : std::true_type {};

template <class T> constexpr bool is_root_class_v = is_root_class<T>::value;

// tuple version
template <class... Ders>
constexpr auto select_root_derivatives(std::tuple<Ders...> /* ids */) {
    return std::tuple_cat(
        std::conditional_t<is_root_class_v<Ders>, std::tuple<Ders>,
                           std::tuple<>>{}...);
}

// tuple version
template <class... Ders>
constexpr auto select_non_root_derivatives(std::tuple<Ders...> /* ids */) {
    return std::tuple_cat(
        std::conditional_t<!is_root_class_v<Ders>, std::tuple<Ders>,
                           std::tuple<>>{}...);
}

// variadic version
template <class... Ders> constexpr auto select_root_derivatives2() {
    return std::tuple_cat(
        std::conditional_t<is_root_class_v<Ders>, std::tuple<Ders>,
                           std::tuple<>>{}...);
}

// variadic version
template <class... Ders> constexpr auto select_non_root_derivatives2() {
    return std::tuple_cat(
        std::conditional_t<!is_root_class_v<Ders>, std::tuple<Ders>,
                           std::tuple<>>{}...);
}

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_SELECT_ROOT_DERIVATIVES_HPP
