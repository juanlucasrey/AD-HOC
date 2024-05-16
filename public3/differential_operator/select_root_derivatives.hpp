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
