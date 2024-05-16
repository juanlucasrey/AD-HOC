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

#ifndef ADHOC3_PARTITION_MULTINOMIAL_COEFFICIENT_INDEX_SEQUENCE_HPP
#define ADHOC3_PARTITION_MULTINOMIAL_COEFFICIENT_INDEX_SEQUENCE_HPP

#include "multinomial_coefficient.hpp"
#include "tuple_utils.hpp"

#include <array>
#include <cstddef>

namespace adhoc3 {

namespace detail {

template <std::size_t Bins, std::size_t Balls, std::size_t... I>
constexpr auto FirstMultinomialIS_aux(std::index_sequence<I...> /* i */) {
    constexpr auto first = FirstMultinomial<Bins, Balls>();
    return std::index_sequence<first[I]...>{};
}

} // namespace detail

template <std::size_t Bins, std::size_t Balls>
constexpr auto FirstMultinomialIS() {
    return detail::FirstMultinomialIS_aux<Bins, Balls>(
        std::make_index_sequence<Bins>{});
}

namespace detail {

template <std::size_t Bins, std::size_t Balls, std::size_t... I>
constexpr auto LastMultinomialIS_aux(std::index_sequence<I...> /* i */) {
    constexpr auto last = LastMultinomial<Bins, Balls>();
    return std::index_sequence<last[I]...>{};
}

} // namespace detail

template <std::size_t Bins, std::size_t Balls>
constexpr auto LastMultinomialIS() {
    return detail::LastMultinomialIS_aux<Bins, Balls>(
        std::make_index_sequence<Bins>{});
}

namespace detail {

template <std::size_t... I, std::size_t... I2>
constexpr auto NextMultinomialIS_aux(std::index_sequence<I...> /* i */,
                                     std::index_sequence<I2...> /* i2 */) {
    constexpr std::array<std::size_t, sizeof...(I)> temp{I...};
    constexpr auto next = NextMultinomial(temp);
    return std::index_sequence<next[I2]...>{};
}

} // namespace detail

template <std::size_t... I>
constexpr auto NextMultinomialIS(std::index_sequence<I...> i) {
    return detail::NextMultinomialIS_aux(
        i, std::make_index_sequence<sizeof...(I)>{});
}

namespace detail {

template <class CurrentIS, class Output>
constexpr auto MultinomialSequences_aux(CurrentIS const current, Output out) {
    auto constexpr next = NextMultinomialIS(current);
    if constexpr (!std::is_same_v<decltype(current), decltype(next)>) {
        return MultinomialSequences_aux(
            next, std::tuple_cat(out, std::make_tuple(current)));
    } else {
        return std::tuple_cat(out, std::make_tuple(current));
    }
}

} // namespace detail

template <std::size_t Bins, std::size_t Balls>
constexpr auto MultinomialSequences() {
    constexpr auto first = FirstMultinomialIS<Bins, Balls>();
    return detail::MultinomialSequences_aux(first, std::tuple<>{});
}

template <std::size_t Order, std::size_t... I, std::size_t FactorialInputs>
constexpr auto
MultinomialCoeff(std::array<std::size_t, FactorialInputs> const &factorials,
                 std::index_sequence<I...> /* i */) -> std::size_t {
    constexpr std::array<std::size_t, sizeof...(I)> temp{I...};
    return MultinomialCoeff<Order>(factorials, temp);
}

template <std::size_t... I>
constexpr auto MultinomialCoeff2(std::index_sequence<I...> idx_seq)
    -> std::size_t {
    constexpr auto order = sum(idx_seq);
    constexpr auto bins = idx_seq.size();
    constexpr std::array<std::size_t, bins> temp{I...};

    constexpr auto factorials = factorialarr<order>();

    std::size_t res = factorials[order];

    std::size_t k = 0;
    while (k < bins) {
        res /= factorials[temp[k]];
        k++;
    }
    return res;
}

} // namespace adhoc3

#endif // ADHOC3_PARTITION_MULTINOMIAL_COEFFICIENT_INDEX_SEQUENCE_HPP
