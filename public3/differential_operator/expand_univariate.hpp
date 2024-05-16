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

#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_EXPAND_UNIVARIATE_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_EXPAND_UNIVARIATE_HPP

#include "../partition/integer_partition_index_sequence.hpp"
#include "../partition/partition_function.hpp"
#include "sorted_mult.hpp"

#include <cstddef>
#include <tuple>
#include <utility>

namespace adhoc3 {

namespace detail {

template <std::size_t Idx, std::size_t End, class Nodes, class... Ids,
          std::size_t... I, class... Out>
constexpr auto expand_univariate_aux(Nodes nodes, std::tuple<Ids...> in,
                                     std::index_sequence<I...> prev,
                                     std::tuple<Out...> out) {
    auto constexpr curent = multiply_ordered(nodes, pow<I>(Ids{})...);
    if constexpr (Idx == End) {
        return std::tuple_cat(out, std::make_tuple(curent));
    } else {
        auto constexpr next = NextPartitionIS(prev);
        return expand_univariate_aux<Idx + 1, End>(
            nodes, in, next, std::tuple_cat(out, std::make_tuple(curent)));
    }
};

template <class Id, std::size_t Order, std::size_t... I>
constexpr auto expand_univariate_initial(std::index_sequence<I...> /* i */) {
    return std::tuple<std::tuple<der2::p<1, der2::d<I + 1, Id>>>...>{};
}

} // namespace detail

template <class Id, std::size_t Order, class Nodes>
constexpr auto expand_univariate(Nodes nodes) {
    constexpr auto in = detail::expand_univariate_initial<Id, Order>(
        std::make_index_sequence<Order>{});
    constexpr auto first = FirstPartitionIS<Order>();
    return detail::expand_univariate_aux<1, partition_function(Order)>(
        nodes, in, first, std::tuple<>{});
}

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_EXPAND_UNIVARIATE_HPP
