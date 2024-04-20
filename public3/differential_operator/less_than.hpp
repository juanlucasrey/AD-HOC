/*
 * AD-HOC, Automatic Differentiation for High Order Calculations
 *
 * This file is part of the AD-HOC distribution
 * (https://github.com/juanlucasrey/adhoc2).
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

#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_LESS_THAN_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_LESS_THAN_HPP

#include "../adhoc.hpp"
#include "../dependency.hpp"
#include "differential_operator.hpp"

namespace adhoc3 {

namespace detail {

template <class Nodes, class Id1, class... Ids1, std::size_t Order1,
          std::size_t... Orders1, std::size_t Power1, std::size_t... Powers1,
          class Id2, class... Ids2, std::size_t Order2, std::size_t... Orders2,
          std::size_t Power2, std::size_t... Powers2>
constexpr auto
less_than_check_first(Nodes nodes,
                      std::tuple<der2::p<Power1, der2::d<Order1, Id1>>,
                                 der2::p<Powers1, der2::d<Orders1, Ids1>>...>
                          in1,
                      std::tuple<der2::p<Power2, der2::d<Order2, Id2>>,
                                 der2::p<Powers2, der2::d<Orders2, Ids2>>...>
                          in2);

template <class Nodes, class... Ids1, std::size_t... Orders1,
          std::size_t... Powers1, class... Ids2, std::size_t... Orders2,
          std::size_t... Powers2>
constexpr auto less_than_check_empty(
    Nodes nodes, std::tuple<der2::p<Powers1, der2::d<Orders1, Ids1>>...> in1,
    std::tuple<der2::p<Powers2, der2::d<Orders2, Ids2>>...> in2) {

    if constexpr (std::is_same_v<decltype(in1), decltype(in2)>) {
        return false;
    } else if constexpr (std::is_same_v<decltype(in1), std::tuple<>>) {
        return true;
    } else if constexpr (std::is_same_v<decltype(in2), std::tuple<>>) {
        return false;
    } else {
        // we can now assume that both ids are different and non empty
        return less_than_check_first(nodes, in1, in2);
    }
}

template <class Nodes, class Id1, class... Ids1, std::size_t Order1,
          std::size_t... Orders1, std::size_t Power1, std::size_t... Powers1,
          class Id2, class... Ids2, std::size_t Order2, std::size_t... Orders2,
          std::size_t Power2, std::size_t... Powers2>
constexpr auto less_than_check_first(
    Nodes nodes,
    std::tuple<der2::p<Power1, der2::d<Order1, Id1>>,
               der2::p<Powers1, der2::d<Orders1, Ids1>>...>
    /* in1 */,
    std::tuple<der2::p<Power2, der2::d<Order2, Id2>>,
               der2::p<Powers2, der2::d<Orders2, Ids2>>...> /* in2 */) {
    constexpr auto idx1 = get_idx_first2<Id1>(nodes);
    constexpr auto idx2 = get_idx_first2<Id2>(nodes);
    if constexpr (idx1 < idx2) {
        return false;
    } else if constexpr (idx1 > idx2) {
        return true;
    } else if constexpr (Order1 < Order2) {
        return true;
    } else if constexpr (Order1 > Order2) {
        return false;
    } else if constexpr (Power1 < Power2) {
        return true;
    } else if constexpr (Power1 > Power2) {
        return false;
    } else {
        static_assert(std::is_same_v<der2::p<Power1, der2::d<Order1, Id1>>,
                                     der2::p<Power2, der2::d<Order2, Id2>>>);
        return less_than_check_empty(
            nodes, std::tuple<der2::p<Powers1, der2::d<Orders1, Ids1>>...>{},
            std::tuple<der2::p<Powers2, der2::d<Orders2, Ids2>>...>{});
    }
}

} // namespace detail

template <class Nodes, class... Ids1, std::size_t... Orders1,
          std::size_t... Powers1, class... Ids2, std::size_t... Orders2,
          std::size_t... Powers2>
constexpr auto
less_than(Nodes nodes,
          std::tuple<der2::p<Powers1, der2::d<Orders1, Ids1>>...> in1,
          std::tuple<der2::p<Powers2, der2::d<Orders2, Ids2>>...> in2) {
    return detail::less_than_check_empty(nodes, in1, in2);
}

template <class Nodes, class Id1, std::size_t Order1, std::size_t Power1,
          class Id2, std::size_t Order2, std::size_t Power2>
constexpr auto less_than(Nodes nodes,
                         der2::p<Power1, der2::d<Order1, Id1>> /* in1 */,
                         der2::p<Power2, der2::d<Order2, Id2>> /* in2 */) {
    constexpr auto idx1 = get_idx_first2<Id1>(nodes);
    constexpr auto idx2 = get_idx_first2<Id2>(nodes);
    if constexpr (idx1 < idx2) {
        return false;
    } else if constexpr (idx1 > idx2) {
        return true;
    } else if constexpr (Order1 < Order2) {
        return true;
    } else if constexpr (Order1 > Order2) {
        return false;
    } else {
        return false;
    }
}

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_LESS_THAN_HPP
