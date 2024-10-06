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

#ifndef ADHOC4_UTILS_BIVARIATE_HPP
#define ADHOC4_UTILS_BIVARIATE_HPP

#include "../adhoc.hpp"
#include "../combinatorics/factorial.hpp"
#include "../differential_operator.hpp"
#include "../utils/tuple.hpp"

#include <cstddef>

namespace adhoc4 {

namespace detail {

template <template <class, class> class Bivariate, class Id1, class Id2,
          class Nodes>
constexpr auto sort_bivariate(Bivariate<Id1, Id2> id, Nodes nodes) {
    constexpr auto id1_less_than_id2 = static_cast<bool>(
        get_first_type_idx(nodes, Id1{}) >= get_first_type_idx(nodes, Id2{}));
    if constexpr (id1_less_than_id2) {
        return Bivariate<Id2, Id1>{};
    } else {
        return id;
    }
}

// default for sum and subs
template <std::size_t Order, class IdFirst, class IdSecond, std::size_t... I>
constexpr auto expand_bivariate_dflt(std::index_sequence<I...> /* i */) {
    return std::make_tuple(std::tuple<der::d<I + 1, IdFirst>,
                                      der::d<Order - I - 1, IdSecond>>{}...,
                           std::tuple<der::d<Order, IdFirst>>{},
                           std::tuple<der::d<Order, IdSecond>>{});
}

// default for sum and subs
template <std::size_t Order, template <class, class> class Bivariate, class Id1,
          class Id2>
constexpr auto expand_bivariate_aux(Bivariate<Id1, Id2> /* id */) {
    constexpr auto seq = std::make_index_sequence<Order - 1>{};
    return expand_bivariate_dflt<Order, Id1, Id2>(seq);
}

template <std::size_t Order, std::size_t Order1, class Id1, class Id2,
          std::size_t... I>
constexpr auto expand_bivariate_mul3(std::index_sequence<I...> /* i */) {
    return std::make_tuple(
        std::tuple<der::d<Order1, Id1>, der::d<I + 1, Id2>>{}...,
        std::tuple<der::d<Order1, Id1>>{});
}

template <std::size_t Order, std::size_t Order1, class Id1, class Id2,
          std::size_t... I>
constexpr auto expand_bivariate_mul2(std::index_sequence<I...> /* i */) {
    if constexpr (Order1) {
        constexpr std::size_t add = Order - Order1;
        if constexpr (add) {
            return std::make_tuple(
                std::tuple<der::d<Order1, Id1>, der::d<I + add, Id2>>{}...);
        } else {
            return expand_bivariate_mul3<Order, Order1, Id1, Id2>(
                std::make_index_sequence<sizeof...(I) - 1>{});
        }

    } else {
        return std::make_tuple(std::tuple<der::d<I + Order, Id2>>{}...);
    }
}

template <std::size_t Order, class Id1, class Id2, std::size_t... I>
constexpr auto expand_bivariate_mul(std::index_sequence<I...> /* i */) {
    return std::tuple_cat(expand_bivariate_mul2<Order, I + 1, Id1, Id2>(
                              std::make_index_sequence<I + 2>{})...,
                          expand_bivariate_mul2<Order, 0, Id1, Id2>(
                              std::make_index_sequence<1>{}));
}

template <std::size_t Order, class Id1, class Id2>
constexpr auto expand_bivariate_aux(mul_t<Id1, Id2> /* id */) {
    constexpr auto seq = std::make_index_sequence<Order>{};
    return expand_bivariate_mul<Order, Id1, Id2>(seq);
}

} // namespace detail

template <template <class, class> class Bivariate, class Id1, class Id2,
          std::size_t Order, class Nodes>
constexpr auto expand_bivariate(der::d<Order, Bivariate<Id1, Id2>>
                                /* id */,
                                Nodes nodes) {
    if constexpr (is_constant(Id1{})) {
        return std::tuple<std::tuple<der::d<Order, Id2>>>{};
    } else if constexpr (is_constant(Id2{})) {
        return std::tuple<std::tuple<der::d<Order, Id1>>>{};
    } /* else if constexpr (std::is_same_v<Id1, Id2>) {
        // todo
    } */
    else {
        constexpr auto ordered_bivariate =
            detail::sort_bivariate(Bivariate<Id1, Id2>{}, nodes);
        return detail::expand_bivariate_aux<Order>(ordered_bivariate);
    }
}

} // namespace adhoc4

#endif // ADHOC4_UTILS_BIVARIATE_HPP
