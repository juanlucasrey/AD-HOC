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

#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_EXPAND_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_EXPAND_HPP

#include "../adhoc.hpp"
#include "../constants_type.hpp"
#include "../dependency.hpp"
#include "differential_operator.hpp"
#include "expand_multinomial.hpp"
#include "expand_univariate.hpp"
#include "is_sorted.hpp"
#include "merge_ordered.hpp"
#include "select_root_derivatives.hpp"
#include "sorted_mult.hpp"

namespace adhoc3 {

namespace detail {

template <typename T> struct is_constant_class : std::false_type {};

template <constants::ArgType N>
struct is_constant_class<constants::CD<N>> : std::true_type {};

template <constants::ArgType N>
struct is_constant_class<inv_t<constants::CD<N>>> : std::true_type {};

template <class T>
constexpr bool is_constant_class_v = is_constant_class<T>::value;

template <class Nodes, template <class, class> class Bivariate, class Id1,
          class Id2>
constexpr auto order_bivariate(Nodes nodes, Bivariate<Id1, Id2> /* id */) {
    constexpr auto id1_less_than_id2 = static_cast<bool>(
        get_idx_first2<Id1>(nodes) >= get_idx_first2<Id2>(nodes));
    if constexpr (id1_less_than_id2) {
        return Bivariate<Id2, Id1>{};
    } else {
        return Bivariate<Id1, Id2>{};
    }
}

template <class Nodes, template <class> class Univariate, class Id,
          std::size_t Order, std::size_t Power>
constexpr auto
expand_single(Nodes nodes,
              der2::p<Power, der2::d<Order, Univariate<Id>>> /* id */) {
    constexpr auto expansion = expand_univariate<Id, Order>(nodes);
    constexpr auto expansion_power =
        expand_multinomial<Power>(nodes, expansion);
    return expansion_power;
}

// default for sum and subs
template <std::size_t Order, template <class, class> class Bivariate, class Id1,
          class Id2>
constexpr auto expand_bivariate(Bivariate<Id1, Id2> /* id */) {
    return std::tuple<std::tuple<der2::p<1, der2::d<Order, Id1>>>,
                      std::tuple<der2::p<1, der2::d<Order, Id2>>>>{};
}

template <std::size_t Order, class IdFirst, class IdSecond, std::size_t... I>
constexpr auto generate_operators_mul(std::index_sequence<I...> /* i */) {
    return std::make_tuple(
        std::tuple<der2::p<1, der2::d<Order, IdFirst>>>{},
        std::tuple<der2::p<1, der2::d<Order - I - 1, IdFirst>>,
                   der2::p<1, der2::d<I + 1, IdSecond>>>{}...,
        std::tuple<der2::p<1, der2::d<Order, IdSecond>>>{});
}

template <std::size_t Order, class Id1, class Id2>
constexpr auto expand_bivariate(mul_t<Id1, Id2> /* id */) {
    constexpr auto seq = std::make_index_sequence<Order - 1>{};
    return generate_operators_mul<Order, Id1, Id2>(seq);
}

template <class Nodes, template <class, class> class Bivariate, class Id1,
          class Id2, std::size_t Order, std::size_t Power>
constexpr auto
expand_single(Nodes nodes,
              der2::p<Power, der2::d<Order, Bivariate<Id1, Id2>>> /* id */) {
    if constexpr (is_constant_class_v<Id1>) {
        return std::tuple<std::tuple<der2::p<Power, der2::d<Order, Id2>>>>{};
    } else if constexpr (is_constant_class_v<Id2>) {
        return std::tuple<std::tuple<der2::p<Power, der2::d<Order, Id1>>>>{};
    } else {
        constexpr auto ordered_bivariate =
            order_bivariate(nodes, Bivariate<Id1, Id2>{});
        constexpr auto expansion = expand_bivariate<Order>(ordered_bivariate);
        return expand_multinomial<Power>(nodes, expansion);
    }
}

template <class Arg1, class Arg2>
constexpr auto expand_aux(Arg1 arg1, Arg2 arg2);

template <class Nodes>
constexpr auto expand_aux(Nodes /* nodes */, std::tuple<> id) {
    return id;
}

template <class Nodes, std::size_t N, std::size_t Order, std::size_t Power,
          class... RestFirstOp, class... Ops>
constexpr auto expand_aux(
    Nodes /* nodes */,
    std::tuple<
        std::tuple<der2::p<Power, der2::d<Order, double_t<N>>>, RestFirstOp...>,
        Ops...>
        id) {
    return id;
}

template <class Nodes, class FirstOp, class... RestFirstOp, class... Ops>
constexpr auto
expand_aux(Nodes nodes,
           std::tuple<std::tuple<FirstOp, RestFirstOp...>, Ops...> /* id */) {
    constexpr auto single_result = expand_single(nodes, FirstOp{});
    constexpr auto append_rest = multiply_ordered_tuple(
        nodes, single_result, std::tuple<RestFirstOp...>{});
    return expand_aux(nodes,
                      merge_ordered(nodes, append_rest, std::tuple<Ops...>{}));
}

} // namespace detail

template <class Nodes, class DersIn>
constexpr auto expand(Nodes nodes, DersIn dersin) {
    static_assert(std::tuple_size<DersIn>{} > 0, "no output derivatives!");
    static_assert(is_ordered(nodes, dersin));
    return detail::expand_aux(nodes, dersin);
}

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_EXPAND_HPP
