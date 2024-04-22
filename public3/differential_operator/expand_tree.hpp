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

#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_EXPAND_TREE_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_EXPAND_TREE_HPP

#include "../adhoc.hpp"
#include "../constants_type.hpp"
#include "../dependency.hpp"
#include "differential_operator.hpp"
#include "expand_multinomial.hpp"
#include "expand_univariate.hpp"
#include "merge_ordered.hpp"
#include "select_root_derivatives.hpp"
#include "tuple_utils.hpp"
#include <differential_operator/is_sorted.hpp>
#include <differential_operator/sorted_mult.hpp>

namespace adhoc3 {

namespace detail {

template <typename T> struct is_constant_class2 : std::false_type {};

template <constants::ArgType N>
struct is_constant_class2<constants::CD<N>> : std::true_type {};

template <constants::ArgType N>
struct is_constant_class2<inv_t<constants::CD<N>>> : std::true_type {};

template <class T>
constexpr bool is_constant_class2_v = is_constant_class2<T>::value;

template <class Nodes, template <class, class> class Bivariate, class Id1,
          class Id2>
constexpr auto order_bivariate2(Nodes nodes, Bivariate<Id1, Id2> /* id */) {
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
expand_tree_single(Nodes nodes,
                   der2::p<Power, der2::d<Order, Univariate<Id>>> /* id */) {
    constexpr auto expansion = expand_univariate<Id, Order>(nodes);
    constexpr auto expansion_power =
        expand_multinomial<Power>(nodes, expansion);
    return expansion_power;
}

// default for sum and subs
template <std::size_t Order, template <class, class> class Bivariate, class Id1,
          class Id2>
constexpr auto expand_tree_bivariate(Bivariate<Id1, Id2> /* id */) {
    return std::tuple<std::tuple<der2::p<1, der2::d<Order, Id1>>>,
                      std::tuple<der2::p<1, der2::d<Order, Id2>>>>{};
}

template <std::size_t Order, class IdFirst, class IdSecond, std::size_t... I>
constexpr auto generate_operators_mul2(std::index_sequence<I...> /* i */) {
    return std::make_tuple(
        std::tuple<der2::p<1, der2::d<Order, IdFirst>>>{},
        std::tuple<der2::p<1, der2::d<Order - I - 1, IdFirst>>,
                   der2::p<1, der2::d<I + 1, IdSecond>>>{}...,
        std::tuple<der2::p<1, der2::d<Order, IdSecond>>>{});
}

template <std::size_t Order, class Id1, class Id2>
constexpr auto expand_tree_bivariate(mul_t<Id1, Id2> /* id */) {
    constexpr auto seq = std::make_index_sequence<Order - 1>{};
    return generate_operators_mul2<Order, Id1, Id2>(seq);
}

// default for sum and subs
template <std::size_t Order, template <class, class> class Bivariate, class Id>
constexpr auto expand_tree_bivariate_same(Bivariate<Id, Id> /* id */) {
    return std::tuple<std::tuple<der2::p<1, der2::d<Order, Id>>>>{};
}

template <std::size_t Order, class Id, std::size_t... I>
constexpr auto generate_operators_mul_same(std::index_sequence<I...> /* i */) {
    if constexpr ((Order / 2) * 2 == Order) {
        return std::make_tuple(
            std::tuple<der2::p<1, der2::d<Order, Id>>>{},
            std::tuple<der2::p<1, der2::d<Order - I - 1, Id>>,
                       der2::p<1, der2::d<I + 1, Id>>>{}...,
            std::tuple<der2::p<2, der2::d<Order / 2, Id>>>{});
    } else {
        return std::make_tuple(
            std::tuple<der2::p<1, der2::d<Order, Id>>>{},
            std::tuple<der2::p<1, der2::d<Order - I - 1, Id>>,
                       der2::p<1, der2::d<I + 1, Id>>>{}...);
    }
}

template <std::size_t Order, class Id>
constexpr auto expand_tree_bivariate_same(mul_t<Id, Id> /* id */) {
    constexpr auto seq = std::make_index_sequence<(Order - 1) / 2>{};
    return generate_operators_mul_same<Order, Id>(seq);
}

template <class Nodes, template <class, class> class Bivariate, class Id1,
          class Id2, std::size_t Order, std::size_t Power>
constexpr auto expand_tree_single(
    Nodes nodes, der2::p<Power, der2::d<Order, Bivariate<Id1, Id2>>> /* id */) {
    if constexpr (is_constant_class2_v<Id1>) {
        return std::tuple<std::tuple<der2::p<Power, der2::d<Order, Id2>>>>{};
    } else if constexpr (is_constant_class2_v<Id2>) {
        return std::tuple<std::tuple<der2::p<Power, der2::d<Order, Id1>>>>{};
    } else if constexpr (std::is_same_v<Id1, Id2>) {
        constexpr auto expansion =
            expand_tree_bivariate_same<Order>(Bivariate<Id1, Id2>{});
        return expand_multinomial<Power>(nodes, expansion);
    } else {
        constexpr auto ordered_bivariate =
            order_bivariate2(nodes, Bivariate<Id1, Id2>{});
        constexpr auto expansion =
            expand_tree_bivariate<Order>(ordered_bivariate);
        return expand_multinomial<Power>(nodes, expansion);
    }
}

template <class Arg1, class Arg2>
constexpr auto expand_tree_aux(Arg1 arg1, Arg2 arg2);

template <class Nodes>
constexpr auto expand_tree_aux(Nodes /* nodes */, std::tuple<> id) {
    return id;
}

template <class Nodes, class DersOut, std::size_t N, std::size_t Order,
          std::size_t Power, class... RestFirstOp, class... Ops>
constexpr auto expand_tree_aux(
    Nodes /* nodes */, DersOut dersout,
    std::tuple<
        std::tuple<der2::p<Power, der2::d<Order, double_t<N>>>, RestFirstOp...>,
        Ops...>
        id) {
    constexpr auto all_included = std::apply(
        [id](auto... single_output) {
            return (contains(id, single_output) && ...);
        },
        dersout);
    static_assert(all_included,
                  "output derivative not reached by input derivatives");
    return dersout;
}

template <class Nodes, class FirstOp, class... RestFirstOp>
constexpr auto expand_single(Nodes nodes,
                             std::tuple<FirstOp, RestFirstOp...> /* id */) {
    constexpr auto single_result = expand_tree_single(nodes, FirstOp{});
    return multiply_ordered_tuple(nodes, single_result,
                                  std::tuple<RestFirstOp...>{});
}

template <class Nodes, class DersOut, class FirstOp, class... RestFirstOp,
          class... Ops>
constexpr auto expand_tree_aux(
    Nodes nodes, DersOut dersout,
    std::tuple<std::tuple<FirstOp, RestFirstOp...>, Ops...> /* id */) {
    constexpr auto append_rest =
        expand_single(nodes, std::tuple<FirstOp, RestFirstOp...>{});
    constexpr auto alive_nodes = expand_tree_aux(
        nodes, dersout,
        merge_ordered(nodes, append_rest, std::tuple<Ops...>{}));

    constexpr auto any_included = std::apply(
        [alive_nodes](auto... single_output) {
            return (contains(alive_nodes, single_output) || ...);
        },
        append_rest);

    if constexpr (any_included) {
        return std::tuple_cat(
            std::make_tuple(std::tuple<FirstOp, RestFirstOp...>{}),
            alive_nodes);
    } else {
        return alive_nodes;
    }
}

} // namespace detail

template <class Nodes, class DersOut, class DersIn>
constexpr auto expand_tree(Nodes nodes, DersOut dersout, DersIn dersin) {
    static_assert(std::tuple_size<DersIn>{} > 0, "no output derivatives!");
    static_assert(is_ordered(nodes, dersin));
    constexpr auto all_nodes = detail::expand_tree_aux(nodes, dersout, dersin);

    constexpr auto all_included = std::apply(
        [all_nodes](auto... single_input) {
            return (contains(all_nodes, single_input) && ...);
        },
        dersin);
    static_assert(all_included,
                  "input derivative doesn't reach any output derivatives");

    return select_non_root_derivatives(all_nodes);
}

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_EXPAND_TREE_HPP
