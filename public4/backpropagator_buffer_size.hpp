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

#ifndef ADHOC4_BACKPROPAGATOR_BUFFER_SIZE_HPP
#define ADHOC4_BACKPROPAGATOR_BUFFER_SIZE_HPP

#include "backpropagator_tools.hpp"
#include "combinatorics/trinomial.hpp"
#include "dependency.hpp"
#include "differential_operator.hpp"
#include "monomial_included.hpp"
#include "utils/tuple.hpp"
#include "utils/univariate.hpp"

#include <cstddef>
#include <tuple>
#include <type_traits>

#define LOG false

#if LOG
#include "../tests4/type_name.hpp"
#include <iostream>
#endif

namespace adhoc4::detail {

template <std::size_t Last, class CalcTree, std::size_t SizeMax,
          std::size_t SizeCur, std::size_t MaxOrder, std::size_t N = 0,
          class PrimalSubNodeOrdered1, class PrimalSubNodeOrdered2,
          class DerivativeNodes, class DerivativeNodeNew,
          class DerivativeNodeInputs>
constexpr auto treat_nodes_mul_buffer_size(
    std::tuple<PrimalSubNodeOrdered1, PrimalSubNodeOrdered2> pn,
    DerivativeNodes dn, DerivativeNodeNew dnn, DerivativeNodeInputs dnin) {

    if constexpr (N == Last) {
        return std::make_tuple(dnn,
                               std::integral_constant<std::size_t, SizeMax>{},
                               std::integral_constant<std::size_t, SizeCur>{});
    } else {
        constexpr auto current_node_der = std::get<0>(std::get<N>(dn));
        constexpr auto current_node_loc = std::get<1>(std::get<N>(dn));

#if LOG
        std::cout << "treating derivative tree node" << std::endl;
        std::cout << type_name2<decltype(current_node_der)>() << std::endl;
#endif

        constexpr auto current_derivative_subnode_rest = tail(current_node_der);
        constexpr auto pow = get_power(head(current_node_der));
        constexpr auto rest_power = power(current_derivative_subnode_rest);

        constexpr auto multinomial_sequences =
            TrinomialSequencesMult<pow, MaxOrder - rest_power>();

        using NodesValue = CalcTree::ValuesTupleInverse;
        constexpr auto diff_ops = std::make_tuple(
            d(PrimalSubNodeOrdered1{}) * d(PrimalSubNodeOrdered2{}),
            d(PrimalSubNodeOrdered1{}), d(PrimalSubNodeOrdered2{}));

        constexpr auto next_derivatives = std::apply(
            [current_derivative_subnode_rest, diff_ops](const auto... seq) {
                return std::make_tuple(multiply_ordered(
                    create_differential_operator(diff_ops, seq),
                    current_derivative_subnode_rest, NodesValue{})...);
            },
            multinomial_sequences);

        constexpr auto flags_next_derivatives = std::apply(
            [dnin](auto... next_derivative) {
                return std::integer_sequence<
                    bool, monomial_included(next_derivative, dnin)...>{};
            },
            next_derivatives);

        constexpr auto next_derivatives_filtered =
            filter(next_derivatives, flags_next_derivatives);

        constexpr auto multinomial_sequences_filtered =
            filter(multinomial_sequences, flags_next_derivatives);

        static_assert(size(next_derivatives_filtered));

        constexpr auto locations =
            locate_new_vals_update_buffer_types_size<Last>(
                next_derivatives_filtered, dn, dnn, dnin);

        constexpr auto flags_only_new = std::apply(
            [](auto... location) {
                return std::integer_sequence<
                    bool,
                    std::is_same_v<on_buffer_t, decltype(location.first)>...>{};
            },
            locations);

        constexpr auto next_derivatives_new =
            filter(next_derivatives_filtered, flags_only_new);

        constexpr auto locations_new = filter(locations, flags_only_new);
        constexpr auto next_derivatives_new_with_pos =
            add_position(next_derivatives_new, locations_new);
        constexpr auto dnn_new =
            merge_sorted(next_derivatives_new_with_pos, dnn, NodesValue{});

        constexpr std::size_t size_cur_new =
            SizeCur + size(next_derivatives_new_with_pos) -
            static_cast<std::size_t>(free_on_buffer_size(current_node_loc));
        constexpr std::size_t size_max_new = std::max(SizeMax, size_cur_new);

        return treat_nodes_mul_buffer_size<Last, CalcTree, size_max_new,
                                           size_cur_new, MaxOrder, N + 1>(
            pn, dn, dnn_new, dnin);
    }
}

template <std::size_t Last, class CalcTree, std::size_t SizeMax,
          std::size_t SizeCur, class PrimalSubNode1, class PrimalSubNode2,
          class DerivativeNodes, class DerivativeNodeInputs>
constexpr auto treat_nodes_specialized_buffer_size(
    mul_t<PrimalSubNode1, PrimalSubNode2> /* pn */, DerivativeNodes dn,
    DerivativeNodeInputs dnin) {

    using NodesValue = CalcTree::ValuesTupleInverse;
    constexpr auto ordered_pair = detail::sort_pair(
        std::tuple<PrimalSubNode1, PrimalSubNode2>{}, NodesValue{});

    constexpr auto MaxOrder = detail::max_orders(dnin);

    return treat_nodes_mul_buffer_size<Last, CalcTree, SizeMax, SizeCur,
                                       MaxOrder>(ordered_pair, dn,
                                                 std::tuple<>{}, dnin);
}

template <std::size_t Last, class CalcTree, std::size_t SizeMax,
          std::size_t SizeCur, std::size_t N = 0, class PrimalSubNodeOrdered1,
          class PrimalSubNodeOrdered2, class DerivativeNodes,
          class DerivativeNodeNew, class DerivativeNodeInputs>
constexpr auto treat_nodes_add_buffer_size(
    std::tuple<PrimalSubNodeOrdered1, PrimalSubNodeOrdered2> pn,
    DerivativeNodes dn, DerivativeNodeNew dnn, DerivativeNodeInputs dnin) {

    if constexpr (N == Last) {
        return std::make_tuple(dnn,
                               std::integral_constant<std::size_t, SizeMax>{},
                               std::integral_constant<std::size_t, SizeCur>{});

    } else {
        constexpr auto current_node_der = std::get<0>(std::get<N>(dn));
        constexpr auto current_node_loc = std::get<1>(std::get<N>(dn));

#if LOG
        std::cout << "treating derivative tree node" << std::endl;
        std::cout << type_name2<decltype(current_node_der)>() << std::endl;
#endif

        constexpr auto current_derivative_subnode_rest = tail(current_node_der);
        constexpr auto pow = get_power(head(current_node_der));

        constexpr auto multinomial_sequences = MultinomialSequences<2, pow>();
        constexpr auto diff_ops = std::make_tuple(d(PrimalSubNodeOrdered1{}),
                                                  d(PrimalSubNodeOrdered2{}));

        using NodesValue = CalcTree::ValuesTupleInverse;

        constexpr auto next_derivatives = std::apply(
            [current_derivative_subnode_rest, diff_ops](const auto... seq) {
                return std::make_tuple(multiply_ordered(
                    create_differential_operator(diff_ops, seq),
                    current_derivative_subnode_rest, NodesValue{})...);
            },
            multinomial_sequences);

        constexpr auto flags_next_derivatives = std::apply(
            [dnin](auto... next_derivative) {
                return std::integer_sequence<
                    bool, monomial_included(next_derivative, dnin)...>{};
            },
            next_derivatives);

        constexpr auto next_derivatives_filtered =
            filter(next_derivatives, flags_next_derivatives);

        constexpr auto multinomial_sequences_filtered =
            filter(multinomial_sequences, flags_next_derivatives);

        static_assert(size(next_derivatives_filtered));

        constexpr auto locations =
            locate_new_vals_update_buffer_types_size<Last>(
                next_derivatives_filtered, dn, dnn, dnin);

        constexpr auto flags_only_new = std::apply(
            [](auto... location) {
                return std::integer_sequence<
                    bool,
                    std::is_same_v<on_buffer_t, decltype(location.first)>...>{};
            },
            locations);

        constexpr auto next_derivatives_new =
            filter(next_derivatives_filtered, flags_only_new);

        constexpr auto locations_new = filter(locations, flags_only_new);
        constexpr auto next_derivatives_new_with_pos =
            add_position(next_derivatives_new, locations_new);
        constexpr auto dnn_new =
            merge_sorted(next_derivatives_new_with_pos, dnn, NodesValue{});

        constexpr std::size_t size_cur_new =
            SizeCur + size(next_derivatives_new_with_pos) -
            static_cast<std::size_t>(free_on_buffer_size(current_node_loc));
        constexpr std::size_t size_max_new = std::max(SizeMax, size_cur_new);

        return treat_nodes_add_buffer_size<Last, CalcTree, size_max_new,
                                           size_cur_new, N + 1>(pn, dn, dnn_new,
                                                                dnin);
    }
}

template <std::size_t Last, class CalcTree, std::size_t SizeMax,
          std::size_t SizeCur, class PrimalSubNode1, class PrimalSubNode2,
          class DerivativeNodes, class DerivativeNodeInputs>
constexpr auto treat_nodes_specialized_buffer_size(
    add_t<PrimalSubNode1, PrimalSubNode2> /* pn */, DerivativeNodes dn,
    DerivativeNodeInputs dnin) {

    using NodesValue = CalcTree::ValuesTupleInverse;
    constexpr auto ordered_pair = detail::sort_pair(
        std::tuple<PrimalSubNode1, PrimalSubNode2>{}, NodesValue{});

    return treat_nodes_add_buffer_size<Last, CalcTree, SizeMax, SizeCur>(
        ordered_pair, dn, std::tuple<>{}, dnin);
}

template <std::size_t Last, std::size_t N, class CalcTree, std::size_t SizeMax,
          std::size_t SizeCur, std::size_t MaxOrder, std::size_t PrevOrder = 0,
          template <class> class Univariate, class PrimalSubNode,
          class DerivativeNodes, class DerivativeNodeNew,
          class DerivativeNodeInputs>
constexpr auto treat_nodes_univariate_buffer_size(Univariate<PrimalSubNode> pn,
                                                  DerivativeNodes dn,
                                                  DerivativeNodeNew dnn,
                                                  DerivativeNodeInputs dnin) {

    if constexpr (N == 0) {
        return std::make_tuple(dnn,
                               std::integral_constant<std::size_t, SizeMax>{},
                               std::integral_constant<std::size_t, SizeCur>{});
    } else {
        constexpr std::size_t currentN = N - 1;
        constexpr auto current_node_der = std::get<0>(std::get<currentN>(dn));
        constexpr auto current_node_loc = std::get<1>(std::get<currentN>(dn));

#if LOG
        std::cout << "treating derivative tree node" << std::endl;
        std::cout << type_name2<decltype(current_node_der)>() << std::endl;
#endif

        constexpr auto current_derivative_subnode = head(current_node_der);
        constexpr auto current_derivative_subnode_rest = tail(current_node_der);

        constexpr auto pow = get_power(current_derivative_subnode);
        constexpr auto rest_power = power(current_derivative_subnode_rest);
        constexpr auto expansion_types =
            expand_univariate<PrimalSubNode, MaxOrder - rest_power, pow>();

        using NodesValue = CalcTree::ValuesTupleInverse;

        constexpr auto next_derivatives = std::apply(
            [current_derivative_subnode_rest](const auto... singletype) {
                return std::make_tuple(multiply_differential_operator(
                    singletype, current_derivative_subnode_rest,
                    NodesValue{})...);
            },
            expansion_types);

        constexpr auto flags_next_derivatives = std::apply(
            [dnin](auto... next_derivative) {
                return std::integer_sequence<
                    bool, monomial_included(next_derivative, dnin)...>{};
            },
            next_derivatives);

        constexpr auto next_derivatives_filtered =
            filter(next_derivatives, flags_next_derivatives);

        static_assert(size(next_derivatives_filtered));

        constexpr auto locations =
            locate_new_vals_update_buffer_types_size<Last>(
                next_derivatives_filtered, dn, dnn, dnin);

        constexpr auto flags_only_new = std::apply(
            [](auto... location) {
                return std::integer_sequence<
                    bool,
                    std::is_same_v<on_buffer_t, decltype(location.first)>...>{};
            },
            locations);

        constexpr auto next_derivatives_new =
            filter(next_derivatives_filtered, flags_only_new);

        constexpr auto locations_new = filter(locations, flags_only_new);
        constexpr auto next_derivatives_new_with_pos =
            add_position(next_derivatives_new, locations_new);
        constexpr auto dnn_new =
            merge_sorted(next_derivatives_new_with_pos, dnn, NodesValue{});

        constexpr std::size_t size_cur_new =
            SizeCur + size(next_derivatives_new_with_pos) -
            static_cast<std::size_t>(free_on_buffer_size(current_node_loc));
        constexpr std::size_t size_max_new = std::max(SizeMax, size_cur_new);

        return treat_nodes_univariate_buffer_size<
            Last, currentN, CalcTree, size_max_new, size_cur_new, MaxOrder>(
            pn, dn, dnn_new, dnin);
    }
}

template <std::size_t Last, class CalcTree, std::size_t SizeMax,
          std::size_t SizeCur, template <class> class Univariate,
          class PrimalSubNode, class DerivativeNodes,
          class DerivativeNodeInputs>
constexpr auto treat_nodes_specialized_buffer_size(Univariate<PrimalSubNode> pn,
                                                   DerivativeNodes dn,
                                                   DerivativeNodeInputs dnin) {
    using PrimalNode = Univariate<PrimalSubNode>;
    constexpr auto MaxOrder = detail::max_orders(dnin);

    // we go over univariate derivatives in inverse lexicographic order.
    // why? because it makes sense to calculate the coefficients increasing in
    // power. the lexicographic order is decreasing in powers.
    return treat_nodes_univariate_buffer_size<Last, Last, CalcTree, SizeMax,
                                              SizeCur, MaxOrder>(
        pn, dn, std::tuple<>{}, dnin);
}

template <class CalcTree, std::size_t SizeMax, std::size_t SizeCur,
          class PrimalNode, class DerivativeNodes, class DerivativeNodeInputs>
constexpr auto treat_node_buffer_size(PrimalNode nd, DerivativeNodes dn,
                                      DerivativeNodeInputs dnin) {

    constexpr auto Last = find_first_type_not<DerivativeNodes, PrimalNode>();

    constexpr auto res =
        treat_nodes_specialized_buffer_size<Last, CalcTree, SizeMax, SizeCur>(
            nd, dn, dnin);

    constexpr std::tuple_element_t<0, decltype(res)> new_dn;
    constexpr std::tuple_element_t<1, decltype(res)> new_size_max;
    constexpr std::tuple_element_t<2, decltype(res)> new_size_cur;

    using NodesValue = CalcTree::ValuesTupleInverse;

    constexpr auto dn_remaining = sub_tuple<Last, size(dn) - 1>(dn);
    constexpr auto dn_new_and_remaining =
        merge_sorted(new_dn, dn_remaining, NodesValue{});

    return std::make_tuple(
        dn_new_and_remaining,
        std::integral_constant<std::size_t, new_size_max>{},
        std::integral_constant<std::size_t, new_size_cur()>{});
}

template <class CalcTree, std::size_t SizeMax = 0, std::size_t SizeCur = 0,
          std::size_t N = 0, class DerivativeNode, class DerivativeNodeInputs>
constexpr auto backpropagate_buffer_size(DerivativeNode dn,
                                         DerivativeNodeInputs dnin)
    -> std::size_t {
    using PrimalNodes = CalcTree::ValuesTupleInverse;
    constexpr auto current_primal_node = std::get<N>(PrimalNodes{});

    if constexpr (!is_input(current_primal_node)) {
#if LOG
        std::cout << "treating calc tree node" << std::endl;
        std::cout << type_name2<decltype(current_primal_node)>() << std::endl;
#endif

        constexpr auto res = treat_node_buffer_size<CalcTree, SizeMax, SizeCur>(
            current_primal_node, dn, dnin);

        constexpr std::tuple_element_t<0, decltype(res)> new_dn;
        constexpr std::tuple_element_t<1, decltype(res)> new_size_max;
        constexpr std::tuple_element_t<2, decltype(res)> new_size_cur;

        return backpropagate_buffer_size<CalcTree, new_size_max(),
                                         new_size_cur(), N + 1>(new_dn, dnin);
    } else {
        return SizeMax;
    }
}

} // namespace adhoc4::detail

#endif // ADHOC4_BACKPROPAGATOR_BUFFER_SIZE_HPP
