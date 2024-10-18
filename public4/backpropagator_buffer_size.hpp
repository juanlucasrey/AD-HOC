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

template <class CalcTree, std::size_t Result, class BufferTypes,
          std::size_t MaxOrder, std::size_t N = 0, class PrimalSubNodeOrdered1,
          class PrimalSubNodeOrdered2, class DerivativeNodeLocations,
          class DerivativeNodes, class InterfaceTypes, class DerivativeNodeNew,
          class DerivativeNodeInputs>
constexpr auto treat_nodes_mul_buffer_size(
    std::tuple<PrimalSubNodeOrdered1, PrimalSubNodeOrdered2> pn,
    DerivativeNodeLocations dnl, DerivativeNodes dn, InterfaceTypes it,
    DerivativeNodeNew dnn, DerivativeNodeInputs dnin) {

    if constexpr (N == std::tuple_size_v<DerivativeNodeLocations>) {
        return std::make_tuple(BufferTypes{}, dnn,
                               std::integral_constant<std::size_t, Result>{});
    } else {
        constexpr auto current_node_der = std::get<N>(dn);

#if LOG
        std::cout << "treating derivative tree node" << std::endl;
        std::cout << type_name2<decltype(current_node_der)>() << std::endl;
#endif

        constexpr auto current_node_der_loc = std::get<N>(dnl);
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

#if LOG
        std::cout << "next_derivatives_filtered" << std::endl;
        std::cout << type_name2<decltype(next_derivatives_filtered)>()
                  << std::endl;
#endif

        constexpr auto multinomial_sequences_filtered =
            filter(multinomial_sequences, flags_next_derivatives);

        static_assert(size(next_derivatives_filtered));

        constexpr auto bt_free = free_on_buffer_size(
            current_node_der_loc, current_node_der, BufferTypes{});
        constexpr auto ResultFree = std::max(Result, size(bt_free));

#if LOG
        std::cout << "treat_nodes_univariate_buffer_size" << std::endl;
        std::cout << ResultFree << std::endl;
        std::cout << "bt_free" << std::endl;
        std::cout << type_name2<decltype(bt_free)>() << std::endl;
#endif

        constexpr auto locations =
            locate_new_vals(next_derivatives_filtered, it, bt_free);

#if LOG
        std::cout << "locations" << std::endl;
        std::cout << type_name2<decltype(locations)>() << std::endl;
#endif

        constexpr auto bt_new = update_buffer_types_size(
            next_derivatives_filtered, locations, bt_free);

#if LOG
        std::cout << "bt_new" << std::endl;
        std::cout << type_name2<decltype(bt_new)>() << std::endl;
#endif

        constexpr auto ResultAdd = std::max(ResultFree, size(bt_new));

#if LOG
        std::cout << "treat_nodes_univariate_buffer_size" << std::endl;
        std::cout << ResultAdd << std::endl;
#endif

        constexpr auto flags_only_new = std::apply(
            [](auto... location) {
                return std::integer_sequence<
                    bool,
                    std::is_same_v<on_buffer_new_t, decltype(location)>...>{};
            },
            locations);

#if LOG
        std::cout << "flags_only_new" << std::endl;
        std::cout << type_name2<decltype(flags_only_new)>() << std::endl;
#endif

        constexpr auto next_derivatives_new =
            filter(next_derivatives_filtered, flags_only_new);

#if LOG
        std::cout << "next_derivatives_new" << std::endl;
        std::cout << type_name2<decltype(next_derivatives_new)>() << std::endl;
#endif

        constexpr auto dnn_new =
            merge_sorted2(next_derivatives_new, dnn, NodesValue{});

#if LOG
        std::cout << "dnn_new" << std::endl;
        std::cout << type_name2<decltype(dnn_new)>() << std::endl;
#endif

        return treat_nodes_mul_buffer_size<CalcTree, ResultAdd,
                                           decltype(bt_new), MaxOrder, N + 1>(
            pn, dnl, dn, it, dnn_new, dnin);
    }
}

template <class CalcTree, std::size_t Result, class BufferTypes,
          class PrimalSubNode1, class PrimalSubNode2,
          class DerivativeNodeLocations, class DerivativeNodes,
          class InterfaceTypes, class DerivativeNodeInputs>
constexpr auto treat_nodes_specialized_buffer_size(
    mul_t<PrimalSubNode1, PrimalSubNode2> /* pn */, DerivativeNodeLocations dnl,
    DerivativeNodes dn, InterfaceTypes it, DerivativeNodeInputs dnin) {

#if LOG
    std::cout << "treat_nodes_specialized_buffer_size" << std::endl;
    std::cout << Result << std::endl;
#endif

    using NodesValue = CalcTree::ValuesTupleInverse;
    constexpr auto ordered_pair = detail::sort_pair(
        std::tuple<PrimalSubNode1, PrimalSubNode2>{}, NodesValue{});

    constexpr auto MaxOrder = detail::max_orders(dnin);

    return treat_nodes_mul_buffer_size<CalcTree, Result, BufferTypes, MaxOrder>(
        ordered_pair, dnl, dn, it, std::tuple<>{}, dnin);
}

template <class CalcTree, std::size_t Result, class BufferTypes,
          std::size_t N = 0, class PrimalSubNodeOrdered1,
          class PrimalSubNodeOrdered2, class DerivativeNodeLocations,
          class DerivativeNodes, class InterfaceTypes, class DerivativeNodeNew,
          class DerivativeNodeInputs>
constexpr auto treat_nodes_add_buffer_size(
    std::tuple<PrimalSubNodeOrdered1, PrimalSubNodeOrdered2> pn,
    DerivativeNodeLocations dnl, DerivativeNodes dn, InterfaceTypes it,
    DerivativeNodeNew dnn, DerivativeNodeInputs dnin) {

    if constexpr (N == std::tuple_size_v<DerivativeNodeLocations>) {
        return std::make_tuple(BufferTypes{}, dnn,
                               std::integral_constant<std::size_t, Result>{});

    } else {
        constexpr auto current_node_der = std::get<N>(dn);

#if LOG
        std::cout << "treating derivative tree node" << std::endl;
        std::cout << type_name2<decltype(current_node_der)>() << std::endl;
#endif

        constexpr auto current_node_der_loc = std::get<N>(dnl);
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

#if LOG
        std::cout << "next_derivatives" << std::endl;
        std::cout << type_name2<decltype(next_derivatives)>() << std::endl;
#endif

        constexpr auto flags_next_derivatives = std::apply(
            [dnin](auto... next_derivative) {
                return std::integer_sequence<
                    bool, monomial_included(next_derivative, dnin)...>{};
            },
            next_derivatives);

#if LOG
        std::cout << "flags_next_derivatives" << std::endl;
        std::cout << type_name2<decltype(flags_next_derivatives)>()
                  << std::endl;
#endif

        constexpr auto next_derivatives_filtered =
            filter(next_derivatives, flags_next_derivatives);

#if LOG
        std::cout << "next_derivatives_filtered" << std::endl;
        std::cout << type_name2<decltype(next_derivatives_filtered)>()
                  << std::endl;
#endif

        constexpr auto multinomial_sequences_filtered =
            filter(multinomial_sequences, flags_next_derivatives);

        static_assert(size(next_derivatives_filtered));

        constexpr auto bt_free = free_on_buffer_size(
            current_node_der_loc, current_node_der, BufferTypes{});
        constexpr auto ResultFree = std::max(Result, size(bt_free));

#if LOG
        std::cout << "treat_nodes_univariate_buffer_size" << std::endl;
        std::cout << ResultFree << std::endl;
        std::cout << "bt_free" << std::endl;
        std::cout << type_name2<decltype(bt_free)>() << std::endl;
#endif

        constexpr auto locations =
            locate_new_vals(next_derivatives_filtered, it, bt_free);

#if LOG
        std::cout << "locations" << std::endl;
        std::cout << type_name2<decltype(locations)>() << std::endl;
#endif

        constexpr auto bt_new = update_buffer_types_size(
            next_derivatives_filtered, locations, bt_free);

#if LOG
        std::cout << "bt_new" << std::endl;
        std::cout << type_name2<decltype(bt_new)>() << std::endl;
#endif

        constexpr auto ResultAdd = std::max(ResultFree, size(bt_new));

#if LOG
        std::cout << "treat_nodes_univariate_buffer_size" << std::endl;
        std::cout << ResultAdd << std::endl;
#endif

        constexpr auto flags_only_new = std::apply(
            [](auto... location) {
                return std::integer_sequence<
                    bool,
                    std::is_same_v<on_buffer_new_t, decltype(location)>...>{};
            },
            locations);

        constexpr auto next_derivatives_new =
            filter(next_derivatives_filtered, flags_only_new);

        constexpr auto dnn_new =
            merge_sorted2(next_derivatives_new, dnn, NodesValue{});

        return treat_nodes_add_buffer_size<CalcTree, ResultAdd,
                                           decltype(bt_new), N + 1>(
            pn, dnl, dn, it, dnn_new, dnin);
    }
}

template <class CalcTree, std::size_t Result, class BufferTypes,
          class PrimalSubNode1, class PrimalSubNode2,
          class DerivativeNodeLocations, class DerivativeNodes,
          class InterfaceTypes, class DerivativeNodeInputs>
constexpr auto treat_nodes_specialized_buffer_size(
    add_t<PrimalSubNode1, PrimalSubNode2> /* pn */, DerivativeNodeLocations dnl,
    DerivativeNodes dn, InterfaceTypes it, DerivativeNodeInputs dnin) {

#if LOG
    std::cout << "treat_nodes_specialized_buffer_size" << std::endl;
    std::cout << Result << std::endl;
#endif

    using NodesValue = CalcTree::ValuesTupleInverse;
    constexpr auto ordered_pair = detail::sort_pair(
        std::tuple<PrimalSubNode1, PrimalSubNode2>{}, NodesValue{});

    return treat_nodes_add_buffer_size<CalcTree, Result, BufferTypes>(
        ordered_pair, dnl, dn, it, std::tuple<>{}, dnin);
}

template <std::size_t N, class CalcTree, std::size_t Result, class BufferTypes,
          std::size_t MaxOrder, std::size_t PrevOrder = 0,
          template <class> class Univariate, class PrimalSubNode,
          class DerivativeNodeLocations, class DerivativeNodes,
          class InterfaceTypes, class DerivativeNodeNew,
          class DerivativeNodeInputs>
constexpr auto treat_nodes_univariate_buffer_size(Univariate<PrimalSubNode> pn,
                                                  DerivativeNodeLocations dnl,
                                                  DerivativeNodes dn,
                                                  InterfaceTypes it,
                                                  DerivativeNodeNew dnn,
                                                  DerivativeNodeInputs dnin) {

#if LOG
    std::cout << "treat_nodes_univariate_buffer_size" << std::endl;
    std::cout << Result << std::endl;
#endif

    if constexpr (N == 0) {
        return std::make_tuple(BufferTypes{}, dnn,
                               std::integral_constant<std::size_t, Result>{});
    } else {
        constexpr std::size_t currentN = N - 1;
        constexpr auto current_node_der = std::get<currentN>(dn);

#if LOG
        std::cout << "treating derivative tree node" << std::endl;
        std::cout << type_name2<decltype(current_node_der)>() << std::endl;
#endif

        constexpr auto current_node_der_loc = std::get<currentN>(dnl);
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

#if LOG
        std::cout << "next_derivatives" << std::endl;
        std::cout << type_name2<decltype(next_derivatives)>() << std::endl;
#endif

        constexpr auto flags_next_derivatives = std::apply(
            [dnin](auto... next_derivative) {
                return std::integer_sequence<
                    bool, monomial_included(next_derivative, dnin)...>{};
            },
            next_derivatives);

#if LOG
        std::cout << "flags_next_derivatives" << std::endl;
        std::cout << type_name2<decltype(flags_next_derivatives)>()
                  << std::endl;
#endif

        constexpr auto next_derivatives_filtered =
            filter(next_derivatives, flags_next_derivatives);

#if LOG
        std::cout << "next_derivatives_filtered" << std::endl;
        std::cout << type_name2<decltype(next_derivatives_filtered)>()
                  << std::endl;
#endif

        static_assert(size(next_derivatives_filtered));

        constexpr auto bt_free = free_on_buffer_size(
            current_node_der_loc, current_node_der, BufferTypes{});
        constexpr auto ResultFree = std::max(Result, size(bt_free));

#if LOG
        std::cout << "treat_nodes_univariate_buffer_size" << std::endl;
        std::cout << ResultFree << std::endl;
        std::cout << "bt_free" << std::endl;
        std::cout << type_name2<decltype(bt_free)>() << std::endl;
#endif

        constexpr auto locations =
            locate_new_vals(next_derivatives_filtered, it, bt_free);

#if LOG
        std::cout << "locations" << std::endl;
        std::cout << type_name2<decltype(locations)>() << std::endl;
#endif

        constexpr auto bt_new = update_buffer_types_size(
            next_derivatives_filtered, locations, bt_free);

#if LOG
        std::cout << "bt_new" << std::endl;
        std::cout << type_name2<decltype(bt_new)>() << std::endl;
#endif

        constexpr auto ResultAdd = std::max(ResultFree, size(bt_new));

#if LOG
        std::cout << "treat_nodes_univariate_buffer_size" << std::endl;
        std::cout << ResultAdd << std::endl;
#endif

        constexpr auto flags_only_new = std::apply(
            [](auto... location) {
                return std::integer_sequence<
                    bool,
                    std::is_same_v<on_buffer_new_t, decltype(location)>...>{};
            },
            locations);

#if LOG
        std::cout << "flags_only_new" << std::endl;
        std::cout << type_name2<decltype(flags_only_new)>() << std::endl;
#endif

        constexpr auto next_derivatives_new =
            filter(next_derivatives_filtered, flags_only_new);

#if LOG
        std::cout << "next_derivatives_new" << std::endl;
        std::cout << type_name2<decltype(next_derivatives_new)>() << std::endl;
#endif

        constexpr auto dnn_new =
            merge_sorted2(next_derivatives_new, dnn, NodesValue{});

#if LOG
        std::cout << "dnn_new" << std::endl;
        std::cout << type_name2<decltype(dnn_new)>() << std::endl;
#endif

        return treat_nodes_univariate_buffer_size<currentN, CalcTree, ResultAdd,
                                                  decltype(bt_new), MaxOrder>(
            pn, dnl, dn, it, dnn_new, dnin);
    }
}

template <class CalcTree, std::size_t Result, class BufferTypes,
          template <class> class Univariate, class PrimalSubNode,
          class DerivativeNodeLocations, class DerivativeNodes,
          class InterfaceTypes, class DerivativeNodeInputs>
constexpr auto treat_nodes_specialized_buffer_size(Univariate<PrimalSubNode> pn,
                                                   DerivativeNodeLocations dnl,
                                                   DerivativeNodes dn,
                                                   InterfaceTypes it,
                                                   DerivativeNodeInputs dnin) {

#if LOG
    std::cout << "treat_nodes_specialized_buffer_size" << std::endl;
    std::cout << Result << std::endl;
#endif

    using PrimalNode = Univariate<PrimalSubNode>;
    constexpr auto MaxOrder = detail::max_orders(dnin);

    // we go over univariate derivatives in inverse lexicographic order.
    // why? because it makes sense to calculate the coefficients increasing in
    // power. the lexicographic order is decreasing in powers.
    constexpr std::size_t End = std::tuple_size_v<DerivativeNodeLocations>;
    return treat_nodes_univariate_buffer_size<End, CalcTree, Result,
                                              BufferTypes, MaxOrder>(
        pn, dnl, dn, it, std::tuple<>{}, dnin);
}

template <class CalcTree, std::size_t Result, class BufferTypes,
          class PrimalNode, class DerivativeNodeLocation, class DerivativeNodes,
          class InterfaceTypes, class DerivativeNodeInputs>
constexpr auto treat_node_buffer_size(PrimalNode nd, DerivativeNodeLocation dnl,
                                      DerivativeNodes dn, InterfaceTypes it,
                                      DerivativeNodeInputs dnin) {

#if LOG
    std::cout << "treat_node_buffer_size" << std::endl;
    std::cout << Result << std::endl;
    std::cout << "dnl" << std::endl;
    std::cout << type_name2<decltype(dnl)>() << std::endl;
    std::cout << "dn" << std::endl;
    std::cout << type_name2<decltype(dn)>() << std::endl;
#endif

    constexpr auto flags_derivative_nodes = std::apply(
        [nd](auto... type) {
            return std::integer_sequence<bool,
                                         detail::first_type_is(type, nd)...>{};
        },
        dn);

#if LOG
    std::cout << "flags_derivative_nodes" << std::endl;
    std::cout << type_name2<decltype(flags_derivative_nodes)>() << std::endl;
#endif

    constexpr auto separated_derivative_nodes_loc =
        separate(dnl, flags_derivative_nodes);

#if LOG
    std::cout << "separated_derivative_nodes_loc" << std::endl;
    std::cout << type_name2<decltype(separated_derivative_nodes_loc)>()
              << std::endl;
#endif

    constexpr auto separated_derivative_nodes =
        separate(dn, flags_derivative_nodes);

#if LOG
    std::cout << "separated_derivative_nodes" << std::endl;
    std::cout << type_name2<decltype(separated_derivative_nodes)>()
              << std::endl;
#endif

    constexpr auto return_pair =
        treat_nodes_specialized_buffer_size<CalcTree, Result, BufferTypes>(
            nd, std::get<0>(separated_derivative_nodes_loc),
            std::get<0>(separated_derivative_nodes), it, dnin);

    constexpr auto bt_new = std::get<0>(return_pair);
    constexpr auto dn_new = std::get<1>(return_pair);

#if LOG
    std::cout << "dn_new" << std::endl;
    std::cout << type_name2<decltype(dn_new)>() << std::endl;
#endif

    constexpr auto new_size = std::get<2>(return_pair);

#if LOG
    std::cout << "treat_node_buffer_size" << std::endl;
    std::cout << new_size() << std::endl;
#endif

    constexpr auto ResultNew = std::max(Result, new_size());

#if LOG
    std::cout << "treat_node_buffer_size" << std::endl;
    std::cout << ResultNew << std::endl;
#endif

    constexpr auto dn_remaining = std::get<1>(separated_derivative_nodes);

#if LOG
    std::cout << "dn_remaining" << std::endl;
    std::cout << type_name2<decltype(dn_remaining)>() << std::endl;
#endif

    using NodesValue = CalcTree::ValuesTupleInverse;

    constexpr auto dn_new_and_remaining =
        merge_sorted2(dn_new, dn_remaining, NodesValue{});

#if LOG
    std::cout << "dn_new_and_remaining" << std::endl;
    std::cout << type_name2<decltype(dn_new_and_remaining)>() << std::endl;
#endif

    // todo better for multiple output nodes
    constexpr auto dnl_new =
        make_tuple_same<size(dn_new_and_remaining)>(on_buffer_t{});

    return std::make_tuple(bt_new, dnl_new, dn_new_and_remaining,
                           std::integral_constant<std::size_t, ResultNew>{});
}

template <class CalcTree, std::size_t Result = 0, std::size_t N = 0,
          class BufferTypes = std::tuple<>, class DerivativeNodeLocation,
          class DerivativeNodes, class InterfaceTypes,
          class DerivativeNodeInputs>
constexpr auto backpropagate_buffer_size(DerivativeNodeLocation const dnl,
                                         DerivativeNodes dn, InterfaceTypes it,
                                         DerivativeNodeInputs dnin)
    -> std::size_t {

#if LOG
    std::cout << "backpropagate_buffer_size" << std::endl;
    std::cout << Result << std::endl;
#endif

    using PrimalNodes = CalcTree::ValuesTupleInverse;
    constexpr auto current_primal_node = std::get<N>(PrimalNodes{});

    if constexpr (!is_input(current_primal_node)) {
#if LOG
        std::cout << "treating calc tree node" << std::endl;
        std::cout << type_name2<decltype(current_primal_node)>() << std::endl;
#endif

        constexpr auto res =
            treat_node_buffer_size<CalcTree, Result, BufferTypes>(
                current_primal_node, dnl, dn, it, dnin);

        constexpr auto bt_new = std::get<0>(res);
        constexpr auto dnl_new = std::get<1>(res);
        constexpr auto dn_new = std::get<2>(res);
        constexpr auto new_size = std::get<3>(res)();
#if LOG
        std::cout << "dnl_new" << std::endl;
        std::cout << type_name2<decltype(dnl_new)>() << std::endl;

        std::cout << "dn_new" << std::endl;
        std::cout << type_name2<decltype(dn_new)>() << std::endl;

        std::cout << "backpropagate_buffer_size" << std::endl;
        std::cout << new_size << std::endl;
#endif

        constexpr auto ResultNew = std::max(Result, new_size);

#if LOG
        std::cout << "backpropagate_buffer_size" << std::endl;
        std::cout << ResultNew << std::endl;
#endif

        return backpropagate_buffer_size<CalcTree, ResultNew, N + 1,
                                         decltype(bt_new)>(dnl_new, dn_new, it,
                                                           dnin);
    } else {
        return Result;
    }
}

} // namespace adhoc4::detail

#endif // ADHOC4_BACKPROPAGATOR_BUFFER_SIZE_HPP
