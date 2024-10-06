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

#ifndef ADHOC4_BACKPROPAGATOR_AUX_HPP
#define ADHOC4_BACKPROPAGATOR_AUX_HPP

#include "backpropagator_tools.hpp"
#include "combinatorics/combinations.hpp"
#include "combinatorics/multinomial_coefficient_index_sequence.hpp"
#include "differential_operator.hpp"
#include "monomial_included.hpp"
#include "sort.hpp"
#include "utils/bivariate.hpp"
#include "utils/index_sequence.hpp"
#include "utils/univariate.hpp"

#include <array>
#include <tuple>

#include "../tests4/type_name.hpp"
#include <iostream>

namespace adhoc4::detail {

template <std::size_t Power, std::size_t N = 0, std::size_t OutSize,
          std::size_t MaxOrder, class ExpansionTypes>
auto multiply_mult(std::array<double, OutSize> &arrayout,
                   std::array<double, MaxOrder> const &powers_val1,
                   std::array<double, MaxOrder> const &powers_val2,
                   ExpansionTypes et) {
    static_assert(OutSize == std::tuple_size_v<ExpansionTypes>);
    if constexpr (N < OutSize) {
        constexpr auto current_type = std::get<N>(et);

        //     if constexpr (MultinomialCoeff2(current_sequence) != 1) {
        //         arrayout[N] *=
        //             static_cast<double>(MultinomialCoeff2(current_sequence));
        //     }
        multiply_mult<Power, N + 1>(arrayout, powers_val1, powers_val2, et);
    }
}

template <std::size_t MaxOrder, std::size_t PrevOrder, class PrimalSubNode1,
          class PrimalSubNode2, class DerivativeNodeLocations,
          class DerivativeNodes, class CalcTree, class InterfaceTypes,
          class InterfaceArray, class BufferTypes, class BufferArray,
          class DerivativeNodeNew, class DerivativeNodeInputs>
auto treat_derivative_nodes_mul(mul_t<PrimalSubNode1, PrimalSubNode2> pn,
                                DerivativeNodeLocations dnl, DerivativeNodes dn,
                                CalcTree const &ct, InterfaceTypes it,
                                InterfaceArray &ia, BufferTypes bt,
                                BufferArray &ba,
                                std::array<double, MaxOrder> const &powers_val1,
                                std::array<double, MaxOrder> const &powers_val2,
                                DerivativeNodeNew dnn,
                                DerivativeNodeInputs dnin) {

    std::cout << "dnl" << std::endl;
    std::cout << type_name2<decltype(dnl)>() << std::endl;

    std::cout << "dn" << std::endl;
    std::cout << type_name2<decltype(dn)>() << std::endl;

    if constexpr (size(dnl)) {
        constexpr auto current_node_der_loc = head(dnl);
        constexpr auto current_node_der = head(dn);
        constexpr auto current_derivative_subnode = head(current_node_der);
        constexpr auto current_derivative_subnode_rest = tail(current_node_der);
        constexpr auto pow = get_power(current_derivative_subnode);

        // std::cout << "pow" << std::endl;
        // std::cout << pow << std::endl;

        // std::cout << "current_derivative_subnode" << std::endl;
        // std::cout << type_name2<decltype(current_derivative_subnode)>()
        //           << std::endl;

        using NodesValue = CalcTree::ValuesTupleInverse;
        constexpr auto expansion_types =
            expand_bivariate(current_derivative_subnode, NodesValue{});

        // std::cout << "expansion_types" << std::endl;
        // std::cout << type_name2<decltype(expansion_types)>() << std::endl;

        constexpr auto next_derivatives = std::apply(
            [current_derivative_subnode_rest](const auto... singletype) {
                return std::make_tuple(multiply_ordered_tuple(
                    singletype, current_derivative_subnode_rest,
                    NodesValue{})...);
            },
            expansion_types);

        // std::cout << "next_derivatives" << std::endl;
        // std::cout << type_name2<decltype(next_derivatives)>() << std::endl;

        // std::cout << "dnin" << std::endl;
        // std::cout << type_name2<decltype(dnin)>() << std::endl;

        constexpr auto flags_next_derivatives = std::apply(
            [dnin](auto... next_derivative) {
                return std::tuple_cat(
                    std::conditional_t<monomial_included(next_derivative, dnin),
                                       std::tuple<std::true_type>,
                                       std::tuple<std::false_type>>{}...);
            },
            next_derivatives);

        // std::cout << "flags_next_derivatives" << std::endl;
        // std::cout << type_name2<decltype(flags_next_derivatives)>()
        //           << std::endl;

        constexpr auto next_derivatives_filtered =
            filter(next_derivatives, flags_next_derivatives);

        // std::cout << "next_derivatives_filtered" << std::endl;
        // std::cout << type_name2<decltype(next_derivatives_filtered)>()
        //           << std::endl;

        static_assert(size(next_derivatives_filtered));

        double const this_val_derivative = get_differential_operator_value(
            current_node_der_loc, current_node_der, it, ia, bt, ba);

        constexpr auto bt_free =
            free_on_buffer(current_node_der_loc, current_node_der, bt);

        constexpr auto next_derivatives_size =
            std::tuple_size_v<decltype(next_derivatives_filtered)>;

        std::array<double, next_derivatives_size> next_derivatives_values;
        next_derivatives_values.fill(this_val_derivative);

        constexpr auto expansion_types_filtered =
            filter(expansion_types, flags_next_derivatives);

        multiply_mult<pow>(next_derivatives_values, powers_val1, powers_val2,
                           expansion_types_filtered);

        // constexpr std::size_t pow = 0;
        constexpr auto bt_new = bt;
        constexpr auto dnn_new = dnn;

        return treat_derivative_nodes_mul<MaxOrder, pow>(
            pn, tail(dnl), tail(dn), ct, it, ia, bt_new, ba, powers_val1,
            powers_val2, dnn_new, dnin);
    } else {
        return std::make_tuple(bt, dnn);
    }

    // return std::make_tuple(bt, dnn);
}

template <class PrimalSubNode1, class PrimalSubNode2,
          class DerivativeNodeLocations, class DerivativeNodes, class CalcTree,
          class InterfaceTypes, class InterfaceArray, class BufferTypes,
          class BufferArray, std::size_t MaxOrder, class DerivativeNodeInputs>
auto treat_derivative_nodes(mul_t<PrimalSubNode1, PrimalSubNode2> pn,
                            DerivativeNodeLocations dnl, DerivativeNodes dn,
                            CalcTree const &ct, InterfaceTypes it,
                            InterfaceArray &ia, BufferTypes bt, BufferArray &ba,
                            std::array<double, MaxOrder> & /* ua */,
                            DerivativeNodeInputs dnin) {
    using NodesValue = CalcTree::ValuesTupleInverse;

    // are nodes inverted??
    constexpr auto id1_less_than_id2 =
        static_cast<bool>(get_first_type_idx(NodesValue{}, PrimalSubNode1{}) >=
                          get_first_type_idx(NodesValue{}, PrimalSubNode2{}));
    double const val1 =
        id1_less_than_id2 ? ct.get(PrimalSubNode2{}) : ct.get(PrimalSubNode1{});
    double const val2 =
        id1_less_than_id2 ? ct.get(PrimalSubNode1{}) : ct.get(PrimalSubNode2{});

    auto const powers_val1 = detail::powers<MaxOrder>(val1);
    auto const powers_val2 = detail::powers<MaxOrder>(val2);

    return treat_derivative_nodes_mul<MaxOrder, 0>(pn, dnl, dn, ct, it, ia, bt,
                                                   ba, powers_val1, powers_val2,
                                                   std::tuple<>{}, dnin);
}

template <std::size_t MaxOrder, std::size_t PrevOrder, class PrimalSubNode1,
          class PrimalSubNode2, class DerivativeNodeLocations,
          class DerivativeNodes, class CalcTree, class InterfaceTypes,
          class InterfaceArray, class BufferTypes, class BufferArray,
          class DerivativeNodeNew, class DerivativeNodeInputs>
auto treat_derivative_nodes_sum(add_t<PrimalSubNode1, PrimalSubNode2> pn,
                                DerivativeNodeLocations dnl, DerivativeNodes dn,
                                CalcTree const &ct, InterfaceTypes it,
                                InterfaceArray &ia, BufferTypes bt,
                                BufferArray &ba, DerivativeNodeNew dnn,
                                DerivativeNodeInputs dnin) {

    if constexpr (size(dnl)) {
        constexpr auto current_node_der_loc = head(dnl);
        constexpr auto current_node_der = head(dn);
        constexpr auto current_derivative_subnode = head(current_node_der);
        constexpr auto current_derivative_subnode_rest = tail(current_node_der);
        constexpr auto pow = get_power(current_derivative_subnode);

        using NodesValue = CalcTree::ValuesTupleInverse;
        constexpr auto expansion_types =
            expand_bivariate(current_derivative_subnode, NodesValue{});

        constexpr auto next_derivatives = std::apply(
            [current_derivative_subnode_rest](const auto... singletype) {
                return std::make_tuple(multiply_ordered_tuple(
                    singletype, current_derivative_subnode_rest,
                    NodesValue{})...);
            },
            expansion_types);

        constexpr auto flags_next_derivatives = std::apply(
            [dnin](auto... next_derivative) {
                return std::tuple_cat(
                    std::conditional_t<monomial_included(next_derivative, dnin),
                                       std::tuple<std::true_type>,
                                       std::tuple<std::false_type>>{}...);
            },
            next_derivatives);

        constexpr auto next_derivatives_filtered =
            filter(next_derivatives, flags_next_derivatives);

        static_assert(size(next_derivatives_filtered));

        double const this_val_derivative = get_differential_operator_value(
            current_node_der_loc, current_node_der, it, ia, bt, ba);

        constexpr auto bt_free =
            free_on_buffer(current_node_der_loc, current_node_der, bt);

        constexpr auto next_derivatives_size =
            std::tuple_size_v<decltype(next_derivatives_filtered)>;

        std::array<double, next_derivatives_size> next_derivatives_values;
        next_derivatives_values.fill(this_val_derivative);

        constexpr auto multinomial_sequences = MultinomialSequences<2, pow>();
        constexpr auto multinomial_sequences_filtered =
            filter(multinomial_sequences, flags_next_derivatives);
        calc_add(multinomial_sequences_filtered, next_derivatives_values);

        constexpr auto locations =
            locate_new_vals(next_derivatives_filtered, it, bt_free);

        constexpr auto bt_new =
            update_buffer_types(next_derivatives_filtered, locations, bt_free);

        write_results(next_derivatives_filtered, next_derivatives_values,
                      locations, bt_new, ba, it, ia);

        constexpr auto flags_only_new = std::apply(
            [](auto... location) {
                return std::tuple_cat(
                    std::conditional_t<
                        std::is_same_v<on_buffer_new_t, decltype(location)>,
                        std::tuple<std::true_type>,
                        std::tuple<std::false_type>>{}...);
            },
            locations);

        constexpr auto next_derivatives_new =
            filter(next_derivatives_filtered, flags_only_new);

        constexpr auto dnn_new =
            merge_sorted(next_derivatives_new, dnn, NodesValue{});

        return treat_derivative_nodes_sum<MaxOrder, pow>(
            pn, tail(dnl), tail(dn), ct, it, ia, bt_new, ba, dnn_new, dnin);
    } else {
        return std::make_tuple(bt, dnn);
    }
}

template <class PrimalSubNode1, class PrimalSubNode2,
          class DerivativeNodeLocations, class DerivativeNodes, class CalcTree,
          class InterfaceTypes, class InterfaceArray, class BufferTypes,
          class BufferArray, std::size_t MaxOrder, class DerivativeNodeInputs>
auto treat_derivative_nodes(add_t<PrimalSubNode1, PrimalSubNode2> pn,
                            DerivativeNodeLocations dnl, DerivativeNodes dn,
                            CalcTree const &ct, InterfaceTypes it,
                            InterfaceArray &ia, BufferTypes bt, BufferArray &ba,
                            std::array<double, MaxOrder> & /* ua */,
                            DerivativeNodeInputs dnin) {
    // using PrimalNode = add_t<PrimalSubNode1, PrimalSubNode2>;

    // todo MaxOrder is never used? remove it?
    return treat_derivative_nodes_sum<MaxOrder, 0>(pn, dnl, dn, ct, it, ia, bt,
                                                   ba, std::tuple<>{}, dnin);
}

template <std::size_t PrevOrder, template <class> class Univariate,
          class PrimalSubNode, class DerivativeNodeLocations,
          class DerivativeNodes, class CalcTree, class InterfaceTypes,
          class InterfaceArray, class BufferTypes, class BufferArray,
          std::size_t MaxOrder, class DerivativeNodeNew,
          class DerivativeNodeInputs>
auto treat_derivative_nodes_univariate(
    Univariate<PrimalSubNode> pn, DerivativeNodeLocations dnl,
    DerivativeNodes dn, CalcTree const &ct, InterfaceTypes it,
    InterfaceArray &ia, BufferTypes bt, BufferArray &ba,
    std::array<double, MaxOrder> const &ua,
    std::array<double, MaxOrder> &ua_elevated, DerivativeNodeNew dnn,
    DerivativeNodeInputs dnin) {

    // std::cout << "it" << std::endl;
    // std::cout << type_name2<decltype(it)>() << std::endl;

    if constexpr (size(dnl)) {
        constexpr auto current_node_der_loc = head(dnl);
        constexpr auto current_node_der = head(dn);
        constexpr auto current_derivative_subnode = head(current_node_der);
        constexpr auto current_derivative_subnode_rest = tail(current_node_der);

        constexpr auto pow = get_power(current_derivative_subnode);
        // std::cout << pow << std::endl;
        // if (pow > 1) {
        //     std::cout << pow << std::endl;
        // }

        constexpr auto expansion_types =
            expand_univariate<PrimalSubNode, pow, MaxOrder>();
        std::cout << "expansion_types" << std::endl;
        std::cout << type_name2<decltype(expansion_types)>() << std::endl;

        using NodesValue = CalcTree::ValuesTupleInverse;

        constexpr auto next_derivatives = std::apply(
            [current_derivative_subnode_rest](const auto... singletype) {
                return std::make_tuple(multiply_differential_operator(
                    singletype, current_derivative_subnode_rest,
                    NodesValue{})...);
            },
            expansion_types);
        std::cout << "next_derivatives" << std::endl;
        std::cout << type_name2<decltype(next_derivatives)>() << std::endl;

        std::cout << "derivs inputs" << std::endl;
        std::cout << type_name2<decltype(dnin)>() << std::endl;

        constexpr auto flags_next_derivatives = std::apply(
            [dnin](auto... next_derivative) {
                return std::tuple_cat(
                    std::conditional_t<monomial_included(next_derivative, dnin),
                                       std::tuple<std::true_type>,
                                       std::tuple<std::false_type>>{}...);
            },
            next_derivatives);
        std::cout << "flags_next_derivatives" << std::endl;
        std::cout << type_name2<decltype(flags_next_derivatives)>()
                  << std::endl;

        constexpr auto next_derivatives_filtered =
            filter(next_derivatives, flags_next_derivatives);
        std::cout << "next_derivatives_filtered" << std::endl;
        std::cout << type_name2<decltype(next_derivatives_filtered)>()
                  << std::endl;
        static_assert(size(next_derivatives_filtered));

        constexpr auto next_derivatives_size =
            std::tuple_size_v<decltype(next_derivatives_filtered)>;

        std::array<double, next_derivatives_size> next_derivatives_values{0};

        if constexpr (pow > PrevOrder) {
            elevate_univariate<PrevOrder, pow>(ua, ua_elevated);
        }

        copy_filtered(ua_elevated, next_derivatives_values,
                      flags_next_derivatives);

        double const this_val_derivative = get_differential_operator_value(
            current_node_der_loc, current_node_der, it, ia, bt, ba);
        constexpr auto bt_free =
            free_on_buffer(current_node_der_loc, current_node_der, bt);

        for (std::size_t i = 0; i < next_derivatives_size; i++) {
            next_derivatives_values[i] *= this_val_derivative;
        }

        constexpr auto locations =
            locate_new_vals(next_derivatives_filtered, it, bt_free);
        std::cout << "locations" << std::endl;
        std::cout << type_name2<decltype(locations)>() << std::endl;

        constexpr auto bt_new =
            update_buffer_types(next_derivatives_filtered, locations, bt_free);
        std::cout << "bt_new" << std::endl;
        std::cout << type_name2<decltype(bt_new)>() << std::endl;

        write_results(next_derivatives_filtered, next_derivatives_values,
                      locations, bt_new, ba, it, ia);

        constexpr auto flags_only_new = std::apply(
            [](auto... location) {
                return std::tuple_cat(
                    std::conditional_t<
                        std::is_same_v<on_buffer_new_t, decltype(location)>,
                        std::tuple<std::true_type>,
                        std::tuple<std::false_type>>{}...);
            },
            locations);

        std::cout << "next_derivatives" << std::endl;
        std::cout << type_name2<decltype(next_derivatives)>() << std::endl;
        std::cout << "flags_only_new" << std::endl;
        std::cout << type_name2<decltype(flags_only_new)>() << std::endl;
        constexpr auto next_derivatives_new =
            filter(next_derivatives_filtered, flags_only_new);
        // constexpr auto next_derivatives_new = std::tuple<>{};

        constexpr auto dnn_new =
            merge_sorted(next_derivatives_new, dnn, NodesValue{});

        return treat_derivative_nodes_univariate<pow>(
            pn, tail(dnl), tail(dn), ct, it, ia, bt_new, ba, ua, ua_elevated,
            dnn_new, dnin);
    } else {
        return std::make_tuple(bt, dnn);
    }
}

template <template <class> class Univariate, class PrimalSubNode,
          class DerivativeNodeLocations, class DerivativeNodes, class CalcTree,
          class InterfaceTypes, class InterfaceArray, class BufferTypes,
          class BufferArray, std::size_t MaxOrder, class DerivativeNodeInputs>
auto treat_derivative_nodes(Univariate<PrimalSubNode> pn,
                            DerivativeNodeLocations dnl, DerivativeNodes dn,
                            CalcTree const &ct, InterfaceTypes it,
                            InterfaceArray &ia, BufferTypes bt, BufferArray &ba,
                            std::array<double, MaxOrder> &ua,
                            DerivativeNodeInputs dnin) {
    using PrimalNode = Univariate<PrimalSubNode>;

    PrimalNode::template d<MaxOrder>(ct.get(PrimalNode{}),
                                     ct.get(PrimalSubNode{}), ua);

    derivative_to_taylor(ua);

    std::array<double, MaxOrder> ua_elevated;

    return treat_derivative_nodes_univariate<0>(
        pn, dnl, dn, ct, it, ia, bt, ba, ua, ua_elevated, std::tuple<>{}, dnin);
}

template <class PrimalNode, class DerivativeNodeLocation, class DerivativeNodes,
          class CalcTree, class InterfaceTypes, class InterfaceArray,
          class BufferTypes, class BufferArray, class UnivariateArray,
          class DerivativeNodeInputs>
auto treat_node(PrimalNode nd, DerivativeNodeLocation dnl, DerivativeNodes dn,
                CalcTree const &ct, InterfaceTypes it, InterfaceArray &ia,
                BufferTypes bt, BufferArray &ba, UnivariateArray &ua,
                DerivativeNodeInputs dnin) {

    constexpr auto flags_derivative_nodes = std::apply(
        [nd](auto... type) {
            return std::tuple_cat(
                std::conditional_t<detail::first_type_is(type, nd),
                                   std::tuple<std::true_type>,
                                   std::tuple<std::false_type>>{}...);
        },
        dn);

    constexpr auto separated_derivative_nodes_loc =
        separate(dnl, flags_derivative_nodes);

    constexpr auto separated_derivative_nodes =
        separate(dn, flags_derivative_nodes);

    auto return_pair = treat_derivative_nodes(
        nd, std::get<0>(separated_derivative_nodes_loc),
        std::get<0>(separated_derivative_nodes), ct, it, ia, bt, ba, ua, dnin);

    constexpr decltype(return_pair) pair_const;
    constexpr auto bt_new = std::get<0>(pair_const);
    constexpr auto dn_new = std::get<1>(pair_const);

    constexpr auto dn_remaining = std::get<1>(separated_derivative_nodes);

    using NodesValue = CalcTree::ValuesTupleInverse;
    constexpr auto dn_new_and_remaining =
        merge_sorted(dn_new, dn_remaining, NodesValue{});

    // todo better
    constexpr auto dnl_new =
        make_tuple_same<size(dn_new_and_remaining)>(on_buffer_t{});

    return std::make_tuple(bt_new, dnl_new, dn_new_and_remaining);
}

template <class PrimalNodesLoop, class DerivativeNodeLocation,
          class DerivativeNodes, class CalcTree, class InterfaceTypes,
          class InterfaceArray, class BufferTypes, class BufferArray,
          class UnivariateArray, class DerivativeNodeInputs>
void backpropagate_aux(PrimalNodesLoop ndl, DerivativeNodeLocation const dnl,
                       DerivativeNodes dn, CalcTree const &ct,
                       InterfaceTypes it, InterfaceArray &ia, BufferTypes bt,
                       BufferArray &ba, UnivariateArray &ua,
                       DerivativeNodeInputs dnin) {
    constexpr auto current_primal_node = head(ndl);
    if constexpr (!is_input(current_primal_node)) {
        auto res = treat_node(current_primal_node, dnl, dn, ct, it, ia, bt, ba,
                              ua, dnin);

        constexpr decltype(res) pair_const;
        constexpr auto bt_new = std::get<0>(pair_const);
        constexpr auto dnl_new = std::get<1>(pair_const);
        constexpr auto dn_new = std::get<2>(pair_const);
        // std::cout << type_name2<decltype(bt_new)>() << std::endl;
        // std::cout << type_name2<decltype(dnl)>() << std::endl;
        // std::cout << type_name2<decltype(dnl_new)>() << std::endl;
        // std::cout << type_name2<decltype(dn_new)>() << std::endl;

        backpropagate_aux(tail(ndl), dnl_new, dn_new, ct, it, ia, bt_new, ba,
                          ua, dnin);
    }
}

} // namespace adhoc4::detail

#endif // ADHOC4_BACKPROPAGATOR_AUX_HPP
