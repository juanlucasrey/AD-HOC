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

// #include "combinatorics/combinations.hpp"
// #include "combinatorics/multinomial_coefficient_index_sequence.hpp"
// #include "utils/bivariate.hpp"
// #include "utils/index_sequence.hpp"
#include "backpropagator_tools.hpp"
#include "combinatorics/trinomial.hpp"
#include "dependency.hpp"
#include "differential_operator.hpp"
#include "monomial_included.hpp"
#include "utils/tuple.hpp"
#include "utils/univariate.hpp"

#include <array>
#include <cstddef>
#include <tuple>
#include <type_traits>

#include "../tests4/type_name.hpp"
#include <iostream>

#define LOG_LEVEL 0
#define LOG_VALS 0

namespace adhoc4::detail {

template <std::size_t N = 0, class PrimalSubNodeOrdered1,
          class PrimalSubNodeOrdered2, class DerivativeNodeLocations,
          class DerivativeNodes, class CalcTree, class InterfaceTypes,
          class InterfaceArray, class BufferTypes, class BufferArray,
          std::size_t MaxOrder, class DerivativeNodeNew,
          class DerivativeNodeInputs>
auto treat_nodes_mul(
    std::tuple<PrimalSubNodeOrdered1, PrimalSubNodeOrdered2> pn,
    DerivativeNodeLocations dnl, DerivativeNodes dn, CalcTree const &ct,
    InterfaceTypes it, InterfaceArray &ia, BufferTypes bt, BufferArray &ba,
    std::array<double, MaxOrder> const &powers1,
    std::array<double, MaxOrder> const &powers2, DerivativeNodeNew dnn,
    DerivativeNodeInputs dnin) {

    if constexpr (N == std::tuple_size_v<DerivativeNodeLocations>) {
        return std::make_tuple(bt, dnn);
    } else {
        constexpr auto current_node_der = std::get<N>(dn);

        if constexpr (LOG_LEVEL > 0) {
            std::cout << "treating derivative tree node" << std::endl;
            std::cout << type_name2<decltype(current_node_der)>() << std::endl;
        }

        constexpr auto current_node_der_loc = std::get<N>(dnl);
        constexpr auto current_derivative_subnode_rest = tail(current_node_der);
        constexpr auto pow = get_power(head(current_node_der));

        constexpr auto multinomial_sequences =
            TrinomialSequencesMult<pow, MaxOrder>();

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
                return std::tuple_cat(
                    std::conditional_t<monomial_included(next_derivative, dnin),
                                       std::tuple<std::true_type>,
                                       std::tuple<std::false_type>>{}...);
            },
            next_derivatives);

        constexpr auto next_derivatives_filtered =
            filter(next_derivatives, flags_next_derivatives);

        if constexpr (LOG_LEVEL > 0) {
            std::cout << "next_derivatives_filtered" << std::endl;
            std::cout << type_name2<decltype(next_derivatives_filtered)>()
                      << std::endl;
        }

        constexpr auto multinomial_sequences_filtered =
            filter(multinomial_sequences, flags_next_derivatives);

        static_assert(size(next_derivatives_filtered));

        double const this_val_derivative = get_differential_operator_value(
            current_node_der_loc, current_node_der, it, ia, bt, ba);

        constexpr auto bt_free =
            free_on_buffer(current_node_der_loc, current_node_der, bt);

        constexpr auto next_derivatives_size =
            std::tuple_size_v<decltype(next_derivatives_filtered)>;

        std::array<double, next_derivatives_size> next_derivatives_values;
        next_derivatives_values.fill(this_val_derivative);

        calc_mul(next_derivatives_values, powers1, powers2,
                 multinomial_sequences_filtered);

        constexpr auto locations =
            locate_new_vals(next_derivatives_filtered, it, bt_free);

        constexpr auto bt_new =
            update_buffer_types(next_derivatives_filtered, locations, bt_free);

        write_results(next_derivatives_filtered, next_derivatives_values,
                      locations, bt_new, ba, it, ia);

        if constexpr (LOG_VALS > 0) {
            std::cout.precision(std::numeric_limits<double>::max_digits10);

            std::cout << "buffer" << std::endl;
            for (std::size_t i = 0; i < ba.size(); i++) {
                std::cout << ba[i] << ", ";
            }
            std::cout << std::endl;

            std::cout << "interface" << std::endl;
            for (std::size_t i = 0; i < ia.size(); i++) {
                std::cout << ia[i] << ", ";
            }
            std::cout << std::endl;
        }

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
            merge_sorted2(next_derivatives_new, dnn, NodesValue{});

        return treat_nodes_mul<N + 1>(pn, dnl, dn, ct, it, ia, bt_new, ba,
                                      powers1, powers2, dnn_new, dnin);
    }
}

template <class PrimalSubNode1, class PrimalSubNode2,
          class DerivativeNodeLocations, class DerivativeNodes, class CalcTree,
          class InterfaceTypes, class InterfaceArray, class BufferTypes,
          class BufferArray, class DerivativeNodeInputs>
auto treat_nodes_specialized(mul_t<PrimalSubNode1, PrimalSubNode2> /* pn */,
                             DerivativeNodeLocations dnl, DerivativeNodes dn,
                             CalcTree const &ct, InterfaceTypes it,
                             InterfaceArray &ia, BufferTypes bt,
                             BufferArray &ba, DerivativeNodeInputs dnin) {
    using NodesValue = CalcTree::ValuesTupleInverse;
    constexpr auto ordered_pair = detail::sort_pair(
        std::tuple<PrimalSubNode1, PrimalSubNode2>{}, NodesValue{});

    constexpr auto MaxOrder = detail::max_orders(dnin);

    double const val1 = ct.get(std::get<0>(ordered_pair));
    double const val2 = ct.get(std::get<1>(ordered_pair));
    auto const powers_val1 = detail::powers<MaxOrder>(val1);
    auto const powers_val2 = detail::powers<MaxOrder>(val2);

    return treat_nodes_mul(ordered_pair, dnl, dn, ct, it, ia, bt, ba,
                           powers_val1, powers_val2, std::tuple<>{}, dnin);
}

template <std::size_t N = 0, class PrimalSubNodeOrdered1,
          class PrimalSubNodeOrdered2, class DerivativeNodeLocations,
          class DerivativeNodes, class CalcTree, class InterfaceTypes,
          class InterfaceArray, class BufferTypes, class BufferArray,
          class DerivativeNodeNew, class DerivativeNodeInputs>
auto treat_nodes_add(
    std::tuple<PrimalSubNodeOrdered1, PrimalSubNodeOrdered2> pn,
    DerivativeNodeLocations dnl, DerivativeNodes dn, CalcTree const &ct,
    InterfaceTypes it, InterfaceArray &ia, BufferTypes bt, BufferArray &ba,
    DerivativeNodeNew dnn, DerivativeNodeInputs dnin) {

    if constexpr (N == std::tuple_size_v<DerivativeNodeLocations>) {
        return std::make_tuple(bt, dnn);
    } else {
        constexpr auto current_node_der = std::get<N>(dn);

        if constexpr (LOG_LEVEL > 0) {
            std::cout << "treating derivative tree node" << std::endl;
            std::cout << type_name2<decltype(current_node_der)>() << std::endl;
        }

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

        if constexpr (LOG_LEVEL > 0) {
            std::cout << "next_derivatives_filtered" << std::endl;
            std::cout << type_name2<decltype(next_derivatives_filtered)>()
                      << std::endl;
        }

        constexpr auto multinomial_sequences_filtered =
            filter(multinomial_sequences, flags_next_derivatives);

        static_assert(size(next_derivatives_filtered));

        double const this_val_derivative = get_differential_operator_value(
            current_node_der_loc, current_node_der, it, ia, bt, ba);

        constexpr auto bt_free =
            free_on_buffer(current_node_der_loc, current_node_der, bt);

        constexpr auto next_derivatives_size =
            std::tuple_size_v<decltype(next_derivatives_filtered)>;

        std::array<double, next_derivatives_size> next_derivatives_values;
        next_derivatives_values.fill(this_val_derivative);

        calc_add(multinomial_sequences_filtered, next_derivatives_values);

        constexpr auto locations =
            locate_new_vals(next_derivatives_filtered, it, bt_free);

        constexpr auto bt_new =
            update_buffer_types(next_derivatives_filtered, locations, bt_free);

        write_results(next_derivatives_filtered, next_derivatives_values,
                      locations, bt_new, ba, it, ia);

        if constexpr (LOG_VALS > 0) {
            std::cout.precision(std::numeric_limits<double>::max_digits10);

            std::cout << "buffer" << std::endl;
            for (std::size_t i = 0; i < ba.size(); i++) {
                std::cout << ba[i] << ", ";
            }
            std::cout << std::endl;

            std::cout << "interface" << std::endl;
            for (std::size_t i = 0; i < ia.size(); i++) {
                std::cout << ia[i] << ", ";
            }
            std::cout << std::endl;
        }

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
            merge_sorted2(next_derivatives_new, dnn, NodesValue{});

        return treat_nodes_add<N + 1>(pn, dnl, dn, ct, it, ia, bt_new, ba,
                                      dnn_new, dnin);
    }
}

template <std::size_t N = 0, class PrimalSubNode1, class PrimalSubNode2,
          class DerivativeNodeLocations, class DerivativeNodes, class CalcTree,
          class InterfaceTypes, class InterfaceArray, class BufferTypes,
          class BufferArray, class DerivativeNodeInputs>
auto treat_nodes_specialized(add_t<PrimalSubNode1, PrimalSubNode2> pn,
                             DerivativeNodeLocations dnl, DerivativeNodes dn,
                             CalcTree const &ct, InterfaceTypes it,
                             InterfaceArray &ia, BufferTypes bt,
                             BufferArray &ba, DerivativeNodeInputs dnin) {
    using NodesValue = CalcTree::ValuesTupleInverse;
    constexpr auto ordered_pair = detail::sort_pair(
        std::tuple<PrimalSubNode1, PrimalSubNode2>{}, NodesValue{});

    return treat_nodes_add(ordered_pair, dnl, dn, ct, it, ia, bt, ba,
                           std::tuple<>{}, dnin);
}

template <std::size_t N, std::size_t PrevOrder = 0,
          template <class> class Univariate, class PrimalSubNode,
          class DerivativeNodeLocations, class DerivativeNodes, class CalcTree,
          class InterfaceTypes, class InterfaceArray, class BufferTypes,
          class BufferArray, std::size_t MaxOrder, class DerivativeNodeNew,
          class DerivativeNodeInputs>
auto treat_nodes_univariate(Univariate<PrimalSubNode> pn,
                            DerivativeNodeLocations dnl, DerivativeNodes dn,
                            CalcTree const &ct, InterfaceTypes it,
                            InterfaceArray &ia, BufferTypes bt, BufferArray &ba,
                            std::array<double, MaxOrder> const &ua,
                            std::array<double, MaxOrder> &ua_elevated,
                            DerivativeNodeNew dnn, DerivativeNodeInputs dnin) {
    if constexpr (N == 0) {
        return std::make_tuple(bt, dnn);
    } else {
        constexpr std::size_t currentN = N - 1;
        constexpr auto current_node_der = std::get<currentN>(dn);

        if constexpr (LOG_LEVEL > 0) {
            std::cout << "treating derivative tree node" << std::endl;
            std::cout << type_name2<decltype(current_node_der)>() << std::endl;
        }

        constexpr auto current_node_der_loc = std::get<currentN>(dnl);
        constexpr auto current_derivative_subnode = head(current_node_der);
        constexpr auto current_derivative_subnode_rest = tail(current_node_der);

        constexpr auto pow = get_power(current_derivative_subnode);
        constexpr auto expansion_types =
            expand_univariate<PrimalSubNode, MaxOrder, pow>();

        if constexpr (LOG_LEVEL > 1) {
            std::cout << "expansion_types" << std::endl;
            std::cout << type_name2<decltype(expansion_types)>() << std::endl;
        }

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
                return std::tuple_cat(
                    std::conditional_t<monomial_included(next_derivative, dnin),
                                       std::tuple<std::true_type>,
                                       std::tuple<std::false_type>>{}...);
            },
            next_derivatives);

        constexpr auto next_derivatives_filtered =
            filter(next_derivatives, flags_next_derivatives);

        if constexpr (LOG_LEVEL > 1) {
            std::cout << "next_derivatives_filtered" << std::endl;
            std::cout << type_name2<decltype(next_derivatives_filtered)>()
                      << std::endl;
        }

        static_assert(size(next_derivatives_filtered));
        constexpr auto next_derivatives_size =
            std::tuple_size_v<decltype(next_derivatives_filtered)>;

        std::array<double, next_derivatives_size> next_derivatives_values;

        if constexpr (pow > PrevOrder) {
            elevate_univariate<PrevOrder, pow>(ua, ua_elevated);
        }

        copy_filtered_inverted(ua_elevated, next_derivatives_values,
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

        constexpr auto bt_new =
            update_buffer_types(next_derivatives_filtered, locations, bt_free);

        write_results(next_derivatives_filtered, next_derivatives_values,
                      locations, bt_new, ba, it, ia);

        if constexpr (LOG_VALS > 0) {
            std::cout.precision(std::numeric_limits<double>::max_digits10);

            std::cout << "buffer" << std::endl;
            for (std::size_t i = 0; i < ba.size(); i++) {
                std::cout << ba[i] << ", ";
            }
            std::cout << std::endl;

            std::cout << "interface" << std::endl;
            for (std::size_t i = 0; i < ia.size(); i++) {
                std::cout << ia[i] << ", ";
            }
            std::cout << std::endl;
        }

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
            merge_sorted2(next_derivatives_new, dnn, NodesValue{});

        return treat_nodes_univariate<currentN, pow>(
            pn, dnl, dn, ct, it, ia, bt_new, ba, ua, ua_elevated, dnn_new,
            dnin);
    }
}

template <template <class> class Univariate, class PrimalSubNode,
          class DerivativeNodeLocations, class DerivativeNodes, class CalcTree,
          class InterfaceTypes, class InterfaceArray, class BufferTypes,
          class BufferArray, class DerivativeNodeInputs>
auto treat_nodes_specialized(Univariate<PrimalSubNode> pn,
                             DerivativeNodeLocations dnl, DerivativeNodes dn,
                             CalcTree const &ct, InterfaceTypes it,
                             InterfaceArray &ia, BufferTypes bt,
                             BufferArray &ba, DerivativeNodeInputs dnin) {

    using PrimalNode = Univariate<PrimalSubNode>;
    constexpr auto MaxOrder = detail::max_orders(dnin);
    std::array<double, MaxOrder> univariate_array{};
    PrimalNode::template d<MaxOrder>(ct.get(PrimalNode{}),
                                     ct.get(PrimalSubNode{}), univariate_array);
    derivative_to_taylor(univariate_array);
    std::array<double, MaxOrder> univariate_array_elevated;

    // we go over univariate derivatives in inverse lexicographic order.
    // why? because it makes sense to calculate the coefficients increasing in
    // power. the lexicographic order is decreasing in powers.
    constexpr std::size_t End = std::tuple_size_v<DerivativeNodeLocations>;
    return treat_nodes_univariate<End>(
        pn, dnl, dn, ct, it, ia, bt, ba, univariate_array,
        univariate_array_elevated, std::tuple<>{}, dnin);
}

template <class PrimalNode, class DerivativeNodeLocation, class DerivativeNodes,
          class CalcTree, class InterfaceTypes, class InterfaceArray,
          class BufferTypes, class BufferArray, class DerivativeNodeInputs>
auto treat_node(PrimalNode nd, DerivativeNodeLocation dnl, DerivativeNodes dn,
                CalcTree const &ct, InterfaceTypes it, InterfaceArray &ia,
                BufferTypes bt, BufferArray &ba, DerivativeNodeInputs dnin) {

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

    auto return_pair = treat_nodes_specialized(
        nd, std::get<0>(separated_derivative_nodes_loc),
        std::get<0>(separated_derivative_nodes), ct, it, ia, bt, ba, dnin);

    constexpr
        typename std::tuple_element<0, decltype(return_pair)>::type bt_new;
    constexpr
        typename std::tuple_element<1, decltype(return_pair)>::type dn_new;

    constexpr auto dn_remaining = std::get<1>(separated_derivative_nodes);

    using NodesValue = CalcTree::ValuesTupleInverse;

    constexpr auto dn_new_and_remaining =
        merge_sorted2(dn_new, dn_remaining, NodesValue{});

    // todo better for multiple output nodes
    constexpr auto dnl_new =
        make_tuple_same<size(dn_new_and_remaining)>(on_buffer_t{});

    return std::make_tuple(bt_new, dnl_new, dn_new_and_remaining);
}

template <std::size_t N = 0, class DerivativeNodeLocation,
          class DerivativeNodes, class CalcTree, class InterfaceTypes,
          class InterfaceArray, class BufferTypes, class BufferArray,
          class DerivativeNodeInputs>
void backpropagate_aux2(DerivativeNodeLocation const dnl, DerivativeNodes dn,
                        CalcTree const &ct, InterfaceTypes it,
                        InterfaceArray &ia, BufferTypes bt, BufferArray &ba,
                        DerivativeNodeInputs dnin) {

    using PrimalNodes = CalcTree::ValuesTupleInverse;
    constexpr auto current_primal_node = std::get<N>(PrimalNodes{});

    if constexpr (!is_input(current_primal_node)) {
        if constexpr (LOG_LEVEL > 0) {
            std::cout << "treating calc tree node" << std::endl;
            std::cout << type_name2<decltype(current_primal_node)>()
                      << std::endl;
        }

        auto res =
            treat_node(current_primal_node, dnl, dn, ct, it, ia, bt, ba, dnin);

        constexpr typename std::tuple_element<0, decltype(res)>::type bt_new;
        constexpr typename std::tuple_element<1, decltype(res)>::type dnl_new;
        constexpr typename std::tuple_element<2, decltype(res)>::type dn_new;

        backpropagate_aux2<N + 1>(dnl_new, dn_new, ct, it, ia, bt_new, ba,
                                  dnin);
    }
}

} // namespace adhoc4::detail

#endif // ADHOC4_BACKPROPAGATOR_AUX_HPP
