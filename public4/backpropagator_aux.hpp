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
#include "combinatorics/pow.hpp"
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

#define LOG_LEVEL 0
#define LOG_VALS false

#if LOG_LEVEL
#include "../tests4/type_name.hpp"
#include <iostream>
#endif

#if LOG_VALS
#include <iostream>
#endif

namespace adhoc4::detail {

template <std::size_t Last, std::size_t N = 0, class PrimalSubNodeOrdered1,
          class PrimalSubNodeOrdered2, class DerivativeNodes, class CalcTree,
          class InterfaceArray, class BufferFlags, class BufferArray,
          std::size_t MaxOrder, class DerivativeNodeNew,
          class DerivativeNodeInputs>
auto treat_nodes_mul(
    std::tuple<PrimalSubNodeOrdered1, PrimalSubNodeOrdered2> pn,
    DerivativeNodes dn, CalcTree const &ct, InterfaceArray &ia, BufferFlags bf,
    BufferArray &ba, std::array<double, MaxOrder> const &powers1,
    std::array<double, MaxOrder> const &powers2, DerivativeNodeNew dnn,
    DerivativeNodeInputs dnin) {

    if constexpr (N == Last) {
        return std::make_tuple(bf, dnn);
    } else {
        constexpr auto current_node_der = std::get<0>(std::get<N>(dn));
        constexpr auto current_node_loc = std::get<1>(std::get<N>(dn));

#if LOG_LEVEL
        std::cout << "treating derivative tree node" << std::endl;
        std::cout << type_name2<decltype(current_node_der)>() << std::endl;
#endif

        constexpr auto current_derivative_subnode_rest = tail(current_node_der);
        constexpr auto this_power = get_power(head(current_node_der));
        constexpr auto rest_power = power(current_derivative_subnode_rest);

        constexpr auto multinomial_sequences =
            TrinomialSequencesMult<this_power, MaxOrder - rest_power>();

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

#if LOG_LEVEL
        std::cout << "next_derivatives_filtered" << std::endl;
        std::cout << type_name2<decltype(next_derivatives_filtered)>()
                  << std::endl;
#endif

        constexpr auto multinomial_sequences_filtered =
            filter(multinomial_sequences, flags_next_derivatives);

        static_assert(size(next_derivatives_filtered));

        double const this_val_derivative =
            get_differential_operator_value(current_node_loc, ia, ba);

        constexpr auto bf_free = free_on_buffer(current_node_loc, bf);

        constexpr auto next_derivatives_size =
            std::tuple_size_v<decltype(next_derivatives_filtered)>;

        std::array<double, next_derivatives_size> next_derivatives_values;
        next_derivatives_values.fill(this_val_derivative);

        calc_mul(next_derivatives_values, powers1, powers2,
                 multinomial_sequences_filtered);

#if LOG_LEVEL
        std::cout << "locations" << std::endl;
        std::cout << type_name2<decltype(locations)>() << std::endl;
#endif

        constexpr auto bf_new_pair = locate_new_vals_update_buffer_types<Last>(
            next_derivatives_filtered, bf_free, dn, dnn, dnin);
        constexpr auto bf_new = std::get<0>(bf_new_pair);
        constexpr auto locations = std::get<1>(bf_new_pair);

        write_results(next_derivatives_values, locations, ba, ia);

#if LOG_VALS
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
#endif

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

#if LOG_LEVEL
        std::cout << "dnn_new" << std::endl;
        std::cout << type_name2<decltype(dnn_new)>() << std::endl;
#endif

        return treat_nodes_mul<Last, N + 1>(pn, dn, ct, ia, bf_new, ba, powers1,
                                            powers2, dnn_new, dnin);
    }
}

template <std::size_t Last, std::size_t N = 0, class PrimalSubNode,
          class PrimalSubNodeConst, class DerivativeNodes, class CalcTree,
          class InterfaceArray, class BufferFlags, class BufferArray,
          class DerivativeNodeNew, class DerivativeNodeInputs>
auto treat_nodes_mul_const(std::tuple<PrimalSubNode, PrimalSubNodeConst> pn,
                           DerivativeNodes dn, CalcTree const &ct,
                           InterfaceArray &ia, BufferFlags bf, BufferArray &ba,
                           DerivativeNodeNew dnn, DerivativeNodeInputs dnin) {

    if constexpr (N == Last) {
        return std::make_tuple(bf, dnn);
    } else {
        constexpr auto current_node_der = std::get<0>(std::get<N>(dn));
        constexpr auto current_node_loc = std::get<1>(std::get<N>(dn));

#if LOG_LEVEL
        std::cout << "treating derivative tree node" << std::endl;
        std::cout << type_name2<decltype(current_node_der)>() << std::endl;
#endif

        constexpr auto current_derivative_subnode_rest = tail(current_node_der);
        constexpr auto this_power = get_power(head(current_node_der));

        using NodesValue = CalcTree::ValuesTupleInverse;

        constexpr auto next_derivatives_filtered = std::make_tuple(
            multiply_ordered(d<this_power>(PrimalSubNode{}),
                             current_derivative_subnode_rest, NodesValue{}));

        double const this_val_derivative =
            get_differential_operator_value(current_node_loc, ia, ba);

        constexpr auto bf_free = free_on_buffer(current_node_loc, bf);

        constexpr double mult = pow<this_power>(PrimalSubNodeConst::v());
        std::array<double, 1> next_derivatives_values{this_val_derivative *
                                                      mult};

        constexpr auto bf_new_pair = locate_new_vals_update_buffer_types<Last>(
            next_derivatives_filtered, bf_free, dn, dnn, dnin);
        constexpr auto bf_new = std::get<0>(bf_new_pair);
        constexpr auto locations = std::get<1>(bf_new_pair);

        write_results(next_derivatives_values, locations, ba, ia);

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

        return treat_nodes_mul_const<Last, N + 1>(pn, dn, ct, ia, bf_new, ba,
                                                  dnn_new, dnin);
    }
}

template <std::size_t Last, class PrimalSubNode1, class PrimalSubNode2,
          class DerivativeNodes, class CalcTree, class InterfaceArray,
          class BufferFlags, class BufferArray, class DerivativeNodeInputs>
auto treat_nodes_specialized(mul_t<PrimalSubNode1, PrimalSubNode2> /* pn */,
                             DerivativeNodes dn, CalcTree const &ct,
                             InterfaceArray &ia, BufferFlags bf,
                             BufferArray &ba, DerivativeNodeInputs dnin) {
    using NodesValue = CalcTree::ValuesTupleInverse;
    constexpr auto ordered_pair = detail::sort_pair(
        std::tuple<PrimalSubNode1, PrimalSubNode2>{}, NodesValue{});

    static_assert(!is_constant(std::get<0>(ordered_pair)));
    constexpr bool is_const = is_constant(std::get<1>(ordered_pair));

    if constexpr (is_const) {
        return treat_nodes_mul_const<Last>(ordered_pair, dn, ct, ia, bf, ba,
                                           std::tuple<>{}, dnin);

    } else {
        constexpr auto MaxOrder = detail::max_orders(dnin);
        double const val1 = ct.get(std::get<0>(ordered_pair));
        double const val2 = ct.get(std::get<1>(ordered_pair));
        auto const powers_val1 = detail::powers<MaxOrder>(val1);
        auto const powers_val2 = detail::powers<MaxOrder>(val2);

        return treat_nodes_mul<Last>(ordered_pair, dn, ct, ia, bf, ba,
                                     powers_val1, powers_val2, std::tuple<>{},
                                     dnin);
    }
}

template <std::size_t Last, std::size_t N = 0, class PrimalSubNodeOrdered1,
          class PrimalSubNodeOrdered2, class DerivativeNodes, class CalcTree,
          class InterfaceArray, class BufferFlags, class BufferArray,
          class DerivativeNodeNew, class DerivativeNodeInputs>
auto treat_nodes_add(
    std::tuple<PrimalSubNodeOrdered1, PrimalSubNodeOrdered2> pn,
    DerivativeNodes dn, CalcTree const &ct, InterfaceArray &ia, BufferFlags bf,
    BufferArray &ba, DerivativeNodeNew dnn, DerivativeNodeInputs dnin) {

    if constexpr (N == Last) {
        return std::make_tuple(bf, dnn);
    } else {
        constexpr auto current_node_der = std::get<0>(std::get<N>(dn));
        constexpr auto current_node_loc = std::get<1>(std::get<N>(dn));

#if LOG_LEVEL
        std::cout << "treating derivative tree node" << std::endl;
        std::cout << type_name2<decltype(current_node_der)>() << std::endl;
#endif

        constexpr auto current_derivative_subnode_rest = tail(current_node_der);
        constexpr auto this_power = get_power(head(current_node_der));

        constexpr auto multinomial_sequences =
            MultinomialSequences<2, this_power>();
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

#if LOG_LEVEL
        std::cout << "next_derivatives_filtered" << std::endl;
        std::cout << type_name2<decltype(next_derivatives_filtered)>()
                  << std::endl;
#endif

        constexpr auto multinomial_sequences_filtered =
            filter(multinomial_sequences, flags_next_derivatives);

        static_assert(size(next_derivatives_filtered));

        double const this_val_derivative =
            get_differential_operator_value(current_node_loc, ia, ba);

        constexpr auto bf_free = free_on_buffer(current_node_loc, bf);

        constexpr auto next_derivatives_size =
            std::tuple_size_v<decltype(next_derivatives_filtered)>;

        std::array<double, next_derivatives_size> next_derivatives_values;
        next_derivatives_values.fill(this_val_derivative);

        calc_add(multinomial_sequences_filtered, next_derivatives_values);

        constexpr auto bf_new_pair = locate_new_vals_update_buffer_types<Last>(
            next_derivatives_filtered, bf_free, dn, dnn, dnin);
        constexpr auto bf_new = std::get<0>(bf_new_pair);
        constexpr auto locations = std::get<1>(bf_new_pair);

        write_results(next_derivatives_values, locations, ba, ia);

#if LOG_VALS
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
#endif

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

        return treat_nodes_add<Last, N + 1>(pn, dn, ct, ia, bf_new, ba, dnn_new,
                                            dnin);
    }
}

template <std::size_t Last, class PrimalSubNode1, class PrimalSubNode2,
          class DerivativeNodes, class CalcTree, class InterfaceArray,
          class BufferFlags, class BufferArray, class DerivativeNodeInputs>
auto treat_nodes_specialized(add_t<PrimalSubNode1, PrimalSubNode2> /* pn */,
                             DerivativeNodes dn, CalcTree const &ct,
                             InterfaceArray &ia, BufferFlags bf,
                             BufferArray &ba, DerivativeNodeInputs dnin) {
    using NodesValue = CalcTree::ValuesTupleInverse;
    constexpr auto ordered_pair = detail::sort_pair(
        std::tuple<PrimalSubNode1, PrimalSubNode2>{}, NodesValue{});

    return treat_nodes_add<Last>(ordered_pair, dn, ct, ia, bf, ba,
                                 std::tuple<>{}, dnin);
}

template <std::size_t Last, std::size_t N, std::size_t PrevOrder = 0,
          template <class> class Univariate, class PrimalSubNode,
          class DerivativeNodes, class CalcTree, class InterfaceArray,
          class BufferFlags, class BufferArray, std::size_t MaxOrder,
          class DerivativeNodeNew, class DerivativeNodeInputs>
auto treat_nodes_univariate(Univariate<PrimalSubNode> pn, DerivativeNodes dn,
                            CalcTree const &ct, InterfaceArray &ia,
                            BufferFlags bf, BufferArray &ba,
                            std::array<double, MaxOrder> const &ua,
                            std::array<double, MaxOrder> &ua_elevated,
                            DerivativeNodeNew dnn, DerivativeNodeInputs dnin) {
    if constexpr (N == 0) {
        return std::make_tuple(bf, dnn);
    } else {
        constexpr std::size_t currentN = N - 1;
        constexpr auto current_node_der = std::get<0>(std::get<currentN>(dn));
        constexpr auto current_node_loc = std::get<1>(std::get<currentN>(dn));

#if LOG_LEVEL
        std::cout << "treating derivative tree node" << std::endl;
        std::cout << type_name2<decltype(current_node_der)>() << std::endl;
#endif

        constexpr auto current_derivative_subnode = head(current_node_der);
        constexpr auto current_derivative_subnode_rest = tail(current_node_der);

#if LOG_LEVEL
        std::cout << "current_derivative_subnode_rest" << std::endl;
        std::cout << type_name2<decltype(current_derivative_subnode_rest)>()
                  << std::endl;
#endif

        constexpr auto this_power = get_power(current_derivative_subnode);
        constexpr auto rest_power = power(current_derivative_subnode_rest);
        constexpr auto expansion_types =
            expand_univariate<PrimalSubNode, MaxOrder - rest_power,
                              this_power>();

#if LOG_LEVEL
        std::cout << "expansion_types" << std::endl;
        std::cout << type_name2<decltype(expansion_types)>() << std::endl;
#endif

        using NodesValue = CalcTree::ValuesTupleInverse;

        constexpr auto next_derivatives = std::apply(
            [current_derivative_subnode_rest](const auto... singletype) {
                return std::make_tuple(multiply_differential_operator(
                    singletype, current_derivative_subnode_rest,
                    NodesValue{})...);
            },
            expansion_types);

#if LOG_LEVEL
        std::cout << "next_derivatives" << std::endl;
        std::cout << type_name2<decltype(next_derivatives)>() << std::endl;
#endif

        constexpr auto flags_next_derivatives = std::apply(
            [dnin](auto... next_derivative) {
                return std::integer_sequence<
                    bool, monomial_included(next_derivative, dnin)...>{};
            },
            next_derivatives);

#if LOG_LEVEL
        std::cout << "flags_next_derivatives" << std::endl;
        std::cout << type_name2<decltype(flags_next_derivatives)>()
                  << std::endl;
#endif

        constexpr auto next_derivatives_filtered =
            filter(next_derivatives, flags_next_derivatives);

#if LOG_LEVEL
        std::cout << "next_derivatives_filtered" << std::endl;
        std::cout << type_name2<decltype(next_derivatives_filtered)>()
                  << std::endl;
#endif

        static_assert(size(next_derivatives_filtered));
        constexpr auto next_derivatives_size =
            std::tuple_size_v<decltype(next_derivatives_filtered)>;

        std::array<double, next_derivatives_size> next_derivatives_values;

        if constexpr (this_power > PrevOrder) {
            elevate_univariate<PrevOrder, this_power>(ua, ua_elevated);
        }

        copy_filtered_inverted2<this_power - 1, MaxOrder - rest_power>(
            ua_elevated, next_derivatives_values, flags_next_derivatives);

        double const this_val_derivative =
            get_differential_operator_value(current_node_loc, ia, ba);

        constexpr auto bf_free = free_on_buffer(current_node_loc, bf);

#if LOG_LEVEL
        std::cout << "bt_free" << std::endl;
        std::cout << type_name2<decltype(bt_free)>() << std::endl;
#endif

        for (std::size_t i = 0; i < next_derivatives_size; i++) {
            next_derivatives_values[i] *= this_val_derivative;
        }

        constexpr auto bf_new_pair = locate_new_vals_update_buffer_types<Last>(
            next_derivatives_filtered, bf_free, dn, dnn, dnin);
        constexpr auto bf_new = std::get<0>(bf_new_pair);
        constexpr auto locations = std::get<1>(bf_new_pair);

#if LOG_LEVEL
        std::cout << "bf_new" << std::endl;
        std::cout << type_name2<decltype(bf_new)>() << std::endl;
#endif

        write_results(next_derivatives_values, locations, ba, ia);

#if LOG_VALS
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
#endif

        constexpr auto flags_only_new = std::apply(
            [](auto... location) {
                return std::integer_sequence<
                    bool,
                    std::is_same_v<on_buffer_t, decltype(location.first)>...>{};
            },
            locations);

#if LOG_LEVEL
        std::cout << "flags_only_new" << std::endl;
        std::cout << type_name2<decltype(flags_only_new)>() << std::endl;
#endif

        constexpr auto next_derivatives_new =
            filter(next_derivatives_filtered, flags_only_new);

        constexpr auto locations_new = filter(locations, flags_only_new);

        constexpr auto next_derivatives_new_with_pos =
            add_position(next_derivatives_new, locations_new);

#if LOG_LEVEL
        std::cout << "next_derivatives_new" << std::endl;
        std::cout << type_name2<decltype(next_derivatives_new)>() << std::endl;
#endif

        constexpr auto dnn_new =
            merge_sorted(next_derivatives_new_with_pos, dnn, NodesValue{});

        return treat_nodes_univariate<Last, currentN, this_power>(
            pn, dn, ct, ia, bf_new, ba, ua, ua_elevated, dnn_new, dnin);
    }
}

template <std::size_t Last, template <class> class Univariate,
          class PrimalSubNode, class DerivativeNodes, class CalcTree,
          class InterfaceArray, class BufferFlags, class BufferArray,
          class DerivativeNodeInputs>
auto treat_nodes_specialized(Univariate<PrimalSubNode> pn, DerivativeNodes dn,
                             CalcTree const &ct, InterfaceArray &ia,
                             BufferFlags bf, BufferArray &ba,
                             DerivativeNodeInputs dnin) {

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
    return treat_nodes_univariate<Last, Last>(
        pn, dn, ct, ia, bf, ba, univariate_array, univariate_array_elevated,
        std::tuple<>{}, dnin);
}

template <class PrimalNode, class DerivativeNodes, class CalcTree,
          class InterfaceArray, class BufferFlags, class BufferArray,
          class DerivativeNodeInputs>
auto treat_node(PrimalNode nd, DerivativeNodes dn, CalcTree const &ct,
                InterfaceArray &ia, BufferFlags bf, BufferArray &ba,
                DerivativeNodeInputs dnin) {

    constexpr auto Last = find_first_type_not<DerivativeNodes, PrimalNode>();

    auto res = treat_nodes_specialized<Last>(nd, dn, ct, ia, bf, ba, dnin);

    constexpr std::tuple_element_t<0, decltype(res)> bf_new;
    constexpr std::tuple_element_t<1, decltype(res)> dn_new;

    using NodesValue = CalcTree::ValuesTupleInverse;

    constexpr auto dn_remaining = sub_tuple<Last, size(dn) - 1>(dn);
    constexpr auto dn_new_and_remaining =
        merge_sorted(dn_new, dn_remaining, NodesValue{});

    return std::make_tuple(bf_new, dn_new_and_remaining);
}

template <std::size_t N = 0, class DerivativeNodes, class CalcTree,
          class InterfaceArray, class BufferFlags, class BufferArray,
          class DerivativeNodeInputs>
void backpropagate_aux(DerivativeNodes dn, CalcTree const &ct,
                       InterfaceArray &ia, BufferFlags bf, BufferArray &ba,
                       DerivativeNodeInputs dnin) {

    using PrimalNodes = CalcTree::ValuesTupleInverse;
    constexpr auto current_primal_node = std::get<N>(PrimalNodes{});

    if constexpr (!is_input(current_primal_node)) {

#if LOG_LEVEL
        std::cout << "treating calc tree node" << std::endl;
        std::cout << type_name2<decltype(current_primal_node)>() << std::endl;
#endif

        auto res = treat_node(current_primal_node, dn, ct, ia, bf, ba, dnin);

        constexpr std::tuple_element_t<0, decltype(res)> bf_new;
        constexpr std::tuple_element_t<1, decltype(res)> dn_new;

#if LOG_LEVEL
        std::cout << "dn_new" << std::endl;
        std::cout << type_name2<decltype(dn_new)>() << std::endl;
#endif

        backpropagate_aux<N + 1>(dn_new, ct, ia, bf_new, ba, dnin);
    }
}

} // namespace adhoc4::detail

#endif // ADHOC4_BACKPROPAGATOR_AUX_HPP
