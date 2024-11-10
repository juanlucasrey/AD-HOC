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
#include <utility>

#define LOG_LEVEL 0
#define LOG_VALS false

#if LOG_LEVEL
#include "type_name.hpp"
#include <iostream>
#endif

#if LOG_VALS
#include <iostream>
#endif

namespace adhoc4::detail {

template <std::size_t Last, std::size_t Idx1, std::size_t Idx2,
          std::size_t MaxOrder, std::size_t N = 0, class DerivativeNodes,
          class DerivativeNodeInputsI, class DerivativeNodeNew, class CalcTree,
          class InterfaceArray, class BufferFlags, class BufferArray,
          std::size_t ThisMaxOrder, std::size_t CrossedSize, class CrossedTypes,
          class DerivativeNodeInputs>
inline auto
treat_nodes_mul(DerivativeNodes dn, DerivativeNodeInputsI dnin_i,
                DerivativeNodeNew dnn_i, CalcTree const &ct, InterfaceArray &ia,
                BufferFlags bf, BufferArray &ba,
                std::array<double, ThisMaxOrder> const &powers1,
                std::array<double, ThisMaxOrder> const &powers2,
                std::array<double, CrossedSize> const &powers_crossed,
                CrossedTypes xt, DerivativeNodeInputs dnin) {

    if constexpr (N == Last) {
        return std::make_tuple(bf, dnn_i);
    } else {
        using NodesValue = typename CalcTree::ValuesTupleInverse;
        constexpr auto current_dn =
            convert_to_diffop<NodesValue>(std::get<N>(dn));

        constexpr auto current_node_der = std::get<0>(current_dn);
        constexpr auto current_node_loc = std::get<1>(current_dn);

#if LOG_LEVEL
        std::cout << "treating derivative tree node" << std::endl;
        std::cout << type_name2<decltype(current_node_der)>() << std::endl;
#endif

        constexpr auto current_derivative_subnode_rest = tail(current_node_der);
        constexpr auto this_power = get_power(head(current_node_der));
        constexpr auto rest_power = power(current_derivative_subnode_rest);

        constexpr auto multinomial_sequences =
            TrinomialSequencesMult<this_power, MaxOrder - rest_power>();

        using NodesValue = typename CalcTree::ValuesTupleInverse;
        using PrimalSubNodeOrdered1 = std::tuple_element_t<Idx1, NodesValue>;
        using PrimalSubNodeOrdered2 = std::tuple_element_t<Idx2, NodesValue>;
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

        constexpr auto next_derivatives_idx =
            convert_to_index_many2<NodesValue>(next_derivatives);

        constexpr auto next_derivatives_filtered =
            filter(next_derivatives_idx, flags_next_derivatives);

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

        calc_mul(next_derivatives_values, powers1, powers2, powers_crossed, xt,
                 multinomial_sequences_filtered);

        constexpr auto bf_new_pair =
            locate_new_vals_update_buffer_types<NodesValue, Last>(
                next_derivatives_filtered, bf_free, dn, dnn_i, dnin_i);
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
            merge_sorted(next_derivatives_new_with_pos, dnn_i, NodesValue{});

        return treat_nodes_mul<Last, Idx1, Idx2, MaxOrder, N + 1>(
            dn, dnin_i, dnn_new, ct, ia, bf_new, ba, powers1, powers2,
            powers_crossed, xt, dnin);
    }
}

template <std::size_t Idx, class PrimalSubNodeConst, std::size_t Last,
          std::size_t N = 0, class DerivativeNodes, class DerivativeNodeInputsI,
          class DerivativeNodeNew, class CalcTree, class InterfaceArray,
          class BufferFlags, class BufferArray, class DerivativeNodeInputs>
inline auto treat_nodes_mul_const(DerivativeNodes dn,
                                  DerivativeNodeInputsI dnin_i,
                                  DerivativeNodeNew dnn_i, CalcTree const &ct,
                                  InterfaceArray &ia, BufferFlags bf,
                                  BufferArray &ba, DerivativeNodeInputs dnin) {

    if constexpr (N == Last) {
        return std::make_tuple(bf, dnn_i);
    } else {
        using NodesValue = typename CalcTree::ValuesTupleInverse;
        constexpr auto current_dn =
            convert_to_diffop<NodesValue>(std::get<N>(dn));

        constexpr auto current_node_der = std::get<0>(current_dn);
        constexpr auto current_node_loc = std::get<1>(current_dn);

#if LOG_LEVEL
        std::cout << "treating derivative tree node" << std::endl;
        std::cout << type_name2<decltype(current_node_der)>() << std::endl;
#endif

        constexpr auto current_derivative_subnode_rest = tail(current_node_der);
        constexpr auto this_power = get_power(head(current_node_der));

        using NodesValue = typename CalcTree::ValuesTupleInverse;
        using PrimalSubNode = std::tuple_element_t<Idx, NodesValue>;
        constexpr auto next_derivatives_filtered = std::make_tuple(
            multiply_ordered(d<this_power>(PrimalSubNode{}),
                             current_derivative_subnode_rest, NodesValue{}));

        constexpr auto next_derivatives_filtered_idx =
            convert_to_index_many2<NodesValue>(next_derivatives_filtered);

        double const this_val_derivative =
            get_differential_operator_value(current_node_loc, ia, ba);

        constexpr auto bf_free = free_on_buffer(current_node_loc, bf);

        constexpr double mult = pow<this_power>(PrimalSubNodeConst::v());
        std::array<double, 1> next_derivatives_values{this_val_derivative *
                                                      mult};

        constexpr auto bf_new_pair =
            locate_new_vals_update_buffer_types<NodesValue, Last>(
                next_derivatives_filtered_idx, bf_free, dn, dnn_i, dnin_i);
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
            filter(next_derivatives_filtered_idx, flags_only_new);

        constexpr auto locations_new = filter(locations, flags_only_new);
        constexpr auto next_derivatives_new_with_pos =
            add_position(next_derivatives_new, locations_new);

        constexpr auto dnn_new =
            merge_sorted(next_derivatives_new_with_pos, dnn_i, NodesValue{});

        return treat_nodes_mul_const<Idx, PrimalSubNodeConst, Last, N + 1>(
            dn, dnin_i, dnn_new, ct, ia, bf_new, ba, dnin);
    }
}

template <std::size_t Last, class PrimalSubNode1, class PrimalSubNode2,
          class DerivativeNodes, class DerivativeNodeInputsI, class CalcTree,
          class InterfaceArray, class BufferFlags, class BufferArray,
          class DerivativeNodeInputs>
inline auto
treat_nodes_specialized(mul_t<PrimalSubNode1, PrimalSubNode2> /* pn */,
                        DerivativeNodes dn, DerivativeNodeInputsI dnin_i,
                        CalcTree const &ct, InterfaceArray &ia, BufferFlags bf,
                        BufferArray &ba, DerivativeNodeInputs dnin) {
    using NodesValue = typename CalcTree::ValuesTupleInverse;

    constexpr auto is_const1 = is_constant(PrimalSubNode1{});
    constexpr auto is_const2 = is_constant(PrimalSubNode2{});
    // if the 2 variables are constants the multiplication should not be here
    static_assert(!(is_const1 && is_const2));
    if constexpr (is_const1) {
        constexpr auto idx2 = find<NodesValue, PrimalSubNode2>();
        return treat_nodes_mul_const<idx2, PrimalSubNode1, Last>(
            dn, dnin_i, std::tuple<>{}, ct, ia, bf, ba, dnin);
    } else if constexpr (is_const2) {
        constexpr auto idx1 = find<NodesValue, PrimalSubNode1>();
        return treat_nodes_mul_const<idx1, PrimalSubNode2, Last>(
            dn, dnin_i, std::tuple<>{}, ct, ia, bf, ba, dnin);
    } else {

        // this assumes the derivatives nodes are ordered. the first one should
        // have the highest power
        using FirstDerivativeNode = std::tuple_element_t<
            0,
            std::tuple_element_t<0, std::tuple_element_t<0, DerivativeNodes>>>;
        constexpr auto ThisMaxOrder = get_power(FirstDerivativeNode{});
        constexpr auto MaxOrder = detail::max_orders(dnin_i);

        constexpr auto valstypes =
            AllTrinomialSequencesStored<ThisMaxOrder, MaxOrder>();

        double const val1 = ct.get(PrimalSubNode1{});
        double const val2 = ct.get(PrimalSubNode2{});
        auto const powers_val1 = detail::powers<ThisMaxOrder>(val1);
        auto const powers_val2 = detail::powers<ThisMaxOrder>(val2);

        constexpr auto idx1 = find<NodesValue, PrimalSubNode1>();
        constexpr auto idx2 = find<NodesValue, PrimalSubNode2>();
        constexpr auto id1_less_than_id2 = static_cast<bool>(idx1 >= idx2);
        if constexpr (id1_less_than_id2) {

            auto const powers_crossed =
                crossed_powers(valstypes, powers_val2, powers_val1);

            return treat_nodes_mul<Last, idx2, idx1, MaxOrder>(
                dn, dnin_i, std::tuple<>{}, ct, ia, bf, ba, powers_val2,
                powers_val1, powers_crossed, valstypes, dnin);
        } else {
            auto const powers_crossed =
                crossed_powers(valstypes, powers_val1, powers_val2);

            return treat_nodes_mul<Last, idx1, idx2, MaxOrder>(
                dn, dnin_i, std::tuple<>{}, ct, ia, bf, ba, powers_val1,
                powers_val2, powers_crossed, valstypes, dnin);
        }
    }
}

template <std::size_t Last, std::size_t Idx1, std::size_t Idx2,
          bool positiveIdx1 = true, bool positiveIdx2 = true, std::size_t N = 0,
          class DerivativeNodes, class DerivativeNodeInputsI,
          class DerivativeNodeNew, class CalcTree, class InterfaceArray,
          class BufferFlags, class BufferArray, class DerivativeNodeInputs>
inline auto treat_nodes_add(DerivativeNodes dn, DerivativeNodeInputsI dnin_i,
                            DerivativeNodeNew dnn_i, CalcTree const &ct,
                            InterfaceArray &ia, BufferFlags bf, BufferArray &ba,
                            DerivativeNodeInputs dnin) {

    if constexpr (N == Last) {
        return std::make_tuple(bf, dnn_i);
    } else {
        using NodesValue = typename CalcTree::ValuesTupleInverse;
        constexpr auto current_dn =
            convert_to_diffop<NodesValue>(std::get<N>(dn));

        constexpr auto current_node_der = std::get<0>(current_dn);
        constexpr auto current_node_loc = std::get<1>(current_dn);

#if LOG_LEVEL
        std::cout << "treating derivative tree node" << std::endl;
        std::cout << type_name2<decltype(current_node_der)>() << std::endl;
#endif

        constexpr auto current_derivative_subnode_rest = tail(current_node_der);
        constexpr auto this_power = get_power(head(current_node_der));

        constexpr auto multinomial_sequences =
            MultinomialSequences<2, this_power>();

        using NodesValue = typename CalcTree::ValuesTupleInverse;
        using PrimalSubNodeOrdered1 = std::tuple_element_t<Idx1, NodesValue>;
        using PrimalSubNodeOrdered2 = std::tuple_element_t<Idx2, NodesValue>;
        constexpr auto diff_ops = std::make_tuple(d(PrimalSubNodeOrdered1{}),
                                                  d(PrimalSubNodeOrdered2{}));

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

        constexpr auto next_derivatives_idx =
            convert_to_index_many2<NodesValue>(next_derivatives);

        constexpr auto next_derivatives_filtered =
            filter(next_derivatives_idx, flags_next_derivatives);

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

        calc_add<positiveIdx1, positiveIdx2>(multinomial_sequences_filtered,
                                             next_derivatives_values);

        constexpr auto bf_new_pair =
            locate_new_vals_update_buffer_types<NodesValue, Last>(
                next_derivatives_filtered, bf_free, dn, dnn_i, dnin_i);
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
            merge_sorted(next_derivatives_new_with_pos, dnn_i, NodesValue{});

        return treat_nodes_add<Last, Idx1, Idx2, positiveIdx1, positiveIdx2,
                               N + 1>(dn, dnin_i, dnn_new, ct, ia, bf_new, ba,
                                      dnin);
    }
}

template <std::size_t Last, class PrimalSubNode1, class PrimalSubNode2,
          class DerivativeNodesI, class DerivativeNodeInputsI, class CalcTree,
          class InterfaceArray, class BufferFlags, class BufferArray,
          class DerivativeNodeInputs>
inline auto
treat_nodes_specialized(sub_t<PrimalSubNode1, PrimalSubNode2> /* pn */,
                        DerivativeNodesI dn, DerivativeNodeInputsI dnin_i,
                        CalcTree const &ct, InterfaceArray &ia, BufferFlags bf,
                        BufferArray &ba, DerivativeNodeInputs dnin) {
    using NodesValue = typename CalcTree::ValuesTupleInverse;
    constexpr auto idx1 = find<NodesValue, PrimalSubNode1>();
    constexpr auto idx2 = find<NodesValue, PrimalSubNode2>();
    constexpr auto id1_less_than_id2 = static_cast<bool>(idx1 >= idx2);
    if constexpr (id1_less_than_id2) {
        return treat_nodes_add<Last, idx2, idx1, false, true>(
            dn, dnin_i, std::tuple<>{}, ct, ia, bf, ba, dnin);
    } else {
        return treat_nodes_add<Last, idx1, idx2, true, false>(
            dn, dnin_i, std::tuple<>{}, ct, ia, bf, ba, dnin);
    }
}

template <std::size_t Last, class PrimalSubNode1, class PrimalSubNode2,
          class DerivativeNodesI, class DerivativeNodeInputsI, class CalcTree,
          class InterfaceArray, class BufferFlags, class BufferArray,
          class DerivativeNodeInputs>
inline auto
treat_nodes_specialized(add_t<PrimalSubNode1, PrimalSubNode2> /* pn */,
                        DerivativeNodesI dn, DerivativeNodeInputsI dnin_i,
                        CalcTree const &ct, InterfaceArray &ia, BufferFlags bf,
                        BufferArray &ba, DerivativeNodeInputs dnin) {
    using NodesValue = typename CalcTree::ValuesTupleInverse;
    constexpr auto idx1 = find<NodesValue, PrimalSubNode1>();
    constexpr auto idx2 = find<NodesValue, PrimalSubNode2>();
    constexpr auto id1_less_than_id2 = static_cast<bool>(idx1 >= idx2);
    if constexpr (id1_less_than_id2) {
        return treat_nodes_add<Last, idx2, idx1>(dn, dnin_i, std::tuple<>{}, ct,
                                                 ia, bf, ba, dnin);
    } else {
        return treat_nodes_add<Last, idx1, idx2>(dn, dnin_i, std::tuple<>{}, ct,
                                                 ia, bf, ba, dnin);
    }
}

template <std::size_t Last, std::size_t N, std::size_t Idx,
          std::size_t PrevOrder = 0, class DerivativeNodes,
          class DerivativeNodeInputsI, class DerivativeNodeNew, class CalcTree,
          class InterfaceArray, class BufferFlags, class BufferArray,
          std::size_t MaxOrder, class DerivativeNodeInputs>
inline auto treat_nodes_univariate(
    DerivativeNodes dn, DerivativeNodeInputsI dnin_i, DerivativeNodeNew dnn,
    CalcTree const &ct, InterfaceArray &ia, BufferFlags bf, BufferArray &ba,
    std::array<double, MaxOrder> const &ua,
    std::array<double, MaxOrder> &ua_elevated, DerivativeNodeInputs dnin) {
    if constexpr (N == 0) {
        return std::make_tuple(bf, dnn);
    } else {
        constexpr std::size_t currentN = N - 1;
        using NodesValue = typename CalcTree::ValuesTupleInverse;
        constexpr auto current_dn =
            convert_to_diffop<NodesValue>(std::get<currentN>(dn));

        constexpr auto current_node_der = std::get<0>(current_dn);
        constexpr auto current_node_loc = std::get<1>(current_dn);

#if LOG_LEVEL
        std::cout << "treating derivative tree node" << std::endl;
        std::cout << type_name2<decltype(current_node_der)>() << std::endl;
#endif

        constexpr auto current_derivative_subnode = head(current_node_der);
        constexpr auto current_derivative_subnode_rest = tail(current_node_der);

        constexpr auto this_power = get_power(current_derivative_subnode);
        constexpr auto rest_power = power(current_derivative_subnode_rest);

        using PrimalSubNode = std::tuple_element_t<Idx, NodesValue>;
        constexpr auto expansion_types =
            expand_univariate<PrimalSubNode, MaxOrder - rest_power,
                              this_power>();

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

        constexpr auto next_derivatives_idx =
            convert_to_index_many2<NodesValue>(next_derivatives);

        constexpr auto next_derivatives_filtered =
            filter(next_derivatives_idx, flags_next_derivatives);

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

        for (std::size_t i = 0; i < next_derivatives_size; i++) {
            next_derivatives_values[i] *= this_val_derivative;
        }

        constexpr auto bf_new_pair =
            locate_new_vals_update_buffer_types<NodesValue, Last>(
                next_derivatives_filtered, bf_free, dn, dnn, dnin_i);
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

        return treat_nodes_univariate<Last, currentN, Idx, this_power>(
            dn, dnin_i, dnn_new, ct, ia, bf_new, ba, ua, ua_elevated, dnin);
    }
}

template <std::size_t Last, template <class> class Univariate,
          class PrimalSubNode, class DerivativeNodes,
          class DerivativeNodeInputsI, class CalcTree, class InterfaceArray,
          class BufferFlags, class BufferArray, class DerivativeNodeInputs>
inline auto
treat_nodes_specialized(Univariate<PrimalSubNode> /* pn */, DerivativeNodes dn,
                        DerivativeNodeInputsI dnin_i, CalcTree const &ct,
                        InterfaceArray &ia, BufferFlags bf, BufferArray &ba,
                        DerivativeNodeInputs dnin) {

    using PrimalNode = Univariate<PrimalSubNode>;
    constexpr auto MaxOrder = detail::max_orders(dnin_i);
    std::array<double, MaxOrder> univariate_array{};
    PrimalNode::template d<MaxOrder>(ct.get(PrimalNode{}),
                                     ct.get(PrimalSubNode{}), univariate_array);
    derivative_to_taylor(univariate_array);
    std::array<double, MaxOrder> univariate_array_elevated;

    // we go over univariate derivatives in inverse lexicographic order.
    // why? because it makes sense to calculate the coefficients increasing in
    // power. the lexicographic order is decreasing in powers.
    using NodesValue = typename CalcTree::ValuesTupleInverse;
    constexpr auto idx = find<NodesValue, PrimalSubNode>();
    return treat_nodes_univariate<Last, Last, idx>(
        dn, dnin_i, std::tuple<>{}, ct, ia, bf, ba, univariate_array,
        univariate_array_elevated, dnin);
}

template <std::size_t CN, class PrimalNode, class DerivativeNodes,
          class DerivativeNodeInputsI, class CalcTree, class InterfaceArray,
          class BufferFlags, class BufferArray, class DerivativeNodeInputs>
inline auto treat_node(DerivativeNodes dn, DerivativeNodeInputsI dnin_i,
                       PrimalNode nd, CalcTree const &ct, InterfaceArray &ia,
                       BufferFlags bf, BufferArray &ba,
                       DerivativeNodeInputs dnin) {

    constexpr auto Last = find_first_type_not2<DerivativeNodes, CN>();

    auto res =
        treat_nodes_specialized<Last>(nd, dn, dnin_i, ct, ia, bf, ba, dnin);

    constexpr std::tuple_element_t<0, decltype(res)> bf_new;
    constexpr std::tuple_element_t<1, decltype(res)> dn_new;

    using NodesValue = typename CalcTree::ValuesTupleInverse;

    constexpr auto dn_remaining = sub_tuple<Last, size(dn) - 1>(dn);
    constexpr auto dn_new_and_remaining =
        merge_sorted(dn_new, dn_remaining, NodesValue{});

    return std::make_tuple(bf_new, dn_new_and_remaining);
}

template <std::size_t N = 0, class DerivativeNodes, class DerivativeNodeInputsI,
          class CalcTree, class InterfaceArray, class BufferFlags,
          class BufferArray, class DerivativeNodeInputs>
inline void backpropagate_aux(DerivativeNodes dn, DerivativeNodeInputsI dnin_i,
                              CalcTree const &ct, InterfaceArray &ia,
                              BufferFlags bf, BufferArray &ba,
                              DerivativeNodeInputs dnin) {

    using PrimalNodes = CalcTree::ValuesTupleInverse;
    constexpr auto current_primal_node = std::get<N>(PrimalNodes{});

    if constexpr (!is_input(current_primal_node)) {

#if LOG_LEVEL
        std::cout << "treating calc tree node" << std::endl;
        std::cout << type_name2<decltype(current_primal_node)>() << std::endl;
#endif

        auto res = treat_node<N>(dn, dnin_i, current_primal_node, ct, ia, bf,
                                 ba, dnin);

        constexpr std::tuple_element_t<0, decltype(res)> bf_new;
        constexpr std::tuple_element_t<1, decltype(res)> dn_new;

        backpropagate_aux<N + 1>(dn_new, dnin_i, ct, ia, bf_new, ba, dnin);
    }
}

} // namespace adhoc4::detail

#endif // ADHOC4_BACKPROPAGATOR_AUX_HPP
