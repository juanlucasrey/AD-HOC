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

#ifndef ADHOC3_TAPE_BACKPROPAGATOR_HPP
#define ADHOC3_TAPE_BACKPROPAGATOR_HPP

#include <tape_buffer_size.hpp>
#include <tuple>
#include <tuple_utils.hpp>

#include <partition/multinomial_coefficient_index_sequence.hpp>

#include "../tests3/type_name.hpp"

#include <array>
#include <type_traits>

namespace adhoc3::detail {

template <class NodeDerivatives, class CalcTreeValue, class InterfaceTypes,
          class InterfaceArray, class BufferTypes, class BufferArray,
          class UnivariateType, class UnivariateArray>
void backpropagate_aux(NodeDerivatives nd, CalcTreeValue const &ct,
                       InterfaceTypes it, InterfaceArray &ia, BufferTypes bt,
                       BufferArray &ba, UnivariateType ut, UnivariateArray &ua);

class available_t {};

template <class NodesDerivativeRest, class InterfaceTypes, class DersActive,
          class CandidateNode>
constexpr auto conditional_node2(NodesDerivativeRest dernodes,
                                 InterfaceTypes dersin,
                                 CandidateNode candnode) {
    if constexpr (contains(dernodes, candnode) || contains(dersin, candnode)) {
        return std::make_tuple(candnode);
    } else {
        return std::tuple<>{};
    }
}

template <std::size_t I = 0, class IntegerSequence, class InputArray>
constexpr void mult_multinomial(IntegerSequence is, InputArray const &in,
                                double &out) {
    if constexpr (I < is.size()) {
        if constexpr (get<I>(is) == 1) {
            out *= in[I];
        } else if constexpr (get<I>(is) > 1) {
            out *= std::pow(in[I], get<I>(is));
        }

        mult_multinomial<I + 1>(is, in, out);
    }
}

template <std::size_t N = 0, class PartitionIntegerSequences,
          class CalculateFlags, class Array, class ArrayOut>
auto backpropagate_univariate(PartitionIntegerSequences const sequences,
                              CalculateFlags const calc_flags,
                              Array const &vals, ArrayOut &arrayout) {
    if constexpr (N < std::tuple_size_v<PartitionIntegerSequences>) {
        if constexpr (std::get<N>(calc_flags)) {
            constexpr auto current_sequence = std::get<N>(sequences);
            constexpr auto order = sum(current_sequence);
            arrayout[N] = BellCoeff(current_sequence) * vals[order - 1];
        }

        backpropagate_univariate<N + 1>(sequences, calc_flags, vals, arrayout);
    }
}

template <std::size_t N = 0, class PartitionIntegerSequences, class Array,
          class ArrayOut>
auto backpropagate_multivariate(PartitionIntegerSequences const sequences,
                                Array const &vals, ArrayOut &arrayout) {
    if constexpr (N < std::tuple_size_v<PartitionIntegerSequences>) {
        constexpr auto current_sequence = std::get<N>(sequences);
        arrayout[N] = MultinomialCoeff2(current_sequence);
        mult_multinomial(current_sequence, vals, arrayout[N]);

        backpropagate_multivariate<N + 1>(sequences, vals, arrayout);
    }
}

template <class CurrentNodesDerivatives, class InterfaceTypes,
          class InterfaceArray, class BufferTypes, class BufferArray>
auto get_differential_operator_value(
    CurrentNodesDerivatives current_derivative_node, InterfaceTypes it,
    InterfaceArray const &ia, BufferTypes bt, BufferArray const &ba) -> double {

    constexpr bool is_on_interface = contains(it, current_derivative_node);

    if constexpr (is_on_interface) {
        constexpr auto idx = get_idx(it, current_derivative_node);
        return ia[idx];
    } else {
        constexpr bool is_on_buffer = contains(bt, current_derivative_node);
        static_assert(is_on_buffer);
        constexpr auto idx = get_idx(bt, current_derivative_node);
        return ba[idx];
    }
}

class on_buffer_add_t {};
class on_buffer_new_t {};
class on_interface_t {};

template <std::size_t N = 0, class TypesToPlace, class InterfaceTypes,
          class BufferTypes, class Output>
constexpr auto locate_new_vals_aux(TypesToPlace derivative_nodes,
                                   InterfaceTypes it, BufferTypes bt,
                                   Output out) {
    if constexpr (N < std::tuple_size_v<TypesToPlace>) {
        constexpr auto current_derivative_node = std::get<N>(derivative_nodes);
        constexpr bool is_on_interface = contains(it, current_derivative_node);
        constexpr bool is_on_buffer = contains(bt, current_derivative_node);

        using this_result = std::conditional_t<
            is_on_interface, on_interface_t,
            std::conditional_t<is_on_buffer, on_buffer_add_t, on_buffer_new_t>>;

        return locate_new_vals_aux<N + 1>(
            derivative_nodes, it, bt,
            std::tuple_cat(out, std::make_tuple(this_result{})));
    } else {
        return out;
    }
}

template <class TypesToPlace, class InterfaceTypes, class BufferTypes>
constexpr auto locate_new_vals(TypesToPlace derivative_nodes, InterfaceTypes it,
                               BufferTypes bt) {
    return locate_new_vals_aux(derivative_nodes, it, bt, std::tuple<>{});
}

template <class Node, class InterfaceNodesTuple, class BufferNodesTuple>
constexpr auto locate_val(Node node, InterfaceNodesTuple interface_nodes_tuple,
                          BufferNodesTuple buffer_nodes_tuple) {

    constexpr bool is_on_interface = contains(interface_nodes_tuple, node);
    if constexpr (is_on_interface) {
        return on_interface_t{};
    } else {
        constexpr bool is_on_buffer = contains(buffer_nodes_tuple, node);
        if constexpr (is_on_buffer) {
            return on_buffer_add_t{};
        } else {
            return on_buffer_new_t{};
        }
    }
}

template <std::size_t N = 0, class TypesToPlace, class LocationIndicators,
          class BufferTypes>
constexpr auto update_buffer_types(TypesToPlace derivative_nodes,
                                   LocationIndicators location_indicators,
                                   BufferTypes bt) {
    if constexpr (N < std::tuple_size_v<LocationIndicators>) {
        constexpr auto current_indicator = std::get<N>(location_indicators);

        if constexpr (std::is_same_v<const on_buffer_new_t,
                                     decltype(current_indicator)>) {

            constexpr auto current_derivative_node =
                std::get<N>(derivative_nodes);
            constexpr auto new_buffer_type =
                replace(bt, available_t{}, current_derivative_node);
            return update_buffer_types<N + 1>(
                derivative_nodes, location_indicators, new_buffer_type);
        } else {
            return update_buffer_types<N + 1>(derivative_nodes,
                                              location_indicators, bt);
        }
    } else {
        return bt;
    }
}

template <std::size_t N = 0, class ResultsTypes, class ResultsArray,
          class LocationIndicators, class BufferTypes, class BufferArray,
          class InterfaceTypes, class InterfaceArray>
auto write_results(ResultsTypes derivatives_nodes,
                   ResultsArray const &results_array,
                   LocationIndicators location_indicators,
                   BufferTypes buffer_types_updated, BufferArray &buffer_array,
                   InterfaceTypes interface_types,
                   InterfaceArray &interface_array) {
    if constexpr (N < std::tuple_size_v<LocationIndicators>) {

        constexpr auto current_indicator = std::get<N>(location_indicators);
        constexpr auto current_derivative_node = std::get<N>(derivatives_nodes);

        if constexpr (std::is_same_v<const on_interface_t,
                                     decltype(current_indicator)>) {
            constexpr auto idx =
                get_idx(interface_types, current_derivative_node);

            // we always add on the interface
            interface_array[idx] += results_array[N];
        } else {
            constexpr auto idx =
                get_idx(buffer_types_updated, current_derivative_node);

            if constexpr (std::is_same_v<const on_buffer_new_t,
                                         decltype(current_indicator)>) {
                // we place a new value so we override
                buffer_array[idx] = results_array[N];
            } else {
                // there is already a value so we add
                buffer_array[idx] += results_array[N];
            }
        }

        write_results<N + 1>(derivatives_nodes, results_array,
                             location_indicators, buffer_types_updated,
                             buffer_array, interface_types, interface_array);
    }
}

template <class Node, class LocationIndicator, class BufferTypes,
          class BufferArray, class InterfaceTypes, class InterfaceArray>
auto write_result(Node node, double const &node_value,
                  LocationIndicator /* location_indicator */,
                  BufferTypes buffer_types_updated, BufferArray &buffer_array,
                  InterfaceTypes interface_types,
                  InterfaceArray &interface_array) {

    if constexpr (std::is_same_v<on_interface_t, LocationIndicator>) {
        constexpr auto idx = get_idx(interface_types, node);

        // we always add on the interface
        interface_array[idx] += node_value;
    } else {
        // if constexpr (contains(buffer_types_updated, node)) {
        constexpr auto idx = get_idx(buffer_types_updated, node);

        if constexpr (std::is_same_v<on_buffer_new_t, LocationIndicator>) {
            // we place a new value so we override
            buffer_array[idx] = node_value;
        } else {
            // there is already a value so we add
            buffer_array[idx] += node_value;
        }
    }
}

template <std::size_t Power, std::size_t Order,
          template <class> class Univariate, class NodeDerivative,
          class FirstNodesDerivativeRest, class CurrentNodesDerivatives,
          class NextNodesDerivatives, class CalcTreeValue, class InterfaceTypes,
          class InterfaceArray, class BufferTypes, class BufferArray,
          class UnivariateType, class UnivariateArray>
void backpropagate_process(
    der2::p<Power, der2::d<Order, Univariate<NodeDerivative>>> /* nd */,
    FirstNodesDerivativeRest ndr,
    CurrentNodesDerivatives current_derivative_node,
    NextNodesDerivatives next_derivative_nodes, CalcTreeValue const &ct,
    InterfaceTypes it, InterfaceArray &ia, BufferTypes bt, BufferArray &ba,
    UnivariateType /* ut */, UnivariateArray &ua) {

    using OutId = Univariate<NodeDerivative>;
    using InId = NodeDerivative;

    constexpr bool is_repeated_univariate =
        std::is_same_v<OutId, UnivariateType>;
    if constexpr (!is_repeated_univariate) {
        // we need to calculate derivative values
        auto univ_ders =
            OutId::template d2<Order>(ct.val(OutId{}), ct.val(InId{}));

        std::copy(univ_ders.begin(), univ_ders.end(), ua.begin());
    }

    using NodesValue = CalcTreeValue::ValuesTupleInverse;

    constexpr auto univariate_size = partition_function(Order);
    std::array<double, univariate_size> univariate_expansion_values;
    univariate_expansion_values.fill(0);

    constexpr auto univariate_expansion_types =
        expand_univariate<InId, Order>(NodesValue{});

    constexpr auto multinomial_expansion_types =
        expand_multinomial<Power>(NodesValue{}, univariate_expansion_types);

    constexpr auto multinomial_expansion_types_full =
        multiply_ordered_tuple(NodesValue{}, multinomial_expansion_types, ndr);

    constexpr auto multinomial_calc_flags = std::apply(
        [next_derivative_nodes, it](auto... type) {
            return std::tuple_cat(
                std::conditional_t < contains(next_derivative_nodes, type) ||
                    contains(it, type),
                std::tuple<std::true_type>,
                std::tuple < std::false_type >> {}...);
        },
        multinomial_expansion_types_full);

    constexpr auto multinomial_sequences =
        MultinomialSequences<univariate_size, Power>();

    constexpr auto multinomial_sequences_filtered =
        filter(multinomial_sequences, multinomial_calc_flags);
    constexpr auto multinomial_sequences_filtered_size =
        std::tuple_size_v<decltype(multinomial_sequences_filtered)>;

    constexpr auto univariate_calc_flags =
        precedent_required(multinomial_sequences_filtered);

    constexpr auto partition_sequences = PartitionSequences<Order>();

    backpropagate_univariate(partition_sequences, univariate_calc_flags, ua,
                             univariate_expansion_values);

    std::array<double, multinomial_sequences_filtered_size>
        multinomial_expansion_values{0};
    multinomial_expansion_values.fill(0);

    backpropagate_multivariate(multinomial_sequences_filtered,
                               univariate_expansion_values,
                               multinomial_expansion_values);

    constexpr auto multinomial_expansion_types_full_filtered =
        filter(multinomial_expansion_types_full, multinomial_calc_flags);

    double const this_val_derivative = get_differential_operator_value(
        current_derivative_node, it, ia, bt, ba);
    constexpr auto bt_free =
        replace(bt, current_derivative_node, available_t{});

    for (std::size_t i = 0; i < multinomial_sequences_filtered_size; i++)
        multinomial_expansion_values[i] *= this_val_derivative;

    constexpr auto locations =
        locate_new_vals(multinomial_expansion_types_full_filtered, it, bt_free);

    constexpr auto new_bt = update_buffer_types(
        multinomial_expansion_types_full_filtered, locations, bt_free);

    write_results(multinomial_expansion_types_full_filtered,
                  multinomial_expansion_values, locations, new_bt, ba, it, ia);

    backpropagate_aux(next_derivative_nodes, ct, it, ia, new_bt, ba, OutId{},
                      ua);
}

template <std::size_t N = 0, class ResultsTypes, class LocationIndicators,
          class BufferTypes, class BufferArray, class InterfaceTypes,
          class InterfaceArray>
auto write_results2(ResultsTypes derivatives_nodes, double result,
                    LocationIndicators location_indicators,
                    BufferTypes buffer_types_updated, BufferArray &buffer_array,
                    InterfaceTypes interface_types,
                    InterfaceArray &interface_array) {
    if constexpr (N < std::tuple_size_v<LocationIndicators>) {

        constexpr auto current_indicator = std::get<N>(location_indicators);
        constexpr auto current_derivative_node = std::get<N>(derivatives_nodes);

        if constexpr (std::is_same_v<const on_interface_t,
                                     decltype(current_indicator)>) {
            constexpr auto idx =
                get_idx(interface_types, current_derivative_node);

            // we always add on the interface
            interface_array[idx] += result;
        } else {
            constexpr auto idx =
                get_idx(buffer_types_updated, current_derivative_node);

            if constexpr (std::is_same_v<const on_buffer_new_t,
                                         decltype(current_indicator)>) {
                // we place a new value so we override
                buffer_array[idx] = result;
            } else {
                // there is already a value so we add
                buffer_array[idx] += result;
            }
        }

        write_results2<N + 1>(derivatives_nodes, result, location_indicators,
                              buffer_types_updated, buffer_array,
                              interface_types, interface_array);
    }
}

template <std::size_t N = 0, class PartitionIntegerSequences, class ArrayOut>
auto multiply_multivariate(PartitionIntegerSequences const sequences,
                           ArrayOut &arrayout) {
    if constexpr (N < std::tuple_size_v<PartitionIntegerSequences>) {
        constexpr auto current_sequence = std::get<N>(sequences);
        if constexpr (MultinomialCoeff2(current_sequence) != 1) {
            arrayout[N] *= MultinomialCoeff2(current_sequence);
        }
        multiply_multivariate<N + 1>(sequences, arrayout);
    }
}

template <std::size_t Power, std::size_t Order, class NodeDerivative1,
          class NodeDerivative2, class FirstNodesDerivativeRest,
          class CurrentNodesDerivatives, class NextNodesDerivatives,
          class CalcTreeValue, class InterfaceTypes, class InterfaceArray,
          class BufferTypes, class BufferArray, class UnivariateType,
          class UnivariateArray>
void backpropagate_process(
    der2::p<Power, der2::d<Order, add_t<NodeDerivative1, NodeDerivative2>>> nd,
    FirstNodesDerivativeRest ndr,
    CurrentNodesDerivatives current_derivative_node,
    NextNodesDerivatives next_derivative_nodes, CalcTreeValue const &ct,
    InterfaceTypes it, InterfaceArray &ia, BufferTypes bt, BufferArray &ba,
    UnivariateType ut, UnivariateArray &ua) {

    constexpr auto current_node_loc =
        locate_val(current_derivative_node, it, bt);

    using NodesValue = CalcTreeValue::ValuesTupleInverse;
    constexpr auto new_subnodes = expand_tree_single(NodesValue{}, nd);
    constexpr auto new_subnodes_full =
        multiply_ordered_tuple(NodesValue{}, new_subnodes, ndr);

    constexpr auto calc_flags = std::apply(
        [next_derivative_nodes, it](auto... type) {
            return std::tuple_cat(
                std::conditional_t < contains(next_derivative_nodes, type) ||
                    contains(it, type),
                std::tuple<std::true_type>,
                std::tuple < std::false_type >> {}...);
        },
        new_subnodes_full);

    constexpr auto new_subnodes_full_filtered =
        filter(new_subnodes_full, calc_flags);

    constexpr auto bt_free =
        replace(bt, current_derivative_node, available_t{});
    constexpr auto new_subnodes_locations =
        locate_new_vals(new_subnodes_full_filtered, it, bt_free);

    constexpr auto new_bt = update_buffer_types(
        new_subnodes_full_filtered, new_subnodes_locations, bt_free);

    constexpr auto same_nodes =
        std::is_same_v<NodeDerivative1, NodeDerivative2>;

    double const this_val_derivative = get_differential_operator_value(
        current_derivative_node, it, ia, bt, ba);

    std::array<double, size(new_subnodes_full_filtered)>
        multinomial_expansion_values;

    if constexpr (is_constant_class2_v<NodeDerivative1>) {
        static_assert(size(new_subnodes_full_filtered) == 1);
        write_result(std::get<0>(new_subnodes_full_filtered),
                     this_val_derivative, std::get<0>(new_subnodes_locations),
                     new_bt, ba, it, ia);
    } else if constexpr (is_constant_class2_v<NodeDerivative2>) {
        static_assert(size(new_subnodes_full_filtered) == 1);
        write_result(std::get<0>(new_subnodes_full_filtered),
                     this_val_derivative, std::get<0>(new_subnodes_locations),
                     new_bt, ba, it, ia);
    } else if constexpr (same_nodes) {
        static_assert(size(new_subnodes_full_filtered) == 1);
        write_result(std::get<0>(new_subnodes_full_filtered),
                     std::pow(2, Power) * this_val_derivative,
                     std::get<0>(new_subnodes_locations), new_bt, ba, it, ia);
    } else {
        multinomial_expansion_values.fill(this_val_derivative);
        constexpr auto multinomial_sequences = MultinomialSequences<2, Power>();
        constexpr auto multinomial_sequences_filtered =
            filter(multinomial_sequences, calc_flags);
        multiply_multivariate(multinomial_sequences_filtered,
                              multinomial_expansion_values);

        write_results(new_subnodes_full_filtered, multinomial_expansion_values,
                      new_subnodes_locations, new_bt, ba, it, ia);
    }

    backpropagate_aux(next_derivative_nodes, ct, it, ia, new_bt, ba, ut, ua);
}

class minus_t {};
class plus_t {};

template <std::size_t N, std::size_t Total, bool Inverted, class TupleOut>
constexpr auto create_signs(TupleOut out) {
    if constexpr (N < Total) {
        if constexpr (N & 1U) {
            if constexpr (Inverted) {
                return create_signs<N + 1, Total, Inverted>(
                    std::tuple_cat(out, std::make_tuple(plus_t{})));

            } else {
                return create_signs<N + 1, Total, Inverted>(
                    std::tuple_cat(out, std::make_tuple(minus_t{})));
            }
        } else {
            if constexpr (Inverted) {
                return create_signs<N + 1, Total, Inverted>(
                    std::tuple_cat(out, std::make_tuple(minus_t{})));
            } else {
                return create_signs<N + 1, Total, Inverted>(
                    std::tuple_cat(out, std::make_tuple(plus_t{})));
            }
        }
    } else {
        return out;
    }
}

template <class Nodes, class NodeDerivative1, class NodeDerivative2,
          std::size_t Power>
constexpr auto signatures() {
    constexpr auto size = combinations(2, Power);

    constexpr auto id1_less_than_id2 =
        static_cast<bool>(get_idx_first2<NodeDerivative1>(Nodes{}) >=
                          get_idx_first2<NodeDerivative2>(Nodes{}));
    return create_signs<0, size, id1_less_than_id2>(std::tuple<>{});
}

template <std::size_t N = 0, class Array, class Signatures>
auto fill_with_signatures(Array &out, Signatures signs, double val) {
    if constexpr (N < std::tuple_size_v<Array>) {
        constexpr auto current_sign = std::get<N>(signs);
        if constexpr (std::is_same_v<decltype(current_sign), const plus_t>) {
            out[N] = val;
        } else {
            out[N] = -val;
        }
        fill_with_signatures<N + 1>(out, signs, val);
    }
}

template <std::size_t Power, std::size_t Order, class NodeDerivative1,
          class NodeDerivative2, class FirstNodesDerivativeRest,
          class CurrentNodesDerivatives, class NextNodesDerivatives,
          class CalcTreeValue, class InterfaceTypes, class InterfaceArray,
          class BufferTypes, class BufferArray, class UnivariateType,
          class UnivariateArray>
void backpropagate_process(
    der2::p<Power, der2::d<Order, sub_t<NodeDerivative1, NodeDerivative2>>> nd,
    FirstNodesDerivativeRest ndr,
    CurrentNodesDerivatives current_derivative_node,
    NextNodesDerivatives next_derivative_nodes, CalcTreeValue const &ct,
    InterfaceTypes it, InterfaceArray &ia, BufferTypes bt, BufferArray &ba,
    UnivariateType ut, UnivariateArray &ua) {

    constexpr auto current_node_loc =
        locate_val(current_derivative_node, it, bt);

    using NodesValue = CalcTreeValue::ValuesTupleInverse;
    constexpr auto new_subnodes = expand_tree_single(NodesValue{}, nd);
    constexpr auto new_subnodes_full =
        multiply_ordered_tuple(NodesValue{}, new_subnodes, ndr);

    constexpr auto calc_flags = std::apply(
        [next_derivative_nodes, it](auto... type) {
            return std::tuple_cat(
                std::conditional_t < contains(next_derivative_nodes, type) ||
                    contains(it, type),
                std::tuple<std::true_type>,
                std::tuple < std::false_type >> {}...);
        },
        new_subnodes_full);

    constexpr auto new_subnodes_full_filtered =
        filter(new_subnodes_full, calc_flags);

    constexpr auto bt_free =
        replace(bt, current_derivative_node, available_t{});
    constexpr auto new_subnodes_locations =
        locate_new_vals(new_subnodes_full_filtered, it, bt_free);

    constexpr auto new_bt = update_buffer_types(
        new_subnodes_full_filtered, new_subnodes_locations, bt_free);

    constexpr auto same_nodes =
        std::is_same_v<NodeDerivative1, NodeDerivative2>;

    double const this_val_derivative = get_differential_operator_value(
        current_derivative_node, it, ia, bt, ba);

    std::array<double, size(new_subnodes_full_filtered)>
        multinomial_expansion_values;

    if constexpr (is_constant_class2_v<NodeDerivative1>) {
        static_assert(size(new_subnodes_full_filtered) == 1);
        write_result(std::get<0>(new_subnodes_full_filtered),
                     -this_val_derivative, std::get<0>(new_subnodes_locations),
                     new_bt, ba, it, ia);
    } else if constexpr (is_constant_class2_v<NodeDerivative2>) {
        static_assert(size(new_subnodes_full_filtered) == 1);
        write_result(std::get<0>(new_subnodes_full_filtered),
                     this_val_derivative, std::get<0>(new_subnodes_locations),
                     new_bt, ba, it, ia);
    } else if constexpr (same_nodes) {
        static_assert(size(new_subnodes_full_filtered) == 1);
        write_result(std::get<0>(new_subnodes_full_filtered), 0.,
                     std::get<0>(new_subnodes_locations), new_bt, ba, it, ia);
    } else {
        constexpr auto multinomial_sequences = MultinomialSequences<2, Power>();
        constexpr auto multinomial_sequences_filtered =
            filter(multinomial_sequences, calc_flags);

        constexpr auto sign =
            signatures<NodesValue, NodeDerivative1, NodeDerivative2, Power>();

        constexpr auto sign_filtered = filter(sign, calc_flags);
        fill_with_signatures(multinomial_expansion_values, sign_filtered,
                             this_val_derivative);

        multiply_multivariate(multinomial_sequences_filtered,
                              multinomial_expansion_values);

        write_results(new_subnodes_full_filtered, multinomial_expansion_values,
                      new_subnodes_locations, new_bt, ba, it, ia);
    }

    backpropagate_aux(next_derivative_nodes, ct, it, ia, new_bt, ba, ut, ua);
}

template <std::size_t... Seq, std::size_t... Coeffs>
constexpr auto
binomial_mult(std::index_sequence<Seq...> /* sequence */,
              std::index_sequence<Coeffs...> /* binomail_coeffs */)
    -> std::size_t {
    return (pow(Coeffs, Seq) * ...);
}

template <std::size_t N = 0, class PartitionIntegerSequences, class ArrayOut,
          class BinomialCoeffs>
auto fill_with_mult(PartitionIntegerSequences const sequences,
                    ArrayOut &arrayout, BinomialCoeffs binomail_coeffs,
                    double value_derivative_node, double value_input_1,
                    double value_input_2) {
    if constexpr (N < std::tuple_size_v<PartitionIntegerSequences>) {
        constexpr auto current_sequence = std::get<N>(sequences);
        arrayout[N] = value_derivative_node;

        constexpr auto coeff = MultinomialCoeff2(current_sequence) *
                               binomial_mult(current_sequence, binomail_coeffs);
        if constexpr (coeff != 1) {
            arrayout[N] *= coeff;
        }

        constexpr auto first_power = get<0>(current_sequence);
        if constexpr (first_power != 0) {
            arrayout[N] *= std::pow(value_input_2, first_power);
        }

        constexpr auto last_power =
            get<current_sequence.size() - 1>(current_sequence);
        if constexpr (last_power != 0) {
            arrayout[N] *= std::pow(value_input_1, last_power);
        }

        fill_with_mult<N + 1>(sequences, arrayout, binomail_coeffs,
                              value_derivative_node, value_input_1,
                              value_input_2);
    }
}

template <std::size_t N = 0, class PartitionIntegerSequences, class ArrayOut,
          class BinomialCoeffs>
auto fill_with_mult_same(PartitionIntegerSequences const sequences,
                         ArrayOut &arrayout, BinomialCoeffs binomail_coeffs,
                         double value_derivative_node, double value_input) {
    if constexpr (N < std::tuple_size_v<PartitionIntegerSequences>) {
        constexpr auto current_sequence = std::get<N>(sequences);
        arrayout[N] = value_derivative_node;

        constexpr auto coeff = MultinomialCoeff2(current_sequence) *
                               binomial_mult(current_sequence, binomail_coeffs);
        if constexpr (coeff != 1) {
            arrayout[N] *= coeff;
        }

        constexpr auto first_power = get<0>(current_sequence);
        if constexpr (first_power != 0) {
            arrayout[N] *= std::pow(value_input, first_power);
        }

        fill_with_mult_same<N + 1>(sequences, arrayout, binomail_coeffs,
                                   value_derivative_node, value_input);
    }
}

template <std::size_t Power, std::size_t Order, class NodeDerivative1,
          class NodeDerivative2, class FirstNodesDerivativeRest,
          class CurrentNodesDerivatives, class NextNodesDerivatives,
          class CalcTreeValue, class InterfaceTypes, class InterfaceArray,
          class BufferTypes, class BufferArray, class UnivariateType,
          class UnivariateArray>
void backpropagate_process(
    der2::p<Power, der2::d<Order, mul_t<NodeDerivative1, NodeDerivative2>>> nd,
    FirstNodesDerivativeRest ndr,
    CurrentNodesDerivatives current_derivative_node,
    NextNodesDerivatives next_derivative_nodes, CalcTreeValue const &ct,
    InterfaceTypes it, InterfaceArray &ia, BufferTypes bt, BufferArray &ba,
    UnivariateType ut, UnivariateArray &ua) {

    constexpr auto current_node_loc =
        locate_val(current_derivative_node, it, bt);

    using NodesValue = CalcTreeValue::ValuesTupleInverse;
    constexpr auto new_subnodes = expand_tree_single(NodesValue{}, nd);
    constexpr auto new_subnodes_full =
        multiply_ordered_tuple(NodesValue{}, new_subnodes, ndr);

    constexpr auto calc_flags = std::apply(
        [next_derivative_nodes, it](auto... type) {
            return std::tuple_cat(
                std::conditional_t < contains(next_derivative_nodes, type) ||
                    contains(it, type),
                std::tuple<std::true_type>,
                std::tuple < std::false_type >> {}...);
        },
        new_subnodes_full);

    constexpr auto new_subnodes_full_filtered =
        filter(new_subnodes_full, calc_flags);

    constexpr auto bt_free =
        replace(bt, current_derivative_node, available_t{});
    constexpr auto new_subnodes_locations =
        locate_new_vals(new_subnodes_full_filtered, it, bt_free);

    constexpr auto new_bt = update_buffer_types(
        new_subnodes_full_filtered, new_subnodes_locations, bt_free);

    constexpr auto same_nodes =
        std::is_same_v<NodeDerivative1, NodeDerivative2>;

    double const this_val_derivative = get_differential_operator_value(
        current_derivative_node, it, ia, bt, ba);

    std::array<double, size(new_subnodes_full_filtered)>
        multinomial_expansion_values;

    if constexpr (is_constant_class2_v<NodeDerivative2>) {
        static_assert(size(new_subnodes_full_filtered) == 1);
        double const val2 = std::pow(ct.val(NodeDerivative2{}), Power);
        write_result(std::get<0>(new_subnodes_full_filtered),
                     this_val_derivative * val2,
                     std::get<0>(new_subnodes_locations), new_bt, ba, it, ia);
    } else if constexpr (is_constant_class2_v<NodeDerivative1>) {
        static_assert(size(new_subnodes_full_filtered) == 1);
        double const val1 = std::pow(ct.val(NodeDerivative1{}), Power);
        write_result(std::get<0>(new_subnodes_full_filtered),
                     this_val_derivative * val1,
                     std::get<0>(new_subnodes_locations), new_bt, ba, it, ia);
    } else if constexpr (same_nodes) {
        constexpr auto binomial_coeffs = BinomialCoefficientsSquare<Order>();

        constexpr auto multinomial_sequences =
            MultinomialSequences<binomial_coeffs.size(), Power>();
        constexpr auto multinomial_sequences_filtered =
            filter(multinomial_sequences, calc_flags);

        double const val = ct.val(NodeDerivative1{});

        fill_with_mult_same(multinomial_sequences_filtered,
                            multinomial_expansion_values, binomial_coeffs,
                            this_val_derivative, val);

        write_results(new_subnodes_full_filtered, multinomial_expansion_values,
                      new_subnodes_locations, new_bt, ba, it, ia);
    } else {
        constexpr auto multinomial_sequences =
            MultinomialSequences<Order + 1, Power>();
        constexpr auto multinomial_sequences_filtered =
            filter(multinomial_sequences, calc_flags);

        double const val1 = ct.val(NodeDerivative1{});
        double const val2 = ct.val(NodeDerivative2{});

        constexpr auto inverted =
            static_cast<bool>(get_idx_first2<NodeDerivative1>(NodesValue{}) >=
                              get_idx_first2<NodeDerivative2>(NodesValue{}));

        constexpr auto binomial_coeffs = BinomialCoefficients<Order>();

        fill_with_mult(multinomial_sequences_filtered,
                       multinomial_expansion_values, binomial_coeffs,
                       this_val_derivative, inverted ? val2 : val1,
                       inverted ? val1 : val2);

        write_results(new_subnodes_full_filtered, multinomial_expansion_values,
                      new_subnodes_locations, new_bt, ba, it, ia);
    }

    backpropagate_aux(next_derivative_nodes, ct, it, ia, new_bt, ba, ut, ua);
}

template <class NodeDerivatives, class CalcTreeValue, class InterfaceTypes,
          class InterfaceArray, class BufferTypes, class BufferArray,
          class UnivariateType, class UnivariateArray>
void backpropagate_aux(NodeDerivatives nd, CalcTreeValue const &ct,
                       InterfaceTypes it, InterfaceArray &ia, BufferTypes bt,
                       BufferArray &ba, UnivariateType ut,
                       UnivariateArray &ua) {

    if constexpr (std::tuple_size_v<NodeDerivatives>) {
        constexpr auto current_derivative_node = head(nd);
        constexpr auto next_derivative_nodes = tail(nd);

        constexpr auto current_derivative_subnode =
            head(current_derivative_node);
        constexpr auto current_derivative_subnode_rest =
            tail(current_derivative_node);

        backpropagate_process(current_derivative_subnode,
                              current_derivative_subnode_rest,
                              current_derivative_node, next_derivative_nodes,
                              ct, it, ia, bt, ba, ut, ua);
    }
}

} // namespace adhoc3::detail

#endif // ADHOC3_TAPE_BACKPROPAGATOR_HPP
