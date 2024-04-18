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
void backpropagate_aux(NodeDerivatives nd, CalcTreeValue ct, InterfaceTypes it,
                       InterfaceArray &ia, BufferTypes bt, BufferArray &ba,
                       UnivariateType ut, UnivariateArray &ua);

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
        if constexpr (is_on_interface) {
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
    NextNodesDerivatives next_derivative_nodes, CalcTreeValue ct,
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

    using NodesValue = decltype(ct)::ValuesTupleInverse;

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
        arrayout[N] *= MultinomialCoeff2(current_sequence);
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
    NextNodesDerivatives next_derivative_nodes, CalcTreeValue ct,
    InterfaceTypes it, InterfaceArray &ia, BufferTypes bt, BufferArray &ba,
    UnivariateType ut, UnivariateArray &ua) {

    constexpr auto current_node_loc =
        locate_val(current_derivative_node, it, bt);

    using NodesValue = decltype(ct)::ValuesTupleInverse;
    constexpr auto new_subnodes = expand_tree_single(NodesValue{}, nd);
    // std::cout << type_name2<decltype(new_subnodes)>() << std::endl;
    constexpr auto new_subnodes_full =
        multiply_ordered_tuple(NodesValue{}, new_subnodes, ndr);
    // std::cout << type_name2<decltype(new_subnodes)>() << std::endl;

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

    // static_assert(!same_nodes);

    if constexpr (is_constant_class2_v<NodeDerivative1>) {
        write_result(std::get<0>(new_subnodes_full_filtered),
                     this_val_derivative, std::get<0>(new_subnodes_locations),
                     new_bt, ba, it, ia);
    } else if constexpr (is_constant_class2_v<NodeDerivative2>) {
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

template <class NodeDerivatives, class CalcTreeValue, class InterfaceTypes,
          class InterfaceArray, class BufferTypes, class BufferArray,
          class UnivariateType, class UnivariateArray>
void backpropagate_aux(NodeDerivatives nd, CalcTreeValue ct, InterfaceTypes it,
                       InterfaceArray &ia, BufferTypes bt, BufferArray &ba,
                       UnivariateType ut, UnivariateArray &ua) {

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
