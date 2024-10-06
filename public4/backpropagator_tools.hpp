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

#ifndef ADHOC4_BACKPROPAGATOR_TOOLS_HPP
#define ADHOC4_BACKPROPAGATOR_TOOLS_HPP

#include "combinatorics/multinomial_coefficient_index_sequence.hpp"
#include "differential_operator.hpp"
#include "sort.hpp"
#include "utils/tuple.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <tuple>
#include <type_traits>

namespace adhoc4::detail {

template <std::size_t Order, class DiffOpNode, class Node>
constexpr auto first_type_is_aux(der::d<Order, DiffOpNode> /* in */,
                                 Node /* node */) -> bool {
    return std::is_same_v<Node, DiffOpNode>;
}

template <class DiffOp, class Node>
constexpr auto first_type_is(DiffOp diffop, Node node) -> bool {
    return first_type_is_aux(head(diffop), node);
}

template <class... Ids, std::size_t... Orders, std::size_t... Powers>
constexpr auto order(std::tuple<der::d<Orders, Ids>...> /* id1 */) {
    return (Orders + ...);
}

template <class... Ops>
constexpr auto max_orders(std::tuple<Ops...> /* id1 */) {
    constexpr std::array<std::size_t, sizeof...(Ops)> temp{order(Ops{})...};
    auto result = std::max_element(temp.begin(), temp.end());
    return *result;
}

template <std::size_t Offset, std::size_t PosFrom = 0, std::size_t PosTo = 0,
          class From, class To, class FilteredIS>
auto constexpr copy_filtered_aux(From const &from_array, To &to_array,
                                 FilteredIS fis) {
    if constexpr (PosFrom < size(fis)) {
        if constexpr (std::get<PosFrom>(fis)) {
            to_array[PosTo] = from_array[PosFrom + Offset];
            copy_filtered_aux<Offset, PosFrom + 1, PosTo + 1>(from_array,
                                                              to_array, fis);
        } else {
            copy_filtered_aux<Offset, PosFrom + 1, PosTo>(from_array, to_array,
                                                          fis);
        }
    }
}

template <class FromD, std::size_t FromS, class To, class Filter>
auto constexpr copy_filtered(std::array<FromD, FromS> const &from_array,
                             To &to_array, Filter filt) {
    constexpr std::size_t offset = FromS - size(filt);
    copy_filtered_aux<offset>(from_array, to_array, filt);
}

template <std::size_t Offset, std::size_t PosFrom = 0, std::size_t PosTo = 0,
          class From, class To, class FilteredIS>
void inline copy_filtered_inverted_aux(From const &from_array, To &to_array,
                                       FilteredIS fis) {
    if constexpr (PosFrom < size(fis)) {
        if constexpr (std::get<PosFrom>(fis)) {
            to_array[PosTo] = from_array[size(fis) - 1 - PosFrom + Offset];
            copy_filtered_inverted_aux<Offset, PosFrom + 1, PosTo + 1>(
                from_array, to_array, fis);
        } else {
            copy_filtered_inverted_aux<Offset, PosFrom + 1, PosTo>(
                from_array, to_array, fis);
        }
    }
}

template <class FromD, std::size_t FromS, class To, class Filter>
void inline copy_filtered_inverted(std::array<FromD, FromS> const &from_array,
                                   To &to_array, Filter filt) {
    static_assert(FromS >= size(filt));
    constexpr std::size_t offset = FromS - size(filt);
    copy_filtered_inverted_aux<offset>(from_array, to_array, filt);
}

class available_t {};
class on_buffer_t {};
class on_interface_t {};

template <class DerivativeNodeLoc, class DerivativeNode, class InterfaceTypes,
          class InterfaceArray, class BufferTypes, class BufferArray>
auto get_differential_operator_value(DerivativeNodeLoc current_node_der_loc,
                                     DerivativeNode current_node_der,
                                     InterfaceTypes it,
                                     InterfaceArray const &ia, BufferTypes bt,
                                     BufferArray const &ba) -> double {
    if constexpr (std::is_same_v<decltype(current_node_der_loc),
                                 on_interface_t>) {
        constexpr auto idx = get_first_type_idx(it, current_node_der);
        return ia[idx];
    } else {
        static_assert(
            std::is_same_v<decltype(current_node_der_loc), on_buffer_t>);
        constexpr auto idx = get_first_type_idx(bt, current_node_der);
        return ba[idx];
    }
}

template <class DerivativeNodeLoc, class DerivativeNode, class BufferTypes>
auto constexpr free_on_buffer(DerivativeNodeLoc current_node_der_loc,
                              DerivativeNode current_node_der, BufferTypes bt) {
    if constexpr (std::is_same_v<decltype(current_node_der_loc), on_buffer_t>) {
        constexpr auto newtuple =
            replace_first(bt, current_node_der, available_t{});
        return newtuple;
    } else {
        static_assert(
            std::is_same_v<decltype(current_node_der_loc), on_interface_t>);
        return bt;
    }
}

class on_buffer_add_t {};
class on_buffer_new_t {};
class on_interface_add_t {};

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
            is_on_interface, on_interface_add_t,
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
                replace_first(bt, available_t{}, current_derivative_node);
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

        if constexpr (std::is_same_v<const on_interface_add_t,
                                     decltype(current_indicator)>) {
            constexpr auto idx =
                get_first_type_idx(interface_types, current_derivative_node);

            // we always add on the interface
            interface_array[idx] += results_array[N];
        } else {
            constexpr auto idx = get_first_type_idx(buffer_types_updated,
                                                    current_derivative_node);

            if constexpr (std::is_same_v<const on_buffer_new_t,
                                         decltype(current_indicator)>) {
                // we place a new value so we override
                buffer_array[idx] = results_array[N];
            } else {
                static_assert(std::is_same_v<const on_buffer_add_t,
                                             decltype(current_indicator)>);
                // there is already a value so we add
                buffer_array[idx] += results_array[N];
            }
        }

        write_results<N + 1>(derivatives_nodes, results_array,
                             location_indicators, buffer_types_updated,
                             buffer_array, interface_types, interface_array);
    }
}

template <class Tuple, class Nodes>
constexpr auto is_sorted(Tuple in, Nodes nodes) {
    if constexpr (size(in) <= 1) {
        return true;
    } else {
        if constexpr (less_than(std::get<1>(in), std::get<0>(in), nodes)) {
            return is_sorted(tail(in), nodes);
        } else {
            return false;
        }
    }
}

template <class Tuple1, class Tuple2, class Out, class Nodes>
constexpr auto merge_sorted_aux(Tuple1 in1, Tuple2 in2, Out out, Nodes nodes) {
    if constexpr (size(in2) == 0) {
        return std::tuple_cat(out, in1);
    } else if constexpr (size(in1) == 0) {
        return std::tuple_cat(out, in2);
    } else {
        constexpr auto first1 = head(in1);
        constexpr auto first2 = head(in2);
        if constexpr (first1 == first2) {
            return merge_sorted_aux(tail(in1), tail(in2),
                                    std::tuple_cat(out, first1), nodes);
        } else if constexpr (less_than(first1, first2, nodes)) {
            return merge_sorted_aux(in1, tail(in2), std::tuple_cat(out, first2),
                                    nodes);
        } else {
            return merge_sorted_aux(tail(in1), in2, std::tuple_cat(out, first1),
                                    nodes);
        }
    }
}

template <class Tuple1, class Tuple2, class Nodes>
constexpr auto merge_sorted(Tuple1 in1, Tuple2 in2, Nodes nodes) {
    static_assert(is_sorted(in1, nodes));
    static_assert(is_sorted(in2, nodes));
    return merge_sorted_aux(in1, in2, std::tuple<>{}, nodes);
}

template <std::size_t I, class T> using Type = T;

template <class T, std::size_t... Is>
constexpr auto make_tuple_same_aux(std::index_sequence<Is...>) {
    return std::tuple<Type<Is, T>...>{};
}

template <std::size_t N, class T> constexpr auto make_tuple_same(T /* t */) {
    return make_tuple_same_aux<T>(std::make_index_sequence<N>{});
}

template <class Id1, class Id2, class Nodes>
constexpr auto sort_pair(std::tuple<Id1, Id2> id, Nodes nodes) {
    constexpr auto id1_less_than_id2 = static_cast<bool>(
        get_first_type_idx(nodes, Id1{}) >= get_first_type_idx(nodes, Id2{}));
    if constexpr (id1_less_than_id2) {
        return std::tuple<Id2, Id1>{};
    } else {
        return id;
    }
}

template <std::size_t Loc, std::size_t Size>
auto powers_fill(std::array<double, Size> &vals, double val) {
    if constexpr (Loc < Size) {
        vals[Loc] = vals[Loc - 1] * val;
        powers_fill<Loc + 1>(vals, val);
    }
}

template <std::size_t Order> auto powers(double a) {
    static_assert(Order > 0);
    std::array<double, Order> vals;
    vals[0] = a;
    powers_fill<1>(vals, a);
    return vals;
}

template <std::size_t N = 0, class PartitionIntegerSequences, class ArrayOut>
auto calc_add(PartitionIntegerSequences const sequences, ArrayOut &arrayout) {
    if constexpr (N < std::tuple_size_v<PartitionIntegerSequences>) {
        constexpr auto current_sequence = std::get<N>(sequences);
        if constexpr (MultinomialCoeff2(current_sequence) != 1) {
            arrayout[N] *=
                static_cast<double>(MultinomialCoeff2(current_sequence));
        }
        calc_add<N + 1>(sequences, arrayout);
    }
}

template <class DiffOp, class Result, class Nodes>
constexpr auto
multiply_differential_operator_aux(DiffOp single, std::tuple<> /* mult */,
                                   Result res, Nodes /* nodes */) {
    return std::tuple_cat(res, std::make_tuple(single));
}

template <class Id1, class Id2, class Nodes>
constexpr auto less_than_single(Id1 /* in1 */, Id2 /* in2 */, Nodes nodes) {
    constexpr auto idx1 = get_first_type_idx(nodes, Id1{});
    constexpr auto idx2 = get_first_type_idx(nodes, Id2{});
    return (idx1 > idx2);
}

template <class Single, class Many, class Result, class Nodes>
constexpr auto multiply_differential_operator_aux(Single single, Many mult,
                                                  Result res, Nodes nodes) {
    constexpr std::size_t order = get_power(single);
    constexpr auto id = get_id(single);

    constexpr auto mult_first = head(mult);
    constexpr auto id1 = get_id(mult_first);
    constexpr std::size_t order1 = get_power(mult_first);
    constexpr auto mult_tail = tail(mult);

    if constexpr (id == id1) {
        return std::tuple_cat(
            res, std::tuple<const der::d<order + order1, decltype(id1)>>{},
            mult_tail);
    } else if constexpr (less_than_single(id, id1, nodes)) {
        return multiply_differential_operator_aux(
            single, mult_tail, std::tuple_cat(res, std::make_tuple(mult_first)),
            nodes);
    } else {
        return std::tuple_cat(res, std::tuple<const Single>{}, mult);
    }
}

template <class Single, class Many, class Nodes>
constexpr auto multiply_differential_operator(Single const single,
                                              Many const mult, Nodes nodes) {
    return multiply_differential_operator_aux(single, mult, std::tuple<>{},
                                              nodes);
}

template <class... DerOps1, class DerOps2, class Nodes>
constexpr auto multiply_ordered_tuple(std::tuple<DerOps1...> const in1,
                                      DerOps2 const in2, Nodes nodes) {

    if constexpr (size(in2)) {
        constexpr auto next_derivatives = std::apply(
            [in2, nodes](const auto... singlein1) {
                return std::make_tuple(
                    multiply_differential_operator(singlein1, in2, nodes)...);
            },
            in1);
        return next_derivatives;
    } else {
        return in1;
    }
}

} // namespace adhoc4::detail

#endif // ADHOC4_BACKPROPAGATOR_TOOLS_HPP