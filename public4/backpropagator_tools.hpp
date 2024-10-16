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

#include "adhoc.hpp"
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

template <class Tuple, class T, std::size_t N = 0>
constexpr auto find4() -> std::size_t {
    if constexpr (N == std::tuple_size_v<Tuple>) {
        return N;
    } else if constexpr (std::is_same_v<const std::tuple_element_t<N, Tuple>,
                                        const T>) {
        return N;
    } else {
        return find4<Tuple, T, N + 1>();
    }
}

template <class I> constexpr auto is_derivative_input(I /* in */) -> bool {
    return false;
}

template <StringLiteral... literals, std::size_t... Orders>
constexpr auto
is_derivative_input(std::tuple<der::d<Orders, double_t<literals>>...> /* in */)
    -> bool {
    return true;
}

template <StringLiteral... literals, std::size_t... Orders>
constexpr auto is_derivative_input(
    std::tuple<const der::d<Orders, double_t<literals>>...> /* in */) -> bool {
    return true;
}

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
auto get_differential_operator_value(
    DerivativeNodeLoc /* current_node_der_loc */,
    DerivativeNode /* current_node_der */, InterfaceTypes /* it */,
    InterfaceArray const &ia, BufferTypes /* bt */, BufferArray const &ba)
    -> double {
    if constexpr (std::is_same_v<DerivativeNodeLoc, on_interface_t>) {
        constexpr auto idx = find4<InterfaceTypes, DerivativeNode>();
        return ia[idx];
    } else {
        static_assert(std::is_same_v<DerivativeNodeLoc, on_buffer_t>);
        constexpr auto idx = find4<BufferTypes, DerivativeNode>();
        return ba[idx];
    }
}

template <class Tuple, class New, std::size_t... ISFirst, std::size_t... ISLast>
auto constexpr replace_first3_aux(
    std::index_sequence<ISFirst...> /* is_first */,
    std::index_sequence<ISLast...> /* is_last */) {
    constexpr auto offset = sizeof...(ISFirst) + 1;
    return std::tuple<std::tuple_element_t<ISFirst, Tuple>..., New,
                      std::tuple_element_t<ISLast + offset, Tuple>...>{};
}

template <class Tuple, class Old, class New>
auto constexpr replace_first3(Tuple /* tuple */, Old /* old_value */,
                              New /* new_value */) {
    constexpr auto loc = find4<Tuple, Old>();
    constexpr auto tuple_size = std::tuple_size_v<Tuple>;
    static_assert(loc < tuple_size);

    return replace_first3_aux<Tuple, New>(
        std::make_index_sequence<loc>{},
        std::make_index_sequence<tuple_size - loc - 1>());
}

template <class DerivativeNodeLoc, class DerivativeNode, class BufferTypes>
auto constexpr free_on_buffer(DerivativeNodeLoc current_node_der_loc,
                              DerivativeNode current_node_der, BufferTypes bt) {
    if constexpr (std::is_same_v<decltype(current_node_der_loc), on_buffer_t>) {
        constexpr auto newtuple =
            replace_first3(bt, current_node_der, available_t{});
        return newtuple;
    } else {
        static_assert(
            std::is_same_v<decltype(current_node_der_loc), on_interface_t>);
        return bt;
    }
}

template <class Tuple, std::size_t... ISFirst, std::size_t... ISLast>
auto constexpr remove_first_aux(std::index_sequence<ISFirst...> /* is_first */,
                                std::index_sequence<ISLast...> /* is_last */) {
    constexpr auto offset = sizeof...(ISFirst) + 1;
    return std::tuple<std::tuple_element_t<ISFirst, Tuple>...,
                      std::tuple_element_t<ISLast + offset, Tuple>...>{};
}

template <class Tuple, class Old>
auto constexpr remove_first(Tuple /* tuple */, Old /* old_value */) {
    constexpr auto loc = find4<Tuple, Old>();
    constexpr auto tuple_size = std::tuple_size_v<Tuple>;
    static_assert(loc < tuple_size);

    return remove_first_aux<Tuple>(
        std::make_index_sequence<loc>{},
        std::make_index_sequence<tuple_size - loc - 1>());
}

template <class DerivativeNodeLoc, class DerivativeNode, class BufferTypes>
auto constexpr free_on_buffer_size(DerivativeNodeLoc current_node_der_loc,
                                   DerivativeNode current_node_der,
                                   BufferTypes bt) {
    if constexpr (std::is_same_v<decltype(current_node_der_loc), on_buffer_t>) {
        constexpr auto newtuple = remove_first(bt, current_node_der);
        return newtuple;
    } else {
        static_assert(
            std::is_same_v<decltype(current_node_der_loc), on_interface_t>);
        return bt;
    }
}

template <std::size_t Loc> class on_buffer_add_t {};
class on_buffer_new_t {};
template <std::size_t Loc> class on_interface_add_t {};

template <class T, template <std::size_t> class U>
inline constexpr bool is_instance_of_v = std::false_type{};

template <template <std::size_t> class U, std::size_t V>
inline constexpr bool is_instance_of_v<U<V>, U> = std::true_type{};

template <std::size_t Loc> constexpr std::size_t get_loc(on_buffer_add_t<Loc>) {
    return Loc;
}

template <std::size_t Loc>
constexpr std::size_t get_loc(on_interface_add_t<Loc>) {
    return Loc;
}

template <class Type, class InterfaceTypes, class BufferTypes>
constexpr auto locate_new_vals_aux() {

    if constexpr (is_derivative_input(Type{})) {
        constexpr auto it_loc = find4<InterfaceTypes, Type>();
        constexpr auto it_size = std::tuple_size_v<InterfaceTypes>;
        static_assert(it_loc < it_size);
        return on_interface_add_t<it_loc>{};
    } else {
        constexpr auto bt_loc = find4<BufferTypes, Type>();
        constexpr auto bt_size = std::tuple_size_v<BufferTypes>;
        if constexpr (bt_loc < bt_size) {
            return on_buffer_add_t<bt_loc>{};
        } else {
            constexpr auto it_loc = find4<InterfaceTypes, Type>();
            constexpr auto it_size = std::tuple_size_v<InterfaceTypes>;
            if constexpr (it_loc < it_size) {
                return on_interface_add_t<it_loc>{};
            } else {
                return on_buffer_new_t{};
            }
        }
    }
}

template <class... Types, class InterfaceTypes, class BufferTypes>
constexpr auto locate_new_vals(std::tuple<Types...> /* derivative_nodes */,
                               InterfaceTypes /* it */, BufferTypes /* bt */) {
    return std::make_tuple(
        locate_new_vals_aux<Types, InterfaceTypes, BufferTypes>()...);
}

template <std::size_t N = 0, std::size_t FindOffset = 0, class TypesToPlace,
          class LocationIndicators, class BufferTypes>
constexpr auto update_buffer_types(TypesToPlace derivative_nodes,
                                   LocationIndicators location_indicators,
                                   BufferTypes bt) {
    if constexpr (N < std::tuple_size_v<LocationIndicators>) {
        if constexpr (std::is_same_v<
                          const on_buffer_new_t,
                          const std::tuple_element_t<N, LocationIndicators>>) {
            constexpr auto loc = find4<BufferTypes, available_t, FindOffset>();
            constexpr auto tuple_size = std::tuple_size_v<BufferTypes>;
            static_assert(loc < tuple_size);

            constexpr auto new_buffer_type =
                replace_first3_aux<BufferTypes,
                                   std::tuple_element_t<N, TypesToPlace>>(
                    std::make_index_sequence<loc>{},
                    std::make_index_sequence<tuple_size - loc - 1>());

            return update_buffer_types<N + 1, loc + 1>(
                derivative_nodes, location_indicators, new_buffer_type);
        } else {
            return update_buffer_types<N + 1, FindOffset>(
                derivative_nodes, location_indicators, bt);
        }
    } else {
        return bt;
    }
}

template <class... Types, class... Location, class BufferTypes>
constexpr auto
update_buffer_types_size(std::tuple<Types...> /* derivative_nodes */,
                         std::tuple<Location...> /* location_indicators */,
                         BufferTypes bt) {
    return std::tuple_cat(
        bt, std::conditional_t<
                std::is_same_v<const Location, const on_buffer_new_t>,
                std::tuple<Types>, std::tuple<>>{}...);
}

template <std::size_t N = 0, class ResultsTypes, class ResultsArray,
          class LocationIndicators, class BufferTypes,
          std::size_t BufferArraySize, class InterfaceTypes,
          std::size_t InterfaceArraySize>
auto write_results(ResultsTypes derivatives_nodes,
                   ResultsArray const &results_array,
                   LocationIndicators location_indicators,
                   BufferTypes buffer_types_updated,
                   std::array<double, BufferArraySize> &buffer_array,
                   InterfaceTypes interface_types,
                   std::array<double, InterfaceArraySize> &interface_array) {
    if constexpr (N < std::tuple_size_v<LocationIndicators>) {
        if constexpr (is_instance_of_v<std::decay_t<std::tuple_element_t<
                                           N, LocationIndicators>>,
                                       on_interface_add_t>) {
            constexpr auto idx =
                get_loc(std::tuple_element_t<N, LocationIndicators>{});

            static_assert(idx < InterfaceArraySize);

            // we always add on the interface
            interface_array[idx] += results_array[N];
        } else {
            if constexpr (std::is_same_v<const on_buffer_new_t,
                                         const std::tuple_element_t<
                                             N, LocationIndicators>>) {
                constexpr auto idx =
                    find4<BufferTypes, std::tuple_element_t<N, ResultsTypes>>();
                static_assert(idx < BufferArraySize);
                buffer_array[idx] = results_array[N];
            } else {
                static_assert(
                    is_instance_of_v<std::decay_t<std::tuple_element_t<
                                         N, LocationIndicators>>,
                                     on_buffer_add_t>);
                constexpr auto idx =
                    get_loc(std::tuple_element_t<N, LocationIndicators>{});
                static_assert(idx < BufferArraySize);
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
            res,
            std::tuple<der::d<order + order1, decltype(get_id(head(mult)))>>{},
            mult_tail);
    } else if constexpr (less_than_single(id, id1, nodes)) {
        return multiply_differential_operator_aux(
            single, mult_tail, std::tuple_cat(res, std::make_tuple(head(mult))),
            nodes);
    } else {
        return std::tuple_cat(res, std::tuple<Single>{}, mult);
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

template <class Op1, class Op2, class Nodes>
constexpr auto less_than_check_empty2(Op1 /* in1 */, Op2 /* in2 */,
                                      Nodes /* nodes */);

template <class Id1, std::size_t Order1, class Id2, std::size_t Order2,
          class Id1Remaining, class Id2Remaining, class Nodes>
constexpr auto
less_than_single2(der::d<Order1, Id1> /* in1 */, der::d<Order2, Id2> /* in2 */,
                  Id1Remaining id1rem, Id2Remaining id2rem, Nodes nodes) {
    constexpr auto idx1 = get_first_type_idx(nodes, Id1{});
    constexpr auto idx2 = get_first_type_idx(nodes, Id2{});
    if constexpr (idx1 < idx2) {
        return false;
    } else if constexpr (idx1 > idx2) {
        return true;
    } else if constexpr (Order1 < Order2) {
        return true;
    } else if constexpr (Order1 > Order2) {
        return false;
    } else {
        static_assert(std::is_same_v<der::d<Order1, const Id1>,
                                     der::d<Order2, const Id2>>);
        return less_than_check_empty2(id1rem, id2rem, nodes);
    }
}

template <class In1, class In2, class Nodes>
constexpr auto less_than_check_first2(In1 in1, In2 in2, Nodes nodes) {
    return less_than_single2(head(in1), head(in2), tail(in1), tail(in2), nodes);
}

template <class Op1, class Op2, class Nodes>
constexpr auto less_than_check_empty2(Op1 in1, Op2 in2, Nodes nodes) {
    // if constexpr (size(in2) == 0) {
    if constexpr (std::is_same_v<const Op2, const std::tuple<>>) {
        return false;
        // } else if constexpr (size(in1) == 0) {
    } else if constexpr (std::is_same_v<const Op1, const std::tuple<>>) {
        return true;
    } else {
        // we can now assume that both ids are different and non empty
        return less_than_check_first2(in1, in2, nodes);
    }
}

template <class Op1, class Op2, class Nodes>
constexpr auto less_than2(Op1 in1, Op2 in2, Nodes nodes) {
    if constexpr (std::is_same_v<Op1, Op2>) {
        return false;
    } else {
        return less_than_check_empty2(in1, in2, nodes);
    }
}

template <class Tuple, class Nodes>
constexpr auto is_sorted2(Tuple in, Nodes nodes) {
    if constexpr (size(in) <= 1) {
        return true;
    } else {
        if constexpr (less_than2(std::get<0>(in), std::get<1>(in), nodes)) {
            return false;
        } else {
            return is_sorted2(tail(in), nodes);
        }
    }
}

template <class I1, class I2>
constexpr auto are_same(I1 /* lhs */, I2 /* rhs */) -> bool {
    return std::is_same_v<I1, I2>;
}

template <class Tuple1, class Tuple2, class Out, class Nodes>
constexpr auto merge_sorted_aux2(Tuple1 in1, Tuple2 in2, Out out, Nodes nodes) {
    if constexpr (size(in2) == 0) {
        return std::tuple_cat(out, in1);
    } else if constexpr (size(in1) == 0) {
        return std::tuple_cat(out, in2);
    } else {
        constexpr auto first1 = head(in1);
        constexpr auto first2 = head(in2);
        if constexpr (are_same(first1, first2)) {
            return merge_sorted_aux2(
                tail(in1), tail(in2),
                std::tuple_cat(out, std::make_tuple(first1)), nodes);
        } else if constexpr (less_than2(first1, first2, nodes)) {
            return merge_sorted_aux2(
                in1, tail(in2), std::tuple_cat(out, std::make_tuple(first2)),
                nodes);
        } else {
            return merge_sorted_aux2(
                tail(in1), in2, std::tuple_cat(out, std::make_tuple(first1)),
                nodes);
        }
    }
}

template <class Tuple1, class Tuple2, class Nodes>
constexpr auto merge_sorted2(Tuple1 in1, Tuple2 in2, Nodes nodes) {
    static_assert(is_sorted2(in1, nodes));
    static_assert(is_sorted2(in2, nodes));
    return merge_sorted_aux2(in1, in2, std::tuple<>{}, nodes);
}

template <std::size_t N = 0, std::size_t OutSize, std::size_t MaxOrder,
          class IS>
auto calc_mul(std::array<double, OutSize> &arrayout,
              std::array<double, MaxOrder> const &powers_val1,
              std::array<double, MaxOrder> const &powers_val2, IS et) {
    static_assert(OutSize == std::tuple_size_v<IS>);
    if constexpr (N < OutSize) {
        constexpr auto current_sequence = std::get<N>(et);
        static_assert(current_sequence.size() == 3);
        if constexpr (get<1>(current_sequence)) {
            arrayout[N] *= powers_val2[get<1>(current_sequence) - 1];
        }

        if constexpr (get<2>(current_sequence)) {
            arrayout[N] *= powers_val1[get<2>(current_sequence) - 1];
        }

        if constexpr (MultinomialCoeff2(current_sequence) != 1) {
            arrayout[N] *=
                static_cast<double>(MultinomialCoeff2(current_sequence));
        }

        calc_mul<N + 1>(arrayout, powers_val1, powers_val2, et);
    }
}

template <class DerOp1, class DerOp2, class Result, class Nodes>
constexpr auto multiply_ordered_aux(DerOp1 const in1, DerOp2 const in2,
                                    Result res, Nodes nodes) {
    if constexpr (size(in2) == 0) {
        return std::tuple_cat(res, in1);
    } else if constexpr (size(in1) == 0) {
        return std::tuple_cat(res, in2);
    } else {
        constexpr auto in1_first = head(in1);
        constexpr auto in2_first = head(in2);
        constexpr auto id1 = get_id(in1_first);
        constexpr auto id2 = get_id(in2_first);
        constexpr auto idx1 = find4<Nodes, decltype(id1)>();
        constexpr auto idx2 = find4<Nodes, decltype(id2)>();

        if constexpr (idx1 < idx2) {
            constexpr std::size_t power = get_power(in1_first);
            return multiply_ordered_aux(
                tail(in1), in2, std::tuple_cat(res, d<power>(id1)), nodes);
        } else if constexpr (idx1 > idx2) {
            constexpr std::size_t power = get_power(in2_first);
            return multiply_ordered_aux(
                in1, tail(in2), std::tuple_cat(res, d<power>(id2)), nodes);
        } else {
            constexpr std::size_t power =
                get_power(in1_first) + get_power(in2_first);
            return multiply_ordered_aux(tail(in1), tail(in2),
                                        std::tuple_cat(res, d<power>(id1)),
                                        nodes);
        }
    }
}

template <class Tuple, class Nodes>
constexpr auto is_sorted_diff_op(Tuple in, Nodes nodes) {
    if constexpr (size(in) <= 1) {
        return true;
    } else {
        constexpr auto in_first = std::get<0>(in);
        constexpr auto in_second = std::get<1>(in);
        constexpr auto id1 = get_id(in_first);
        constexpr auto id2 = get_id(in_second);
        constexpr auto idx1 = find4<Nodes, decltype(id1)>();
        constexpr auto idx2 = find4<Nodes, decltype(id2)>();
        if constexpr (idx1 < idx2) {
            return true;
        } else {
            return is_sorted_diff_op(tail(in), nodes);
        }
    }
}

template <class DerOp1, class DerOp2, class Nodes>
constexpr auto multiply_ordered(DerOp1 const in1, DerOp2 const in2,
                                Nodes nodes) {
    static_assert(is_sorted_diff_op(in1, nodes));
    static_assert(is_sorted_diff_op(in2, nodes));
    return multiply_ordered_aux(in1, in2, std::tuple<>{}, nodes);
}

} // namespace adhoc4::detail

#endif // ADHOC4_BACKPROPAGATOR_TOOLS_HPP
