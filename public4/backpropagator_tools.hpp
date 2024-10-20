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
#include "utils/index_sequence.hpp"
#include "utils/tuple.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <tuple>
#include <type_traits>

namespace adhoc4::detail {

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
constexpr auto copy_filtered_aux(From const &from_array, To &to_array,
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
constexpr auto copy_filtered(std::array<FromD, FromS> const &from_array,
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

template <std::size_t From, std::size_t To, std::size_t PosOut = 0,
          std::size_t PosFilter = 0, class In, class Out, class FilteredIS>
void inline copy_filtered_inverted_aux2(In const &from_array, Out &to_array,
                                        FilteredIS fis) {
    constexpr std::size_t ToNext = To - 1;
    if constexpr (get<PosFilter>(fis)) {
        to_array[PosOut] = from_array[ToNext];
        if constexpr (ToNext != From) {
            copy_filtered_inverted_aux2<From, ToNext, PosOut + 1,
                                        PosFilter + 1>(from_array, to_array,
                                                       fis);
        }
    } else {
        if constexpr (ToNext != From) {
            copy_filtered_inverted_aux2<From, ToNext, PosOut, PosFilter + 1>(
                from_array, to_array, fis);
        }
    }
}

template <std::size_t From, std::size_t To, class D, std::size_t InSize,
          std::size_t OutSize, class Filter>
void inline copy_filtered_inverted2(std::array<D, InSize> const &from_array,
                                    std::array<D, OutSize> &to_array,
                                    Filter filt) {
    static_assert(From <= To);
    static_assert((To - From) == filt.size());
    copy_filtered_inverted_aux2<From, To>(from_array, to_array, filt);
}

class on_buffer_t {};
class on_interface_t {};
class on_buffer_not_new_t {};

template <class Where, std::size_t I, class InterfaceArray, class BufferArray>
auto get_differential_operator_value(
    std::pair<Where, std::integral_constant<std::size_t,
                                            I>> /* current_node_der_loc */,
    InterfaceArray const &ia, BufferArray const &ba) -> double {

    if constexpr (std::is_same_v<Where, on_interface_t>) {
        return ia[I];
    } else {
        static_assert(std::is_same_v<Where, on_buffer_t>);
        return ba[I];
    }
}

template <class Tuple, std::size_t... ISFirst, std::size_t... ISLast>
constexpr auto remove_aux(std::index_sequence<ISFirst...> /* is_first */,
                          std::index_sequence<ISLast...> /* is_last */) {
    constexpr auto offset = sizeof...(ISFirst) + 1;
    return std::tuple<std::tuple_element_t<ISFirst, Tuple>...,
                      std::tuple_element_t<ISLast + offset, Tuple>...>{};
}

template <class Tuple, std::size_t N> constexpr auto remove() {
    constexpr auto tuple_size = std::tuple_size_v<Tuple>;
    static_assert(N < tuple_size);
    return remove_aux<Tuple>(std::make_index_sequence<N>{},
                             std::make_index_sequence<tuple_size - N - 1>());
}

template <class Where, std::size_t I, class BufferTypes>
constexpr auto free_on_buffer(
    std::pair<Where, std::integral_constant<std::size_t, I>> /* loc */,
    BufferTypes bt) {
    if constexpr (std::is_same_v<Where, on_buffer_t>) {
        static_assert(get<I>(bt));
        return replace<I, false>(bt);
    } else {
        static_assert(std::is_same_v<Where, on_interface_t>);
        return bt;
    }
}

template <class Tuple, std::size_t... ISFirst, std::size_t... ISLast>
constexpr auto remove_first_aux(std::index_sequence<ISFirst...> /* is_first */,
                                std::index_sequence<ISLast...> /* is_last */) {
    constexpr auto offset = sizeof...(ISFirst) + 1;
    return std::tuple<std::tuple_element_t<ISFirst, Tuple>...,
                      std::tuple_element_t<ISLast + offset, Tuple>...>{};
}

template <class Tuple, class Old>
constexpr auto remove_first(Tuple /* tuple */, Old /* old_value */) {
    constexpr auto loc = find<Tuple, Old>();
    constexpr auto tuple_size = std::tuple_size_v<Tuple>;
    static_assert(loc < tuple_size);

    return remove_first_aux<Tuple>(
        std::make_index_sequence<loc>{},
        std::make_index_sequence<tuple_size - loc - 1>());
}

template <class Where, std::size_t I>
constexpr auto free_on_buffer_size(
    const std::pair<Where, std::integral_constant<std::size_t, I>> /* loc */)
    -> bool {
    if constexpr (std::is_same_v<Where, on_buffer_t>) {
        return true;
    } else {
        static_assert(std::is_same_v<Where, on_interface_t>);
        return false;
    }
}

template <class Tuple, class T, std::size_t N = 0>
constexpr auto find_first() -> std::size_t {
    if constexpr (N == std::tuple_size_v<Tuple>) {
        return N;
    } else if constexpr (std::is_same_v<const std::tuple_element_t<
                                            0, std::tuple_element_t<N, Tuple>>,
                                        const T>) {
        return N;
    } else {
        return find_first<Tuple, T, N + 1>();
    }
}

template <class Tuple, class T, std::size_t N = 0>
constexpr auto find_first_type_not() -> std::size_t {
    if constexpr (N == std::tuple_size_v<Tuple>) {
        return N;
    } else if constexpr (!detail::first_type_is(
                             std::tuple_element_t<
                                 0, std::tuple_element_t<N, Tuple>>{},
                             T{})) {
        return N;
    } else {
        return find_first_type_not<Tuple, T, N + 1>();
    }
}

template <std::size_t Last, std::size_t N = 0, std::size_t FindOffset = 0,
          class DerivativeNodesToPlace, class InterfaceTypes, class BufferTypes,
          class DerivNodeLoc, class DerivativeNodes, class DerivativeNodeNew>
constexpr auto locate_new_vals_update_buffer_types_aux(
    DerivativeNodesToPlace derivative_nodes, InterfaceTypes it, BufferTypes bt,
    DerivNodeLoc dnl, DerivativeNodes dn, DerivativeNodeNew dnn) {

    if constexpr (N < std::tuple_size_v<DerivativeNodesToPlace>) {
        using Type = std::tuple_element_t<N, DerivativeNodesToPlace>;
        if constexpr (is_derivative_input(Type{})) {
            constexpr auto it_loc = find<InterfaceTypes, Type>();
            constexpr auto it_size = std::tuple_size_v<InterfaceTypes>;
            static_assert(it_loc < it_size);
            return locate_new_vals_update_buffer_types_aux<Last, N + 1,
                                                           FindOffset>(
                derivative_nodes, it, bt,
                std::tuple_cat(
                    dnl, std::tuple<std::pair<
                             on_interface_t,
                             std::integral_constant<std::size_t, it_loc>>>{}),
                dn, dnn);
        } else {
            constexpr auto dn_loc = find_first<DerivativeNodes, Type, Last>();
            constexpr auto dn_size = std::tuple_size_v<DerivativeNodes>;
            if constexpr (dn_loc < dn_size) {
                return locate_new_vals_update_buffer_types_aux<Last, N + 1,
                                                               FindOffset>(
                    derivative_nodes, it, bt,
                    std::tuple_cat(
                        dnl,
                        std::tuple<
                            std::pair<on_buffer_not_new_t,
                                      std::integral_constant<
                                          std::size_t, std::get<dn_loc>(dn)
                                                           .second.second>>>{}),
                    dn, dnn);
            } else {
                constexpr auto dnn_loc = find_first<DerivativeNodeNew, Type>();
                constexpr auto dnn_size = std::tuple_size_v<DerivativeNodeNew>;

                if constexpr (dnn_loc < dnn_size) {
                    return locate_new_vals_update_buffer_types_aux<Last, N + 1,
                                                                   FindOffset>(
                        derivative_nodes, it, bt,
                        std::tuple_cat(
                            dnl,
                            std::tuple<std::pair<
                                on_buffer_not_new_t,
                                std::integral_constant<
                                    std::size_t,
                                    std::get<dnn_loc>(dnn).second.second>>>{}),
                        dn, dnn);
                } else {
                    constexpr auto it_loc = find<InterfaceTypes, Type>();
                    constexpr auto it_size = std::tuple_size_v<InterfaceTypes>;
                    if constexpr (it_loc < it_size) {
                        return locate_new_vals_update_buffer_types_aux<
                            Last, N + 1, FindOffset>(
                            derivative_nodes, it, bt,
                            std::tuple_cat(
                                dnl,
                                std::tuple<
                                    std::pair<on_interface_t,
                                              std::integral_constant<
                                                  std::size_t, it_loc>>>{}),
                            dn, dnn);
                    } else {
                        constexpr auto bt_loc = find<false, FindOffset>(bt);
                        constexpr auto tuple_size = bt.size();
                        static_assert(bt_loc < tuple_size);
                        constexpr auto bt_new = replace<bt_loc, true>(bt);

                        return locate_new_vals_update_buffer_types_aux<
                            Last, N + 1, bt_loc + 1>(
                            derivative_nodes, it, bt_new,
                            std::tuple_cat(
                                dnl,
                                std::tuple<std::pair<
                                    on_buffer_t, std::integral_constant<
                                                     std::size_t, bt_loc>>>{}),
                            dn, dnn);
                    }
                }
            }
        }

    } else {
        return std::make_tuple(bt, dnl);
    }
}

template <std::size_t Last, class DerivativeNodesToPlace, class InterfaceTypes,
          class BufferTypes, class DerivativeNodes, class DerivativeNodeNew>
constexpr auto
locate_new_vals_update_buffer_types(DerivativeNodesToPlace derivative_nodes,
                                    InterfaceTypes it, BufferTypes bt,
                                    DerivativeNodes dn, DerivativeNodeNew dnn) {
    return locate_new_vals_update_buffer_types_aux<Last>(
        derivative_nodes, it, bt, std::tuple<>{}, dn, dnn);
}

template <std::size_t Last, std::size_t N = 0, class DerivativeNodesToPlace,
          class InterfaceTypes, class DerivNodeLoc, class DerivativeNodes,
          class DerivativeNodeNew>
constexpr auto locate_new_vals_update_buffer_types_aux_size(
    DerivativeNodesToPlace derivative_nodes, InterfaceTypes it,
    DerivNodeLoc dnl, DerivativeNodes dn, DerivativeNodeNew dnn) {

    if constexpr (N < std::tuple_size_v<DerivativeNodesToPlace>) {
        using Type = std::tuple_element_t<N, DerivativeNodesToPlace>;
        if constexpr (is_derivative_input(Type{})) {
            constexpr auto it_loc = find<InterfaceTypes, Type>();
            constexpr auto it_size = std::tuple_size_v<InterfaceTypes>;
            static_assert(it_loc < it_size);
            return locate_new_vals_update_buffer_types_aux_size<Last, N + 1>(
                derivative_nodes, it,
                std::tuple_cat(
                    dnl, std::tuple<std::pair<
                             on_interface_t,
                             std::integral_constant<std::size_t, it_loc>>>{}),
                dn, dnn);
        } else {
            constexpr auto dn_loc = find_first<DerivativeNodes, Type, Last>();
            constexpr auto dn_size = std::tuple_size_v<DerivativeNodes>;
            if constexpr (dn_loc < dn_size) {
                return locate_new_vals_update_buffer_types_aux_size<Last,
                                                                    N + 1>(
                    derivative_nodes, it,
                    std::tuple_cat(
                        dnl,
                        std::tuple<
                            std::pair<on_buffer_not_new_t,
                                      std::integral_constant<
                                          std::size_t, std::get<dn_loc>(dn)
                                                           .second.second>>>{}),
                    dn, dnn);
            } else {
                constexpr auto dnn_loc = find_first<DerivativeNodeNew, Type>();
                constexpr auto dnn_size = std::tuple_size_v<DerivativeNodeNew>;

                if constexpr (dnn_loc < dnn_size) {
                    return locate_new_vals_update_buffer_types_aux_size<Last,
                                                                        N + 1>(
                        derivative_nodes, it,
                        std::tuple_cat(
                            dnl,
                            std::tuple<std::pair<
                                on_buffer_not_new_t,
                                std::integral_constant<
                                    std::size_t,
                                    std::get<dnn_loc>(dnn).second.second>>>{}),
                        dn, dnn);
                } else {
                    constexpr auto it_loc = find<InterfaceTypes, Type>();
                    constexpr auto it_size = std::tuple_size_v<InterfaceTypes>;
                    if constexpr (it_loc < it_size) {
                        return locate_new_vals_update_buffer_types_aux_size<
                            Last, N + 1>(
                            derivative_nodes, it,
                            std::tuple_cat(
                                dnl,
                                std::tuple<
                                    std::pair<on_interface_t,
                                              std::integral_constant<
                                                  std::size_t, it_loc>>>{}),
                            dn, dnn);
                    } else {
                        return locate_new_vals_update_buffer_types_aux_size<
                            Last, N + 1>(
                            derivative_nodes, it,
                            std::tuple_cat(
                                dnl,
                                std::tuple<std::pair<
                                    on_buffer_t,
                                    std::integral_constant<std::size_t, 0>>>{}),
                            dn, dnn);
                    }
                }
            }
        }
    } else {
        return dnl;
    }
}

template <std::size_t Last, class DerivativeNodesToPlace, class InterfaceTypes,
          class DerivativeNodes, class DerivativeNodeNew>
constexpr auto locate_new_vals_update_buffer_types_size(
    DerivativeNodesToPlace derivative_nodes, InterfaceTypes it,
    DerivativeNodes dn, DerivativeNodeNew dnn) {
    return locate_new_vals_update_buffer_types_aux_size<Last>(
        derivative_nodes, it, std::tuple<>{}, dn, dnn);
}

template <std::size_t N = 0, class ResultsArray, class LocationIndicators,
          std::size_t BufferArraySize, std::size_t InterfaceArraySize>
auto write_results(ResultsArray const &results_array,
                   LocationIndicators location_indicators,
                   std::array<double, BufferArraySize> &buffer_array,
                   std::array<double, InterfaceArraySize> &interface_array) {
    if constexpr (N < std::tuple_size_v<LocationIndicators>) {

        using Indicator = std::tuple_element_t<N, LocationIndicators>;
        using Where = std::tuple_element_t<0, Indicator>;
        using Loc = std::tuple_element_t<1, Indicator>;

        if constexpr (std::is_same_v<Where, on_interface_t>) {
            interface_array[Loc::value] += results_array[N];
        } else if constexpr (std::is_same_v<Where, on_buffer_not_new_t>) {
            buffer_array[Loc::value] += results_array[N];
        } else {
            static_assert(std::is_same_v<Where, on_buffer_t>);
            buffer_array[Loc::value] = results_array[N];
        }

        write_results<N + 1>(results_array, location_indicators, buffer_array,
                             interface_array);
    }
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
constexpr auto sort_pair(std::tuple<Id1, Id2> id, Nodes /* nodes */) {
    constexpr auto idx1 = find<Nodes, Id1>();
    constexpr auto idx2 = find<Nodes, Id2>();
    constexpr auto id1_less_than_id2 = static_cast<bool>(idx1 >= idx2);
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
constexpr auto less_than_single(Id1 /* in1 */, Id2 /* in2 */,
                                Nodes /* nodes */) {
    constexpr auto idx1 = find<Nodes, Id1>();
    constexpr auto idx2 = find<Nodes, Id2>();
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

template <std::size_t Pos1 = 0, std::size_t Pos2 = 0, class Op1, class Op2,
          class Nodes>
constexpr auto less_than3_aux(Op1 in1, Op2 in2, Nodes nodes) {
    if constexpr (Pos2 == std::tuple_size_v<Op2>) {
        return false;
    } else if constexpr (Pos1 == std::tuple_size_v<Op1>) {
        return true;
    } else {
        constexpr std::tuple_element_t<Pos1, Op1> diffop1;
        constexpr std::tuple_element_t<Pos2, Op2> diffop2;
        constexpr auto idx1 = find<Nodes, decltype(get_id(diffop1))>();
        constexpr auto idx2 = find<Nodes, decltype(get_id(diffop2))>();
        constexpr std::size_t Order1 = get_power(diffop1);
        constexpr std::size_t Order2 = get_power(diffop2);
        if constexpr (idx1 < idx2) {
            return false;
        } else if constexpr (idx1 > idx2) {
            return true;
        } else if constexpr (Order1 < Order2) {
            return true;
        } else if constexpr (Order1 > Order2) {
            return false;
        } else {
            static_assert(std::is_same_v<decltype(diffop1), decltype(diffop2)>);
            return less_than3_aux<Pos1 + 1, Pos2 + 1>(in1, in2, nodes);
        }
    }
}

template <class Op1, class Op2, class Nodes>
constexpr auto less_than2(Op1 in1, Op2 in2, Nodes nodes) {
    if constexpr (std::is_same_v<Op1, Op2>) {
        return false;
    } else {
        return less_than3_aux(in1, in2, nodes);
    }
}

template <std::size_t Pos = 0, class Tuple, class Nodes>
constexpr auto is_sorted(Tuple in, Nodes nodes) -> bool {
    if constexpr (std::tuple_size_v<Tuple> <= (Pos + 1)) {
        return true;
    } else {
        using Id1 = std::tuple_element_t<0, std::tuple_element_t<Pos, Tuple>>;
        using Id2 =
            std::tuple_element_t<0, std::tuple_element_t<Pos + 1, Tuple>>;

        if constexpr (less_than2(Id1{}, Id2{}, nodes)) {
            return false;
        } else {
            return is_sorted<Pos + 1>(in, nodes);
        }
    }
}

template <class I1, class I2>
constexpr auto are_same(I1 /* lhs */, I2 /* rhs */) -> bool {
    return std::is_same_v<I1, I2>;
}

template <class Tuple1, class Tuple2, std::size_t Pos1 = 0,
          std::size_t Pos2 = 0, class... Out, class Nodes>
constexpr auto merge_sorted_aux2(std::tuple<Out...> out, Nodes nodes) {
    if constexpr (std::tuple_size_v<Tuple1> == Pos1 &&
                  std::tuple_size_v<Tuple2> == Pos2) {
        return out;
    } else if constexpr (std::tuple_size_v<Tuple2> == Pos2) {
        return merge_sorted_aux2<Tuple1, Tuple2, Pos1 + 1, Pos2>(
            std::tuple<Out...,
                       std::pair<std::true_type,
                                 std::integral_constant<std::size_t, Pos1>>>{},
            nodes);
    } else if constexpr (std::tuple_size_v<Tuple1> == Pos1) {
        return merge_sorted_aux2<Tuple1, Tuple2, Pos1, Pos2 + 1>(
            std::tuple<Out...,
                       std::pair<std::false_type,
                                 std::integral_constant<std::size_t, Pos2>>>{},
            nodes);
    } else {
        using Id1 = std::tuple_element_t<0, std::tuple_element_t<Pos1, Tuple1>>;
        using Id2 = std::tuple_element_t<0, std::tuple_element_t<Pos2, Tuple2>>;

        static_assert(!are_same(Id1{}, Id2{}));

        if constexpr (less_than2(Id1{}, Id2{}, nodes)) {
            return merge_sorted_aux2<Tuple1, Tuple2, Pos1, Pos2 + 1>(
                std::tuple<
                    Out...,
                    std::pair<std::false_type,
                              std::integral_constant<std::size_t, Pos2>>>{},
                nodes);
        } else {
            return merge_sorted_aux2<Tuple1, Tuple2, Pos1 + 1, Pos2>(
                std::tuple<Out...,
                           std::pair<std::true_type, std::integral_constant<
                                                         std::size_t, Pos1>>>{},
                nodes);
        }
    }
}

template <class Tuple1, class Tuple2, class... Flag, class... Position>
constexpr auto
merge_sorted_aux4(std::tuple<std::pair<Flag, Position>...> /* positions */) {
    // why do we do a conditional on the tuple element index?
    // because conditional_t is eager on evaluation and we therefore have to
    // refer to a type that always exists..
    return std::tuple<std::conditional_t<
        Flag::value,
        std::tuple_element_t<
            std::conditional_t<Flag::value, Position,
                               std::integral_constant<std::size_t, 0>>::value,
            Tuple1>,
        std::tuple_element_t<
            std::conditional_t<Flag::value,
                               std::integral_constant<std::size_t, 0>,
                               Position>::value,
            Tuple2>>...>{};
}

template <class Tuple1, class Tuple2, class Nodes>
constexpr auto merge_sorted(Tuple1 in1, Tuple2 in2, Nodes nodes) {
    static_assert(is_sorted(in1, nodes));
    static_assert(is_sorted(in2, nodes));

    if constexpr (std::tuple_size_v<Tuple1> == 0) {
        return in2;
    } else if constexpr (std::tuple_size_v<Tuple2> == 0) {
        return in1;
    } else {
        constexpr auto position =
            merge_sorted_aux2<Tuple1, Tuple2>(std::tuple<>{}, nodes);
        return merge_sorted_aux4<Tuple1, Tuple2>(position);
    }
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

template <std::size_t Pos1 = 0, std::size_t Pos2 = 0, class Tuple1,
          class Tuple2, class Out, class Nodes>
constexpr auto multiply_ordered_aux(Tuple1 in1, Tuple2 in2, Out out,
                                    Nodes nodes) {
    if constexpr (std::tuple_size_v<Tuple1> == Pos1 &&
                  std::tuple_size_v<Tuple2> == Pos2) {
        return out;
    } else if constexpr (std::tuple_size_v<Tuple2> == Pos2) {
        return std::tuple_cat(
            out, sub_tuple<Pos1, std::tuple_size_v<Tuple1> - 1>(in1));
    } else if constexpr (std::tuple_size_v<Tuple1> == Pos1) {
        return std::tuple_cat(
            out, sub_tuple<Pos2, std::tuple_size_v<Tuple2> - 1>(in2));
    } else {
        constexpr auto element1 = std::get<Pos1>(in1);
        constexpr auto element2 = std::get<Pos2>(in2);
        constexpr auto id1 = get_id(element1);
        constexpr auto id2 = get_id(element2);
        constexpr auto idx1 = find<Nodes, decltype(id1)>();
        constexpr auto idx2 = find<Nodes, decltype(id2)>();

        if constexpr (idx1 < idx2) {
            constexpr std::size_t power = get_power(element1);
            return multiply_ordered_aux<Pos1 + 1, Pos2>(
                in1, in2, std::tuple_cat(out, d<power>(id1)), nodes);
        } else if constexpr (idx1 > idx2) {
            constexpr std::size_t power = get_power(element2);
            return multiply_ordered_aux<Pos1, Pos2 + 1>(
                in1, in2, std::tuple_cat(out, d<power>(id2)), nodes);
        } else {
            constexpr std::size_t power =
                get_power(element1) + get_power(element2);
            return multiply_ordered_aux<Pos1 + 1, Pos2 + 1>(
                in1, in2, std::tuple_cat(out, d<power>(id1)), nodes);
        }
    }
}

template <std::size_t Pos = 0, class Tuple, class Nodes>
constexpr auto is_sorted_diff_op(Tuple in, Nodes nodes) -> bool {
    if constexpr (std::tuple_size_v<Tuple> <= (Pos + 1)) {
        return true;
    } else {
        constexpr auto in_first = std::get<Pos>(in);
        constexpr auto in_second = std::get<Pos + 1>(in);
        constexpr auto id1 = get_id(in_first);
        constexpr auto id2 = get_id(in_second);
        constexpr auto idx1 = find<Nodes, decltype(id1)>();
        constexpr auto idx2 = find<Nodes, decltype(id2)>();
        if constexpr (idx1 < idx2) {
            return true;
        } else {
            return is_sorted_diff_op<Pos + 1>(in, nodes);
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

template <class... DiffOps, std::size_t... I>
constexpr auto add_position_aux(std::tuple<DiffOps...> /* tuple */,
                                std::index_sequence<I...> /* idxs */) {
    return std::tuple<std::pair<
        DiffOps, std::pair<on_interface_t,
                           std::integral_constant<std::size_t, I>>>...>{};
}

template <class TupleDiffOps> constexpr auto add_position(TupleDiffOps tuple) {
    return add_position_aux(
        tuple, std::make_index_sequence<std::tuple_size_v<TupleDiffOps>>{});
}

template <class... DiffOps, class... Position>
constexpr auto add_position(std::tuple<DiffOps...> /* diff_op */,
                            std::tuple<Position...> /* pos */) {
    static_assert(sizeof...(DiffOps) == sizeof...(Position));
    return std::tuple<std::pair<DiffOps, Position>...>{};
}

} // namespace adhoc4::detail

#endif // ADHOC4_BACKPROPAGATOR_TOOLS_HPP
