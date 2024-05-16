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

#ifndef ADHOC3_TUPLE_UTILS_HPP
#define ADHOC3_TUPLE_UTILS_HPP

#include "utils.hpp"

#include <iostream>
#include <tuple>
#include <type_traits>
#include <utility>

namespace adhoc3 {

namespace detail {

template <std::size_t IdxToReplace, class Replacement, class... TupleValues,
          std::size_t... IndicesStart, std::size_t... IndicesEnd>
auto constexpr replace_first_aux2(
    std::tuple<TupleValues...> /* in */,
    std::index_sequence<IndicesStart...> /* idx1 */,
    std::index_sequence<IndicesEnd...> /* idx2 */) {
    return std::make_tuple(
        std::tuple_element_t<IndicesStart, std::tuple<TupleValues...>>{}...,
        Replacement{},
        std::tuple_element_t<IndicesEnd + IdxToReplace + 1,
                             std::tuple<TupleValues...>>{}...);
}

template <std::size_t IdxToReplace, class Replacement, class... TupleValues,
          class IndicesStart = std::make_index_sequence<IdxToReplace>,
          class IndicesEnd = std::make_index_sequence<sizeof...(TupleValues) -
                                                      IdxToReplace - 1>>
auto constexpr replace_first_aux(std::tuple<TupleValues...> /* in */) {
    return replace_first_aux2<IdxToReplace, Replacement>(
        std::tuple<TupleValues...>{}, IndicesStart{}, IndicesEnd{});
}

} // namespace detail

template <class Type, class... TupleValues>
auto constexpr get_idx_first(std::tuple<TupleValues...> /* in */) {
    static_assert(has_type<Type, TupleValues...>(),
                  "type to replace not on this tuple");
    return detail::get_index<Type, TupleValues...>::value;
}

template <class TypeToReplace, class Replacement, class... TupleValues>
auto constexpr replace_first(std::tuple<TupleValues...> /* in */) {
    static_assert(has_type<TypeToReplace, TupleValues...>(),
                  "type to replace not on this tuple");
    constexpr auto idx =
        detail::get_index<TypeToReplace, TupleValues...>::value;
    return detail::replace_first_aux<idx, Replacement>(
        std::tuple<TupleValues...>{});
}

template <typename T, size_t N> class generate_tuple_type {
    template <typename = std::make_index_sequence<N>> struct impl;

    template <size_t... Is> struct impl<std::index_sequence<Is...>> {
        template <size_t> using wrap = T;

        using type = std::tuple<wrap<Is>...>;
    };

  public:
    using type = typename impl<>::type;
};

template <class Type, class... TupleValues>
auto constexpr get_idx_first2(std::tuple<TupleValues...> /* in */) {
    if constexpr (has_type<Type, TupleValues...>()) {
        return detail::get_index<Type, TupleValues...>::value;
    } else {
        // if the type is not found we go beyond the end
        return sizeof...(TupleValues);
    }
}

template <std::size_t IdxToReplace, class Replacement, class... TupleValues>
auto constexpr replace_first2(std::tuple<TupleValues...> /* in */) {
    if constexpr (IdxToReplace >= sizeof...(TupleValues)) {
        return std::tuple<TupleValues...>{};
    } else {
        return detail::replace_first_aux<IdxToReplace, Replacement>(
            std::tuple<TupleValues...>{});
    }
}

template <typename T, typename Tuple> struct has_type_tuple;

template <typename T, typename... Us>
struct has_type_tuple<T, std::tuple<Us...>>
    : std::disjunction<std::is_same<T, Us>...> {};

template <typename T, typename Tuple>
constexpr bool has_type_tuple_v = has_type_tuple<T, Tuple>::value;

template <class... Ts, class T>
constexpr auto contains(std::tuple<Ts...> /* tuple */, T /* value */) -> bool {
    return std::disjunction<std::is_same<Ts, T>...>::value;
}

template <class... Ts, class T>
auto constexpr get_idx(std::tuple<Ts...> tuple, T value) -> std::size_t {
    static_assert(contains(tuple, value));
    // constexpr std::array<bool, sizeof...(Ts)> a{{std::is_same_v<T, Ts>...}};
    // return std::distance(a.begin(), std::find(a.begin(), a.end(), true));
    return detail::get_index<T, Ts...>::value;
}

template <class... Ts, class Old, class New>
auto constexpr replace(std::tuple<Ts...> tuple, Old old_value,
                       New /* new_value */) {
    if constexpr (contains(tuple, old_value)) {
        constexpr auto idx_to_replace = get_idx(tuple, old_value);
        return detail::replace_first_aux<idx_to_replace, New>(
            std::tuple<Ts...>{});
    } else {
        return tuple;
    }
}

namespace detail {

template <class ThisVal, class RemoveVal>
constexpr auto remove_aux(ThisVal val, RemoveVal /* rem */) {

    if constexpr (std::is_same_v<ThisVal, RemoveVal>) {
        return std::tuple<>{};
    } else {
        return std::make_tuple(val);
    }
}

} // namespace detail

template <typename T, typename... Us>
constexpr auto remove(std::tuple<Us...> tuple, T value) {
    return std::apply(
        [value](auto... val) {
            return std::tuple_cat(detail::remove_aux(val, value)...);
        },
        tuple);
}

template <typename First, typename... Rest>
constexpr auto head(std::tuple<First, Rest...> /* tuple */) {
    return First{};
}

template <typename First, typename... Rest>
constexpr auto tail(std::tuple<First, Rest...> /* tuple */) {
    return std::tuple<Rest...>{};
}

template <typename Tuple>
constexpr auto size(Tuple /* tuple */) -> std::size_t {
    return std::tuple_size_v<Tuple>;
}

namespace detail {

template <std::size_t N, class UnfilteredTuple, class CalculateFlags,
          class Output>
constexpr auto filter_aux(UnfilteredTuple in, CalculateFlags calc_flags,
                          Output out) {
    if constexpr (N == std::tuple_size_v<UnfilteredTuple>) {
        return out;
    } else {
        if constexpr (std::get<N>(calc_flags)) {
            return filter_aux<N + 1>(
                in, calc_flags,
                std::tuple_cat(out, std::make_tuple(std::get<N>(in))));

        } else {
            return filter_aux<N + 1>(in, calc_flags, out);
        }
    }
}

} // namespace detail

template <class UnfilteredTuple, class CalculateFlags>
constexpr auto filter(UnfilteredTuple in, CalculateFlags calc_flags) {
    return detail::filter_aux<0>(in, calc_flags, std::tuple<>{});
}

template <std::size_t Idx, std::size_t... I>
constexpr auto get(std::index_sequence<I...>) {
    return std::array<std::size_t, sizeof...(I)>{I...}[Idx];
}

namespace detail {

template <std::size_t N, class UnfilteredIS, class CalculateFlags,
          std::size_t... Out>
constexpr auto filter_is_aux(UnfilteredIS in, CalculateFlags calc_flags,
                             std::index_sequence<Out...> out) {
    if constexpr (N == in.size()) {
        return out;
    } else {
        if constexpr (std::get<N>(calc_flags)) {
            return filter_is_aux<N + 1>(
                in, calc_flags, std::index_sequence<Out..., get<N>(in)>{});

        } else {
            return filter_is_aux<N + 1>(in, calc_flags, out);
        }
    }
}

} // namespace detail

template <class UnfilteredIS, class CalculateFlags>
constexpr auto filter_is(UnfilteredIS in, CalculateFlags calc_flags) {
    return detail::filter_is_aux<0>(in, calc_flags, std::index_sequence<>{});
}

template <std::size_t N, class IntegerSequences, class Output>
constexpr auto precedent_required_aux(IntegerSequences integer_sequences,
                                      Output out) {
    if constexpr (N == std::get<0>(integer_sequences).size()) {
        return out;
    } else {
        constexpr auto result_filtered = std::apply(
            [](auto... sequence) { return ((get<N>(sequence) + ...) != 0); },
            integer_sequences);

        using this_result_type =
            std::conditional<result_filtered, std::true_type,
                             std::false_type>::type;

        return precedent_required_aux<N + 1>(
            integer_sequences,
            std::tuple_cat(out, std::make_tuple(this_result_type{})));
    }
}

template <class IntegerSequences>
constexpr auto precedent_required(IntegerSequences integer_sequences) {
    return precedent_required_aux<0>(integer_sequences, std::tuple<>{});
}

template <std::size_t... I>
constexpr auto sum(std::index_sequence<I...> /* idx_seq */) -> std::size_t {
    return (I + ...);
}
} // namespace adhoc3

#endif // ADHOC3_TUPLE_UTILS_HPP
