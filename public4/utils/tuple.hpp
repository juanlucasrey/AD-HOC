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

#ifndef ADHOC4_UTILS_TUPLE_HPP
#define ADHOC4_UTILS_TUPLE_HPP

#include <array>
#include <cmath>
#include <concepts>
#include <tuple>

namespace adhoc4 {

namespace detail {

template <std::size_t Idx = 0, class T>
auto constexpr get_first_type_idx_aux(std::tuple<> /* tuple */, T /* value */)
    -> std::size_t {
    return Idx;
}

template <std::size_t Idx = 0, class FirstT, class... Ts, class T>
auto constexpr get_first_type_idx_aux(std::tuple<FirstT, Ts...> /* tuple */,
                                      T /* value */) -> std::size_t {
    if constexpr (std::same_as<const FirstT, const T>) {
        return Idx;
    } else {
        return get_first_type_idx_aux<Idx + 1>(std::tuple<Ts...>{}, T{});
    }
}

} // namespace detail

template <class... Ts, class T>
auto constexpr get_first_type_idx(std::tuple<Ts...> /* tuple */, T /* value */)
    -> std::size_t {
    return detail::get_first_type_idx_aux(std::tuple<Ts...>{}, T{});
}

template <class... Ts, class T>
auto constexpr contains(std::tuple<Ts...> tuple, T value) -> bool {
    std::size_t constexpr pos = get_first_type_idx(tuple, value);
    return pos < sizeof...(Ts);
}

template <typename T, typename... Us>
constexpr auto remove(std::tuple<Us...> /* tuple */, T /* value */) {
    return std::tuple_cat(
        std::conditional_t<std::is_same_v<const T, const Us>, std::tuple<>,
                           std::tuple<Us>>{}...);
}

template <typename First, typename... Rest>
constexpr auto head(std::tuple<First, Rest...> /* tuple */) {
    return First{};
}

template <typename First, typename... Rest>
constexpr auto tail(std::tuple<First, Rest...> /* tuple */) {
    return std::tuple<Rest...>{};
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

namespace detail {

template <std::size_t N, class UnfilteredTuple, class CalculateFlags,
          class Output1, class Output2>
constexpr auto separate_aux(UnfilteredTuple in, CalculateFlags calc_flags,
                            Output1 out1, Output2 out2) {
    if constexpr (N == std::tuple_size_v<UnfilteredTuple>) {
        return std::make_pair(out1, out2);
    } else {
        if constexpr (std::get<N>(calc_flags)) {
            return separate_aux<N + 1>(
                in, calc_flags,
                std::tuple_cat(out1, std::make_tuple(std::get<N>(in))), out2);

        } else {
            return separate_aux<N + 1>(
                in, calc_flags, out1,
                std::tuple_cat(out2, std::make_tuple(std::get<N>(in))));
        }
    }
}

} // namespace detail

template <class UnfilteredTuple, class CalculateFlags>
constexpr auto separate(UnfilteredTuple in, CalculateFlags calc_flags) {
    return detail::separate_aux<0>(in, calc_flags, std::tuple<>{},
                                   std::tuple<>{});
}

namespace detail {
template <class Out> constexpr auto invert_aux(std::tuple<> /* in */, Out out) {
    return out;
}

template <class Input, class... Inputs, class... Out>
constexpr auto invert_aux(std::tuple<Input, Inputs...> /* in */,
                          std::tuple<Out...> /* out */) {
    return invert_aux(std::tuple<Inputs...>{}, std::tuple<Input, Out...>{});
}
} // namespace detail

template <class Tuple> constexpr auto invert(Tuple in) {
    return detail::invert_aux(in, std::tuple<>{});
}

template <typename Tuple>
constexpr auto size(Tuple /* tuple */) -> std::size_t {
    return std::tuple_size_v<Tuple>;
}

namespace detail {

template <class Old, class New, class Res>
auto constexpr replace_first_aux(std::tuple<> /* tuple */, Old /* old_value */,
                                 New /* new_value */, Res res) {
    return res;
}

template <class Input, class Old, class New, class Res>
auto constexpr replace_first_aux(Input tuple, Old old_value, New new_value,
                                 Res res) {
    if constexpr (std::is_same_v<decltype(head(tuple)), Old>) {
        return std::tuple_cat(res, std::make_tuple(new_value), tail(tuple));
    } else {
        return replace_first_aux(
            tail(tuple), old_value, new_value,
            std::tuple_cat(res, std::make_tuple(head(tuple))));
    }
}

} // namespace detail

template <class Tuple, class Old, class New>
auto constexpr replace_first(Tuple tuple, Old old_value, New new_value) {
    return detail::replace_first_aux(tuple, old_value, new_value,
                                     std::tuple<>{});
}

} // namespace adhoc4

#endif // ADHOC4_UTILS_TUPLE_HPP
