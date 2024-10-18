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

#include <concepts>
#include <cstddef>
#include <tuple>
#include <type_traits>

namespace adhoc4 {

template <class Tuple, class T, std::size_t N = 0>
constexpr auto find() -> std::size_t {
    if constexpr (N == std::tuple_size_v<Tuple>) {
        return N;
    } else if constexpr (std::is_same_v<const std::tuple_element_t<N, Tuple>,
                                        const T>) {
        return N;
    } else {
        return find<Tuple, T, N + 1>();
    }
}

template <class... Ts, class T>
constexpr auto contains(std::tuple<Ts...> /* tuple */, T /* value */) -> bool {
    return (std::is_same_v<const Ts, const T> || ...);
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

template <class... Types, bool... Flags>
constexpr auto filter(std::tuple<Types...> /* in */,
                      std::integer_sequence<bool, Flags...> /* calc_flags */) {
    return std::tuple_cat(
        std::conditional_t<Flags, std::tuple<Types>, std::tuple<>>{}...);
}

template <class... Types, bool... Flags>
constexpr auto
separate(std::tuple<Types...> /* in */,
         std::integer_sequence<bool, Flags...> /* calc_flags */) {
    return std::make_pair(
        std::tuple_cat(
            std::conditional_t<Flags, std::tuple<Types>, std::tuple<>>{}...),
        std::tuple_cat(
            std::conditional_t<!Flags, std::tuple<Types>, std::tuple<>>{}...));
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
constexpr auto replace_first_aux(std::tuple<> /* tuple */, Old /* old_value */,
                                 New /* new_value */, Res res) {
    return res;
}

template <class Input, class Old, class New, class Res>
constexpr auto replace_first_aux(Input tuple, Old old_value, New new_value,
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
constexpr auto replace_first(Tuple tuple, Old old_value, New new_value) {
    return detail::replace_first_aux(tuple, old_value, new_value,
                                     std::tuple<>{});
}

namespace detail {

template <class Tuple, std::size_t... I>
constexpr auto sub_tuple_aux(std::index_sequence<I...> /* i */) {
    return std::tuple<std::tuple_element_t<I, Tuple>...>{};
}

} // namespace detail

template <std::size_t From, std::size_t To, class Tuple>
constexpr auto sub_tuple(Tuple /* tuple */) {
    static_assert(std::tuple_size_v<Tuple> > From);
    return detail::sub_tuple_aux<Tuple>(make_index_sequence<From, To>());
}

} // namespace adhoc4

#endif // ADHOC4_UTILS_TUPLE_HPP
