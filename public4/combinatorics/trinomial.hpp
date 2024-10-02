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

#ifndef ADHOC4_PARTITION_TRINOMIAL_HPP
#define ADHOC4_PARTITION_TRINOMIAL_HPP

#include <array>
#include <cstddef>
#include <tuple>

namespace adhoc4 {

namespace detail {

constexpr auto NextTrinomial(std::array<std::size_t, 3> const &prev) {
    // check the previous trinomial coefficient is not the last
    std::array<std::size_t, 3> arrlast{0};
    arrlast.back() = prev.back();
    if (arrlast == prev) {
        return prev;
    }

    std::array<std::size_t, 3> arr = prev;

    if (arr[0] == 0) {
        if (arr[1] != 0) {
            std::size_t rem_val = arr[1];
            arr[1] = 0;
            arr[0] = rem_val - 1;
            arr[2]++;
        }
    } else {
        arr[0]--;
        arr[1]++;
    }

    return arr;
}

template <std::size_t... I, std::size_t... I2>
constexpr auto NextTrinomialIS_aux(std::index_sequence<I...> /* i */,
                                   std::index_sequence<I2...> /* i2 */) {
    constexpr std::array<std::size_t, sizeof...(I)> temp{I...};
    constexpr auto next = NextTrinomial(temp);
    return std::index_sequence<next[I2]...>{};
}

template <std::size_t... I>
constexpr auto NextTrinomialIS(std::index_sequence<I...> i) {
    return detail::NextTrinomialIS_aux(
        i, std::make_index_sequence<sizeof...(I)>{});
}

template <class CurrentIS, class Output>
constexpr auto TrinomialSequences_aux(CurrentIS const current, Output out) {
    auto constexpr next = NextTrinomialIS(current);
    if constexpr (!std::is_same_v<decltype(current), decltype(next)>) {
        return TrinomialSequences_aux(
            next, std::tuple_cat(out, std::make_tuple(current)));
    } else {
        return std::tuple_cat(out, std::make_tuple(current));
    }
}

} // namespace detail

template <std::size_t Order> constexpr auto TrinomialSequences() {
    constexpr auto first = std::index_sequence<Order, 0, 0>{};
    return detail::TrinomialSequences_aux(first, std::tuple<>{});
}

namespace detail {

template <std::size_t DoublePower>
constexpr auto NextTrinomial2(std::array<std::size_t, 3> const &prev) {
    // check the previous trinomial coefficient is not the last
    std::array<std::size_t, 3> arrlast{0};
    arrlast.back() = prev.back();
    if (arrlast == prev) {
        return prev;
    }

    std::array<std::size_t, 3> arr = prev;

    if (arr[0] == 0) {
        if (arr[1] != 0) {
            std::size_t rem_val = arr[1];
            rem_val--;
            arr[2]++;

            std::size_t move = rem_val >= DoublePower ? DoublePower : rem_val;
            arr[0] = move;
            arr[1] = rem_val - move;
        }
    } else {
        arr[0]--;
        arr[1]++;
    }

    return arr;
}

template <std::size_t DoublePower, std::size_t... I, std::size_t... I2>
constexpr auto NextTrinomialIS_aux2(std::index_sequence<I...> /* i */,
                                    std::index_sequence<I2...> /* i2 */) {
    constexpr std::array<std::size_t, sizeof...(I)> temp{I...};
    constexpr auto next = NextTrinomial2<DoublePower>(temp);
    return std::index_sequence<next[I2]...>{};
}

template <std::size_t DoublePower, std::size_t... I>
constexpr auto NextTrinomialIS2(std::index_sequence<I...> i) {
    return NextTrinomialIS_aux2<DoublePower>(
        i, std::make_index_sequence<sizeof...(I)>{});
}

template <std::size_t DoublePower, class CurrentIS, class Output>
constexpr auto TrinomialSequences_aux2(CurrentIS const current, Output out) {
    auto constexpr next = NextTrinomialIS2<DoublePower>(current);
    if constexpr (!std::is_same_v<decltype(current), decltype(next)>) {
        return TrinomialSequences_aux2<DoublePower>(
            next, std::tuple_cat(out, std::make_tuple(current)));
    } else {
        return std::tuple_cat(out, std::make_tuple(current));
    }
}

} // namespace detail

// this is a highly tailored method that gives the valid monomial expansions of
// (e(x)*e(y) + e(x) + e(y))^Order that do not have a power going beyond
// MaxOrder
template <std::size_t Order, std::size_t MaxOrder>
constexpr auto TrinomialSequencesMult() {
    static_assert(MaxOrder >= Order);
    constexpr auto AllowedDoublePower = MaxOrder - Order;
    constexpr auto DoublePower =
        AllowedDoublePower > Order ? Order : AllowedDoublePower;

    constexpr auto first =
        std::index_sequence<DoublePower, Order - DoublePower, 0>{};
    return detail::TrinomialSequences_aux2<DoublePower>(first, std::tuple<>{});
}

} // namespace adhoc4

#endif // ADHOC4_PARTITION_TRINOMIAL_HPP
