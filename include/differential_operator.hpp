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

#ifndef DIFFERENTIAL_OPERATOR_HPP
#define DIFFERENTIAL_OPERATOR_HPP

#include <cstddef>
#include <tuple>

#include "combinatorics/factorial.hpp"
#include "utils/index_sequence.hpp"
#include "utils/tuple.hpp"

namespace adhoc {

namespace der {

template <std::size_t Order, class Id> struct d {};

} // namespace der

template <std::size_t Order = 1UL, class Id> constexpr auto d(Id /* id */) {
    if constexpr (Order == 0) {
        return std::tuple<>{};
    } else {
        return std::tuple<der::d<Order, Id>>{};
    }
}

namespace detail {

template <class Id, std::size_t Order1, std::size_t Order2>
constexpr auto add_duplicate(der::d<Order1, Id> /* in1 */,
                             der::d<Order2, Id> /* in2 */) {
    return der::d<Order1 + Order2, Id>{};
}

template <class Id1, std::size_t Order1, class... Ids2, std::size_t... Orders2>
constexpr auto add_duplicates(der::d<Order1, Id1> in,
                              std::tuple<der::d<Orders2, Ids2>...> id2) {
    constexpr auto idx = find<std::tuple<Ids2...>, Id1>();

    if constexpr (idx < sizeof...(Orders2)) {
        constexpr auto v = std::get<idx>(id2);
        constexpr auto added = add_duplicate(in, v);
        return std::make_tuple(added);
    } else {
        return std::make_tuple(in);
    }
};

template <class Id1, std::size_t Order1, class... Ids2, std::size_t... Orders2>
constexpr auto
remove_duplicates(der::d<Order1, Id1> in,
                  std::tuple<der::d<Orders2, Ids2>...> /* id2 */) {
    constexpr auto idx = find<std::tuple<Ids2...>, Id1>();

    if constexpr (idx < sizeof...(Orders2)) {
        return std::tuple{};
    } else {
        return std::make_tuple(in);
    }
};

} // namespace detail

template <class... Ids1, std::size_t... Orders1, class... Ids2,
          std::size_t... Orders2>
constexpr auto operator*(std::tuple<der::d<Orders1, Ids1>...> id1,
                         std::tuple<der::d<Orders2, Ids2>...> id2) {
    auto augmented_in1 = std::apply(
        [id2](auto... id1s) {
            return std::tuple_cat(detail::add_duplicates(id1s, id2)...);
        },
        id1);

    auto reduced_in2 = std::apply(
        [id1](auto... id2s) {
            return std::tuple_cat(detail::remove_duplicates(id2s, id1)...);
        },
        id2);

    return std::tuple_cat(augmented_in1, reduced_in2);
}

template <std::size_t Order, class Id>
constexpr auto get_power(const der::d<Order, Id> /* in */) -> std::size_t {
    return Order;
}

template <std::size_t Order, class Id>
constexpr auto get_id(const der::d<Order, Id> /* in */) {
    return Id{};
}

namespace detail {
template <std::size_t P, class... Ids1, std::size_t... Orders1>
constexpr auto
elevate_by_power(std::tuple<der::d<Orders1, Ids1>...> /* id1 */) {
    if constexpr (P == 0) {
        return std::tuple<>{};
    } else {
        return std::tuple<der::d<Orders1 * P, Ids1>...>{};
    }
}
} // namespace detail

template <class... DiffOps, std::size_t... I>
constexpr auto
create_differential_operator(std::tuple<DiffOps...> /* vars */,
                             std::index_sequence<I...> /* is */) {
    return (detail::elevate_by_power<I>(DiffOps{}) * ...);
}

template <class... Ids, std::size_t... Orders>
constexpr auto multiplicity(std::tuple<der::d<Orders, Ids>...> /* id1 */)
    -> double {
    return (factorial(Orders) * ...);
}

template <class... Ids, std::size_t... Orders>
constexpr auto power(std::tuple<der::d<Orders, Ids>...> /* id1 */)
    -> std::size_t {
    if constexpr (sizeof...(Orders)) {
        return (Orders + ...);
    } else {
        return 0;
    }
}

template <class... Ids1, std::size_t... Orders1>
constexpr auto is_one(std::tuple<der::d<Orders1, Ids1>...> /* id1 */)
    -> double {
    return all_equal_to<1>(std::index_sequence<Orders1...>{});
}

} // namespace adhoc

#endif // DIFFERENTIAL_OPERATOR_HPP
