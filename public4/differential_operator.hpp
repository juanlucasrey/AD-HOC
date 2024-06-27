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

#ifndef ADHOC4_DIFFERENTIAL_OPERATOR_HPP
#define ADHOC4_DIFFERENTIAL_OPERATOR_HPP

#include <cstddef>
#include <tuple>

#include "utils/tuple.hpp"

namespace adhoc4 {

namespace der {

template <std::size_t Order, class Id> struct d {};

} // namespace der

template <std::size_t Order = 1UL, class Id> constexpr auto d(Id /* id */) {
    return std::tuple<der::d<Order, Id>>{};
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
    constexpr auto idx = get_first_type_idx(id2, in);

    if constexpr (idx < sizeof...(Orders2)) {
        constexpr auto v = std::get<idx>(id2);
        constexpr auto added = add_duplicate(in, v);
        return std::make_tuple(added);
    } else {
        return std::make_tuple(in);
    }
};

template <class Id1, std::size_t Order1, class... Ids2, std::size_t... Orders2>
constexpr auto remove_duplicates(der::d<Order1, Id1> in,
                                 std::tuple<der::d<Orders2, Ids2>...> id2) {

    constexpr auto idx = get_first_type_idx(id2, in);

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
} // namespace adhoc4

#endif // ADHOC4_DIFFERENTIAL_OPERATOR_HPP