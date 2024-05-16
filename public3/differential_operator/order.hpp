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

#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_ORDER_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_ORDER_HPP

#include "differential_operator.hpp"

namespace adhoc3 {

constexpr auto max(auto const &value, auto const &...args) {
    if constexpr (sizeof...(args) == 0U) {
        return value;
    } else {
        const auto max_val = max(args...);
        return value > max_val ? value : max_val;
    }
}

template <class... Ids, std::size_t... Orders, std::size_t... Powers>
constexpr auto
max_order(std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...> /* id1 */) {
    return max(Orders...);
}

template <class... Ops>
constexpr auto max_orders(std::tuple<Ops...> /* id1 */) {
    return max(max_order(Ops{})...);
}

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_ORDER_HPP
