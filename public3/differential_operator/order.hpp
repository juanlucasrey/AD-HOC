/*
 * AD-HOC, Automatic Differentiation for High Order Calculations
 *
 * This file is part of the AD-HOC distribution
 * (https://github.com/juanlucasrey/AD-HOC).
 * Copyright (c) 2024 Juan Lucas Rey
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
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
