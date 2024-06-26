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

#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_HPP

// #include "adhoc.hpp"
#include "dependency.hpp"
#include "utils.hpp"

namespace adhoc3 {

namespace der {

template <class Id, std::size_t Order = 1> struct d {};
template <class Id, std::size_t Power = 1> struct p {};

} // namespace der

namespace detail {

template <class IdInput, std::size_t... PowersNode, std::size_t... OrdersNode,
          class... IdsNode>
constexpr auto
der_non_null_aux(std::tuple<der::p<der::d<IdsNode, OrdersNode>, PowersNode>...>)
    -> std::size_t {
    return (detail::sum_no_overflow<order<IdsNode, IdInput>::call()...>());
}

} // namespace detail

template <std::size_t... PowersNode, std::size_t... OrdersNode,
          class... IdsNode, std::size_t... PowersInput, class... IdsInput>
constexpr auto
der_non_null(std::tuple<der::p<der::d<IdsNode, OrdersNode>, PowersNode>...> n,
             std::tuple<der::p<der::d<IdsInput, 1>, PowersInput>...>) -> bool {
    return ((PowersInput <= detail::der_non_null_aux<IdsInput>(n)) && ...);
}

namespace detail {

template <class In, class... Ids, std::size_t... Powers, std::size_t... Orders>
constexpr auto
create_single_tuple(In /* in */,
                    std::tuple<der::p<der::d<Ids, Orders>, Powers>...> o) {
    if constexpr (has_type<In, Ids...>()) {
        constexpr auto idx = idx_type2<In, Ids...>();
        constexpr auto v = std::get<idx>(o);
        return std::tuple<decltype(v)>{};
    } else {
        return std::tuple<>{};
    }
};

} // namespace detail

template <class... Ids, std::size_t... Powers, std::size_t... Orders,
          class... CalculationNodes>
constexpr auto order_differential_operator(
    std::tuple<der::p<der::d<Ids, Orders>, Powers>...> diff_operator,
    std::tuple<CalculationNodes...> nodes) {
    return std::apply(
        [diff_operator](auto &&...args) {
            return std::tuple_cat(
                detail::create_single_tuple(args, diff_operator)...);
        },
        nodes);
}

template <std::size_t Order = 1, class Id> constexpr auto d(Id /* id */) {
    return std::tuple<der::p<der::d<Id, 1>, Order>>{};
}

template <std::size_t... Orders1, class... Ids1, std::size_t... Orders2,
          class... Ids2>
constexpr auto operator*(std::tuple<der::p<der::d<Ids1, 1>, Orders1>...> id1,
                         std::tuple<der::p<der::d<Ids2, 1>, Orders2>...> id2) {
    return std::tuple_cat(id1, id2);
}

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_HPP
