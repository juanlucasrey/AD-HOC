/*
 * AD-HOC, Automatic Differentiation for High Order Calculations
 *
 * This file is part of the AD-HOC distribution
 * (https://github.com/juanlucasrey/adhoc2).
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

#ifndef ADHOC3_DERS_HPP
#define ADHOC3_DERS_HPP

#include "adhoc.hpp"
#include "dependency.hpp"
#include "strict_order.hpp"

namespace adhoc3 {

template <class... Types> struct mvar_t {};
template <class Id, std::size_t Order> struct IdAndOrder_t {};

namespace detail {

template <class... Args> auto duplicate_id(Args...) -> bool;

template <> auto constexpr duplicate_id() -> bool { return false; }

template <std::size_t Order1, class IdDerived, std::size_t... Orders,
          class... IdsDerived>
auto constexpr duplicate_id(mvar_t<IdAndOrder_t<IdDerived, Order1>>,
                            mvar_t<IdAndOrder_t<IdsDerived, Orders>>...)
    -> bool {
    bool thisisduplicate = (std::is_same_v<IdDerived, IdsDerived> || ...);
    bool duplicateinnext =
        duplicate_id(mvar_t<IdAndOrder_t<IdsDerived, Orders>>{}...);
    return thisisduplicate || duplicateinnext;
}

template <std::size_t... Orders, class... IdsDerived>
auto var_aux(mvar_t<mvar_t<IdAndOrder_t<IdsDerived, Orders>>...>) {

    static_assert(
        !detail::duplicate_id(mvar_t<IdAndOrder_t<IdsDerived, Orders>>{}...),
        "duplicate ids on multivariate derivative declaration");

    return mvar_t<IdAndOrder_t<IdsDerived, Orders>...>{};
}

template <class... Ts>
using MyContainer = instantiate_t<mvar_t, sorted_list_t<list<Ts...>>>;

} // namespace detail

template <std::size_t Order = 1, class Derived> auto var(Derived) {
    return mvar_t<IdAndOrder_t<Derived, Order>>{};
}

template <std::size_t... Orders, class... IdsDerived>
auto var(mvar_t<IdAndOrder_t<IdsDerived, Orders>>...) {
    return detail::var_aux(
        detail::MyContainer<mvar_t<IdAndOrder_t<IdsDerived, Orders>>...>{});
}

template <class NumIdDerived, std::size_t... Orders, class... IdsDerived>
constexpr auto der_non_null(NumIdDerived,
                            mvar_t<IdAndOrder_t<IdsDerived, Orders>...>)
    -> bool {
    return ((Orders <= order<NumIdDerived, IdsDerived>::call()) && ...);
}

// namespace der {

// template <class Id, std::size_t Order = 1> struct d {};
// template <class Id, std::size_t Power = 1> struct p {};
// template <class... Ids> struct m {};

// } // namespace der

// namespace detail {

// template <class IdInput, std::size_t... PowersNode, std::size_t...
// OrdersNode,
//           class... IdsNode>
// constexpr auto
// der_non_null_aux(der::m<der::p<der::d<IdsNode, OrdersNode>, PowersNode>...>)
//     -> std::size_t {
//     return (detail::sum_no_overflow<order<IdsNode, IdInput>::call()...>());
// }

// } // namespace detail

// template <std::size_t... PowersNode, std::size_t... OrdersNode,
//           class... IdsNode, std::size_t... PowersInput, class... IdsInput>
// constexpr auto
// der_non_null(der::m<der::p<der::d<IdsNode, OrdersNode>, PowersNode>...> n,
//              der::m<der::p<der::d<IdsInput, 1>, PowersInput>...>) -> bool {
//     return ((PowersInput <= detail::der_non_null_aux<IdsInput>(n)) && ...);
// }

} // namespace adhoc3

#endif // ADHOC3_DERS_HPP
