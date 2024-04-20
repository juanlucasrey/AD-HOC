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

#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_IS_SORTED_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_IS_SORTED_HPP

#include "../adhoc.hpp"
#include "../dependency.hpp"
#include "differential_operator.hpp"
#include "less_than.hpp"

namespace adhoc3 {

namespace detail {

template <class Nodes>
constexpr auto is_ordered_aux(Nodes /* nodes */, std::tuple<> /* in */) {
    return true;
}

template <class Nodes, class Id>
constexpr auto is_ordered_aux(Nodes /* nodes */, std::tuple<Id> /* in */) {
    return true;
}

template <class Nodes, class Id1, class Id2, class... Ids>
constexpr auto is_ordered_aux(Nodes nodes,
                              std::tuple<Id1, Id2, Ids...> /* in */) {
    if constexpr (less_than(nodes, Id2{}, Id1{})) {
        return is_ordered_aux(nodes, std::tuple<Id2, Ids...>{});
    } else {
        return false;
    }
}

} // namespace detail

template <class Nodes, class... Ids>
constexpr auto is_ordered(Nodes nodes, std::tuple<Ids...> in) {
    return detail::is_ordered_aux(nodes, in);
}

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_IS_SORTED_HPP
