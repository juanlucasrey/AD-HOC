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
