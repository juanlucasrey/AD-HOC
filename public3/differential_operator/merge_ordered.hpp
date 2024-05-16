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

#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_MERGE_ORDERED_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_MERGE_ORDERED_HPP

#include "../adhoc.hpp"
#include "../dependency.hpp"
#include "differential_operator.hpp"
#include "is_sorted.hpp"
#include "less_than.hpp"

namespace adhoc3 {

namespace detail {

template <class Nodes, class Id1, class... Ids1, class Id2, class... Ids2,
          class... Out>
constexpr auto merge_ordered_compare(Nodes nodes, std::tuple<Id1, Ids1...> in1,
                                     std::tuple<Id2, Ids2...> in2,
                                     std::tuple<Out...> out);

template <class Nodes, class... Ids1, class... Ids2, class... Out>
constexpr auto merge_ordered_empty(Nodes nodes, std::tuple<Ids1...> in1,
                                   std::tuple<Ids2...> in2,
                                   std::tuple<Out...> out) {

    if constexpr (sizeof...(Ids2) == 0) {
        return std::tuple_cat(out, in1);
    } else if constexpr (sizeof...(Ids1) == 0) {
        return std::tuple_cat(out, in2);
    } else {
        return merge_ordered_compare(nodes, in1, in2, out);
    }
}

template <class Nodes, class Id1, class... Ids1, class Id2, class... Ids2,
          class... Out>
constexpr auto merge_ordered_compare(Nodes nodes, std::tuple<Id1, Ids1...> in1,
                                     std::tuple<Id2, Ids2...> in2,
                                     std::tuple<Out...> /* out */) {
    if constexpr (std::is_same_v<Id1, Id2>) {
        return merge_ordered_empty(nodes, std::tuple<Ids1...>{},
                                   std::tuple<Ids2...>{},
                                   std::tuple<Out..., Id1>{});

    } else if constexpr (less_than(nodes, Id1{}, Id2{})) {
        return merge_ordered_empty(nodes, in1, std::tuple<Ids2...>{},
                                   std::tuple<Out..., Id2>{});

    } else {
        return merge_ordered_empty(nodes, std::tuple<Ids1...>{}, in2,
                                   std::tuple<Out..., Id1>{});
    }
}

} // namespace detail

template <class Nodes, class... Ids1, class... Ids2>
constexpr auto merge_ordered(Nodes nodes, std::tuple<Ids1...> in1,
                             std::tuple<Ids2...> in2) {

    static_assert(is_ordered(nodes, in1));
    static_assert(is_ordered(nodes, in2));
    return detail::merge_ordered_empty(nodes, in1, in2, std::tuple<>{});
}

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_MERGE_ORDERED_HPP
