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

#ifndef ADHOC4_SORT_HPP
#define ADHOC4_SORT_HPP

#include "differential_operator.hpp"
#include "utils/tuple.hpp"

#include <cstddef>
#include <type_traits>

namespace adhoc4 {

namespace detail {

template <class Op1, class Op2, class Nodes>
constexpr auto less_than_check_empty(Op1 /* in1 */, Op2 /* in2 */,
                                     Nodes /* nodes */) -> bool;

template <class Id1, std::size_t Order1, class Id2, std::size_t Order2,
          class Id1Remaining, class Id2Remaining, class Nodes>
constexpr auto less_than_single(der::d<Order1, Id1> /* in1 */,
                                der::d<Order2, Id2> /* in2 */,
                                Id1Remaining id1rem, Id2Remaining id2rem,
                                Nodes nodes) -> bool {
    constexpr auto idx1 = find<Nodes, Id1>();
    constexpr auto idx2 = find<Nodes, Id2>();
    if constexpr (idx1 < idx2) {
        return false;
    } else if constexpr (idx1 > idx2) {
        return true;
    } else if constexpr (Order1 < Order2) {
        return false;
    } else if constexpr (Order1 > Order2) {
        return true;
    } else {
        static_assert(std::is_same_v<der::d<Order1, Id1>, der::d<Order2, Id2>>);
        return less_than_check_empty(id1rem, id2rem, nodes);
    }
}

// template <class Op11, class... Ops1, class Op21, class... Ops2, class Nodes>
// constexpr auto less_than_check_first(std::tuple<Op11, Ops1...> in1,
//                                      std::tuple<Op21, Ops2...> in2,
//                                      Nodes nodes) {
//     return less_than_single(Op11{}, Op21{}, std::tuple<Ops1...>{},
//                             std::tuple<Ops2...>{}, nodes);
// }

template <class In1, class In2, class Nodes>
constexpr auto less_than_check_first(In1 in1, In2 in2, Nodes nodes) -> bool {
    return less_than_single(head(in1), head(in2), tail(in1), tail(in2), nodes);
}

template <class Op1, class Op2, class Nodes>
constexpr auto less_than_check_empty(Op1 in1, Op2 in2, Nodes nodes) -> bool {
    if constexpr (std::is_same_v<const Op2, const std::tuple<>>) {
        return false;
    } else if constexpr (std::is_same_v<const Op1, const std::tuple<>>) {
        return true;
    } else {
        // we can now assume that both ids are different and non empty
        return less_than_check_first(in1, in2, nodes);
    }
}

template <class Op1, class Op2, class Nodes>
constexpr auto less_than(Op1 in1, Op2 in2, Nodes nodes) -> bool {
    if constexpr (std::is_same_v<Op1, Op2>) {
        return false;
    } else {
        return less_than_check_empty(in1, in2, nodes);
    }
}

template <class ResultOp, class CandidateOp, class Nodes>
constexpr auto return_largest(const ResultOp res, const CandidateOp candidate,
                              Nodes nodes) {
    if constexpr (less_than(candidate, res, nodes)) {
        return res;
    } else {
        return candidate;
    }
}

template <class ResultOp, class Nodes>
constexpr auto find_largest_aux(ResultOp res, std::tuple<> /* diff_operators */,
                                Nodes /* nodes */) {
    return res;
}

template <class ResultOp, class Op1, class... Ops, class Nodes>
constexpr auto find_largest_aux(ResultOp res,
                                std::tuple<Op1, Ops...> /* diff_operators */,
                                Nodes nodes) {
    return find_largest_aux(return_largest(res, Op1{}, nodes),
                            std::tuple<Ops...>{}, nodes);
}

template <class Op1, class... Ops, class Nodes>
constexpr auto find_largest(std::tuple<Op1, Ops...> /* diff_operators */,
                            Nodes nodes) {
    return find_largest_aux(Op1{}, std::tuple<Ops...>{}, nodes);
}

template <class Result, class Nodes>
constexpr auto
sort_differential_operators_aux(std::tuple<> /* diff_operators */, Result res,
                                Nodes /* nodes */) {
    return res;
}

template <class Ops, class Result, class Nodes>
constexpr auto sort_differential_operators_aux(Ops diff_operators, Result res,
                                               Nodes nodes) {
    constexpr auto largest = find_largest(diff_operators, nodes);
    constexpr auto rem = remove(diff_operators, largest);
    return sort_differential_operators_aux(
        rem, std::tuple_cat(res, std::make_tuple(largest)), nodes);
}

} // namespace detail

template <class Ops, class Nodes>
constexpr auto sort_differential_operators(Ops diff_operators, Nodes nodes) {
    // if constexpr (is_ordered(nodes, diff_operators)) {
    //     return diff_operators;
    // } else {
    return detail::sort_differential_operators_aux(diff_operators,
                                                   std::tuple<>{}, nodes);
    // }
}

} // namespace adhoc4

#endif // ADHOC4_SORT_HPP
