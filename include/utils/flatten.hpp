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

#ifndef UTILS_FLATTEN_HPP
#define UTILS_FLATTEN_HPP

#include <tuple>

namespace adhoc4 {

namespace detail {

template <class TupleRes, class... IndexSequence>
constexpr auto flatten_aux3(TupleRes in1,
                            std::tuple<IndexSequence...> /* in2 */) {
    return std::tuple_cat(std::make_tuple(
        std::tuple_cat(in1, std::make_tuple(IndexSequence{})))...);
}

template <class... TupleRes, class... IndexSequence>
constexpr auto flatten_aux2(std::tuple<TupleRes...> /* res */,
                            std::tuple<IndexSequence...> in) {
    if constexpr (sizeof...(TupleRes) == 0) {
        return std::make_tuple(std::make_tuple(IndexSequence{})...);
    } else {
        return std::tuple_cat(flatten_aux3(TupleRes{}, in)...);
    }
}

template <class Res, class In> constexpr auto flatten_aux(Res res, In in) {
    if constexpr (std::tuple_size_v<In> == 0) {
        return res;
    } else {
        constexpr auto new_res = flatten_aux2(res, head(in));
        return flatten_aux(new_res, tail(in));
    }
}

} // namespace detail

template <class In> constexpr auto flatten(In idx_seqs) {
    return detail::flatten_aux(std::tuple<>{}, idx_seqs);
}

} // namespace adhoc4

#endif // UTILS_FLATTEN_HPP
