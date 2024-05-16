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

#ifndef ADHOC3_TAPE_BUFFER_SIZE_HPP
#define ADHOC3_TAPE_BUFFER_SIZE_HPP

#include "differential_operator/expand_tree.hpp"

#include <cstddef>
#include <tuple>

namespace adhoc3 {

namespace detail {

template <std::size_t MaxWidth, class Nodes, class DersIn, class DersActive>
constexpr auto
tape_buffer_size_aux(Nodes /* nodes */, std::tuple<> /* expandedtree */,
                     DersIn /* dersin */, DersActive /* dersactive */)
    -> std::size_t {
    return MaxWidth;
}

template <class DerivativeNodes, class DersIn, class DersActive,
          class CandidateNode>
constexpr auto conditional_node(DerivativeNodes dernodes, DersIn dersin,
                                DersActive dersactive, CandidateNode candnode) {

    if constexpr (contains(dernodes, candnode) &&
                  !contains(dersactive, candnode) &&
                  !contains(dersin, candnode)) {
        return std::make_tuple(candnode);
    } else {
        return std::tuple<>{};
    }
}

template <std::size_t MaxWidth, class Nodes, class FirstNode,
          class... RestNodes, class DersIn, class DersActive>
constexpr auto
tape_buffer_size_aux(Nodes nodes,
                     std::tuple<FirstNode, RestNodes...> /* expandedtree */,
                     DersIn dersin, DersActive dersactive) -> std::size_t {
    constexpr auto dersactive_current_removed = remove(dersactive, FirstNode{});
    constexpr auto expansion_result = expand_single(nodes, FirstNode{});

    constexpr auto result_filtered = std::apply(
        [dersin, dersactive_current_removed](auto... single_output) {
            return std::tuple_cat(
                conditional_node(std::tuple<RestNodes...>{}, dersin,
                                 dersactive_current_removed, single_output)...);
        },
        expansion_result);

    constexpr auto new_dersactive =
        std::tuple_cat(dersactive_current_removed, result_filtered);

    constexpr std::size_t new_current_widt =
        std::tuple_size_v<decltype(new_dersactive)>;

    constexpr std::size_t new_max_widt =
        new_current_widt > MaxWidth ? new_current_widt : MaxWidth;

    return tape_buffer_size_aux<new_max_widt>(nodes, std::tuple<RestNodes...>{},
                                              dersin, new_dersactive);
}

} // namespace detail

template <class Nodes, class ExpandedTree, class DersIn>
constexpr auto tape_buffer_size(Nodes nodes, ExpandedTree expandedtree,
                                DersIn dersin) -> std::size_t {
    return detail::tape_buffer_size_aux<0U>(nodes, expandedtree, dersin,
                                            std::tuple<>{});
}

} // namespace adhoc3

#endif // ADHOC3_TAPE_BUFFER_SIZE_HPP
