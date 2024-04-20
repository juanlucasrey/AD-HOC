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

#ifndef ADHOC3_TAPE_BUFFER_SIZE_HPP
#define ADHOC3_TAPE_BUFFER_SIZE_HPP

#include <differential_operator/expand_tree.hpp>

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
