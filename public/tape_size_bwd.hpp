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

#ifndef ADHOC_TAPE_SIZE_BWD_HPP
#define ADHOC_TAPE_SIZE_BWD_HPP

#include "calc_order.hpp"
#include "dependency.hpp"

#include <algorithm>
#include <tuple>

namespace adhoc {

namespace detail {

template <std::size_t CurrentWidth, std::size_t MaxWidth,
          class... ActiveLeafsAndRoots, class... StoredNodes>
constexpr auto tape_size_bwd_aux(std::tuple<ActiveLeafsAndRoots...>,
                                 std::tuple<>, std::tuple<StoredNodes...>)
    -> std::size_t {
    return MaxWidth;
}

template <std::size_t CurrentWidth, std::size_t MaxWidth,
          template <class, class> class Bivariate, class Node1, class Node2,
          class... ActiveLeafsAndRoots, class... NodesAlive,
          class... StoredNodes>
constexpr auto
    tape_size_bwd_aux(std::tuple<ActiveLeafsAndRoots...>,
                      std::tuple<Bivariate<Node1, Node2>, NodesAlive...>,
                      std::tuple<StoredNodes...>) -> std::size_t;

template <std::size_t CurrentWidth, std::size_t MaxWidth,
          template <class> class Univariate, class Node,
          class... ActiveLeafsAndRoots, class... NodesAlive,
          class... StoredNodes>
constexpr auto tape_size_bwd_aux(std::tuple<ActiveLeafsAndRoots...>,
                                 std::tuple<Univariate<Node>, NodesAlive...>,
                                 std::tuple<StoredNodes...>) -> std::size_t {
    using this_type = Univariate<Node>;

    constexpr bool is_this_node_root =
        has_type<this_type, ActiveLeafsAndRoots...>();

    constexpr bool is_next_node_new =
        !has_type<Node, StoredNodes..., ActiveLeafsAndRoots...>();

    constexpr std::size_t new_current_widt =
        CurrentWidth - !is_this_node_root + is_next_node_new;

    constexpr std::size_t new_max_widt =
        new_current_widt > MaxWidth ? new_current_widt : MaxWidth;

    return tape_size_bwd_aux<new_current_widt, new_max_widt>(
        std::tuple<ActiveLeafsAndRoots...>{}, std::tuple<NodesAlive...>{},
        std::tuple<Node, StoredNodes...>{});
}

template <std::size_t CurrentWidth, std::size_t MaxWidth,
          template <class, class> class Bivariate, class Node1, class Node2,
          class... ActiveLeafsAndRoots, class... NodesAlive,
          class... StoredNodes>
constexpr auto
tape_size_bwd_aux(std::tuple<ActiveLeafsAndRoots...>,
                  std::tuple<Bivariate<Node1, Node2>, NodesAlive...>,
                  std::tuple<StoredNodes...>) -> std::size_t {
    using this_type = Bivariate<Node1, Node2>;

    constexpr bool is_this_node_root =
        has_type<this_type, ActiveLeafsAndRoots...>();

    constexpr bool is_next_node1_new =
        !has_type<Node1, StoredNodes..., ActiveLeafsAndRoots...>() &&
        equal_or_depends_many<Node1, ActiveLeafsAndRoots...>();

    constexpr bool is_next_node2_new =
        !has_type<Node2, StoredNodes..., ActiveLeafsAndRoots...>() &&
        equal_or_depends_many<Node2, ActiveLeafsAndRoots...>();

    constexpr bool is_same_node = std::is_same_v<Node1, Node2>;

    constexpr std::size_t new_current_widt =
        CurrentWidth - !is_this_node_root + is_next_node1_new +
        (!is_same_node && is_next_node2_new);

    constexpr std::size_t new_max_widt =
        new_current_widt > MaxWidth ? new_current_widt : MaxWidth;

    return tape_size_bwd_aux<new_current_widt, new_max_widt>(
        std::tuple<ActiveLeafsAndRoots...>{}, std::tuple<NodesAlive...>{},
        std::tuple<Node1, Node2, StoredNodes...>{});
}

} // namespace detail

template <class... ActiveLeafsAndRoots>
constexpr auto tape_size_bwd() -> std::size_t {
    auto constexpr nodes =
        detail::calc_order_t<false, ActiveLeafsAndRoots...>();
    return detail::tape_size_bwd_aux<0U, 0U>(
        std::tuple<ActiveLeafsAndRoots...>{}, nodes, std::tuple<>{});
}

} // namespace adhoc

#endif // ADHOC_TAPE_SIZE_BWD_HPP
