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

#ifndef ADHOC_CALC_ORDER_HPP
#define ADHOC_CALC_ORDER_HPP

#include "adhoc.hpp"
#include "constants_type.hpp"
#include "dependency.hpp"
#include "utils.hpp"

#include <tuple>

namespace adhoc {

namespace detail {

template <bool OrderFwd, typename... NodesAlive> struct calc_order_aux_t {
    template <typename... Operations> constexpr static auto call() noexcept;
};

template <bool OrderFwd> struct calc_order_aux_t<OrderFwd> {
    template <typename... Operations> constexpr static auto call() noexcept {
        return std::tuple<Operations...>{};
    }
};

template <bool OrderFwd, constants::ArgType D, typename... NodesAlive>
struct calc_order_aux_t<OrderFwd, constants::CD<D>, NodesAlive...> {
    template <typename... Operations> constexpr static auto call() noexcept {
        // we don't add anything to operations because it's a constant
        return calc_order_aux_t<OrderFwd,
                                NodesAlive...>::template call<Operations...>();
    }
};

template <bool OrderFwd, std::size_t N, typename... NodesAlive>
struct calc_order_aux_t<OrderFwd, double_t<N>, NodesAlive...> {
    template <typename... Operations> constexpr static auto call() noexcept {
        // we don't add anything to operations because it's an input
        return calc_order_aux_t<OrderFwd,
                                NodesAlive...>::template call<Operations...>();
    }
};

template <bool OrderFwd, template <class...> class Xvariate, class... Node,
          typename... NodesAlive>
struct calc_order_aux_t<OrderFwd, Xvariate<Node...>, NodesAlive...> {
    template <typename... Operations> constexpr static auto call() noexcept {
        using this_type = Xvariate<Node...>;
        // it is possible that some intermediate nodes are sent as input
        // static_assert(!has_type<this_type, NodesAlive...>());

        constexpr bool other_types_depend_on_this =
            (equal_or_depends<NodesAlive, this_type>::call() || ...);

        if constexpr (other_types_depend_on_this) {
            // this_type will come up again because it is included on NodesAlive
            return calc_order_aux_t<OrderFwd, NodesAlive...>::template call<
                Operations...>();
        } else {
            if constexpr (OrderFwd) {
                return calc_order_aux_t<OrderFwd, Node..., NodesAlive...>::
                    template call<this_type, Operations...>();
            } else {
                return calc_order_aux_t<OrderFwd, Node..., NodesAlive...>::
                    template call<Operations..., this_type>();
            }
        }
    }
};

template <class... ActiveLeafsAndRoots, class... Nodes>
constexpr auto filter_active_nodes(std::tuple<Nodes...> /* in */) {
    return std::tuple_cat(
        std::conditional_t<(depends_many<Nodes, ActiveLeafsAndRoots...>()),
                           std::tuple<Nodes>, std::tuple<>>{}...);
}

template <bool OrderFwd, class... ActiveLeafsAndRoots>
constexpr auto calc_order_t() {
    if constexpr (OrderFwd) {
        // eventually, forward order should also filter active branches only
        // because those are the only values that we need to store.
        // however that will require to manage memory calculation of non-active
        // branches on an array separately. to do later.
        return detail::calc_order_aux_t<
            OrderFwd, ActiveLeafsAndRoots...>::template call();
    } else {
        auto constexpr nodes_bwd =
            detail::calc_order_aux_t<OrderFwd,
                                     ActiveLeafsAndRoots...>::template call();

        return detail::filter_active_nodes<ActiveLeafsAndRoots...>(nodes_bwd);
    }
}

} // namespace detail

template <bool OrderFwd = true, class... ActiveLeafsAndRoots>
constexpr auto calc_order_t(ActiveLeafsAndRoots const &.../* o */) {
    return detail::calc_order_t<OrderFwd, ActiveLeafsAndRoots...>();
}

} // namespace adhoc

#endif // ADHOC_CALC_ORDER_HPP
