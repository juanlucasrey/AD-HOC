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

#ifndef ADHOC3_CALC_ORDER_HPP
#define ADHOC3_CALC_ORDER_HPP

#include "adhoc.hpp"
#include "constants_type.hpp"
#include "dependency.hpp"
#include "utils.hpp"

#include <tuple>

namespace adhoc3 {

namespace detail {

template <bool IncludeInputs, bool OrderFwd, typename... NodesAlive>
struct calc_order_aux_t {
    template <typename... Operations> constexpr static auto call() noexcept;
};

template <bool IncludeInputs, bool OrderFwd>
struct calc_order_aux_t<IncludeInputs, OrderFwd> {
    template <typename... Operations> constexpr static auto call() noexcept {
        return std::tuple<Operations...>{};
    }
};

template <bool IncludeInputs, bool OrderFwd, constants::ArgType D,
          typename... NodesAlive>
struct calc_order_aux_t<IncludeInputs, OrderFwd, constants::CD<D>,
                        NodesAlive...> {
    template <typename... Operations> constexpr static auto call() noexcept {
        // we don't add anything to operations because it's a constant
        return calc_order_aux_t<IncludeInputs, OrderFwd,
                                NodesAlive...>::template call<Operations...>();
    }
};

template <bool IncludeInputs, bool OrderFwd, std::size_t N,
          typename... NodesAlive>
struct calc_order_aux_t<IncludeInputs, OrderFwd, double_t<N>, NodesAlive...> {
    template <typename... Operations> constexpr static auto call() noexcept {
        if constexpr (IncludeInputs) {
            using this_type = double_t<N>;

            constexpr bool other_types_depend_on_this =
                (equal_or_depends<NodesAlive, this_type>::call() || ...);

            if constexpr (other_types_depend_on_this) {
                // this_type will come up again because it is included on
                // NodesAlive
                return calc_order_aux_t<
                    IncludeInputs, OrderFwd,
                    NodesAlive...>::template call<Operations...>();
            } else {
                if constexpr (OrderFwd) {
                    return calc_order_aux_t<IncludeInputs, OrderFwd,
                                            NodesAlive...>::
                        template call<this_type, Operations...>();
                } else {
                    return calc_order_aux_t<IncludeInputs, OrderFwd,
                                            NodesAlive...>::
                        template call<Operations..., this_type>();
                }
            }
        } else {
            // we don't add anything to operations because it's an input
            return calc_order_aux_t<IncludeInputs, OrderFwd, NodesAlive...>::
                template call<Operations...>();
        }
    }
};

template <bool IncludeInputs, bool OrderFwd, template <class...> class Xvariate,
          class... Node, typename... NodesAlive>
struct calc_order_aux_t<IncludeInputs, OrderFwd, Xvariate<Node...>,
                        NodesAlive...> {
    template <typename... Operations> constexpr static auto call() noexcept {
        using this_type = Xvariate<Node...>;
        // it is possible that some intermediate nodes are sent as input
        // static_assert(!has_type<this_type, NodesAlive...>());

        constexpr bool other_types_depend_on_this =
            (equal_or_depends<NodesAlive, this_type>::call() || ...);

        if constexpr (other_types_depend_on_this) {
            // this_type will come up again because it is included on NodesAlive
            return calc_order_aux_t<IncludeInputs, OrderFwd, NodesAlive...>::
                template call<Operations...>();
        } else {
            if constexpr (OrderFwd) {
                return calc_order_aux_t<
                    IncludeInputs, OrderFwd, Node...,
                    NodesAlive...>::template call<this_type, Operations...>();
            } else {
                return calc_order_aux_t<
                    IncludeInputs, OrderFwd, Node...,
                    NodesAlive...>::template call<Operations..., this_type>();
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

template <bool IncludeInputs, bool OrderFwd, class... ActiveLeafsAndRoots>
constexpr auto calc_order_t() {
    if constexpr (OrderFwd) {
        // eventually, forward order should also filter active branches only
        // because those are the only values that we need to store.
        // however that will require to manage memory calculation of non-active
        // branches on an array separately. to do later.
        return detail::calc_order_aux_t<
            IncludeInputs, OrderFwd, ActiveLeafsAndRoots...>::template call();
    } else {
        auto constexpr nodes_bwd =
            detail::calc_order_aux_t<IncludeInputs, OrderFwd,
                                     ActiveLeafsAndRoots...>::template call();

        return detail::filter_active_nodes<ActiveLeafsAndRoots...>(nodes_bwd);
    }
}

} // namespace detail

template <bool IncludeInputs = false, bool OrderFwd = true,
          class... ActiveLeafsAndRoots>
constexpr auto calc_order_t(ActiveLeafsAndRoots const &.../* o */) {
    return detail::calc_order_t<IncludeInputs, OrderFwd,
                                ActiveLeafsAndRoots...>();
}

} // namespace adhoc3

#endif // ADHOC3_CALC_ORDER_HPP
