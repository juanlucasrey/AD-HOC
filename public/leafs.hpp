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

#ifndef ADHOC_LEAFS_HPP
#define ADHOC_LEAFS_HPP

#include "adhoc.hpp"
#include "constants_type.hpp"
// #include "dependency.hpp"
#include "utils.hpp"

#include <tuple>

namespace adhoc {

namespace detail {

template <typename... NodesAlive> struct leafs_t {
    template <typename... Leafs> constexpr static auto call() noexcept;
};

template <> struct leafs_t<> {
    template <typename... Leafs> constexpr static auto call() noexcept {
        return std::tuple<Leafs...>{};
    }
};

template <constants::ArgType D, typename... NodesAlive>
struct leafs_t<constants::CD<D>, NodesAlive...> {
    template <typename... Leafs> constexpr static auto call() noexcept {
        // we don't add anything because a constant is not a leaf
        return leafs_t<NodesAlive...>::template call<Leafs...>();
    }
};

template <std::size_t N, typename... NodesAlive>
struct leafs_t<double_t<N>, NodesAlive...> {
    template <typename... Leafs> constexpr static auto call() noexcept {
        using this_type = double_t<N>;

        constexpr bool already_included = has_type<this_type, Leafs...>();

        if constexpr (already_included) {
            // It's already been added, so we don't add the leaf
            return leafs_t<NodesAlive...>::template call<Leafs...>();
        } else {
            // we add the leaf
            return leafs_t<NodesAlive...>::template call<Leafs..., this_type>();
        }
    }
};

template <template <class...> class Xvariate, class... Node,
          typename... NodesAlive>
struct leafs_t<Xvariate<Node...>, NodesAlive...> {
    template <typename... Leafs> constexpr static auto call() noexcept {
        return leafs_t<Node..., NodesAlive...>::template call<Leafs...>();
    }
};

} // namespace detail

template <class... Roots> constexpr auto Leafs(Roots const &.../* o */) {
    return detail::leafs_t<Roots...>::template call();
}

} // namespace adhoc

#endif // ADHOC_LEAFS_HPP
