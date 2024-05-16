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
