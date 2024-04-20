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

#ifndef ADHOC_TAPE_SIZE_FWD_HPP
#define ADHOC_TAPE_SIZE_FWD_HPP

#include "calc_order.hpp"
#include "dependency.hpp"
#include "leafs.hpp"
#include "utils.hpp"

#include <tuple>

namespace adhoc {

namespace detail {

// template <class... Roots> struct tape_size_fwd_aux_t {
//     template <class... Leafs, class... Nodes>
//     constexpr static auto call(std::tuple<Leafs...> /* in1 */,
//                                std::tuple<Nodes...> /* in2 */) {
//         return std::tuple_cat(
//             std::conditional_t<has_type<Nodes, Leafs..., Roots...>(),
//                                std::tuple<>, std::tuple<Nodes>>{}...);
//     }
// };

} // namespace detail

// template <class... Roots> constexpr auto tape_size_fwd_t(Roots... /* in */) {
//     auto constexpr leafs = detail::leafs_t<Roots...>::template call();
//     auto constexpr nodes =
//         detail::calc_order_aux_t<true, Roots...>::template call();

//     return detail::tape_size_fwd_aux_t<Roots...>::template call(leafs,
//     nodes);
// }

} // namespace adhoc

#endif // ADHOC_TAPE_SIZE_FWD_HPP
