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
