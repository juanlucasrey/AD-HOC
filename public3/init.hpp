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

#ifndef ADHOC3_INIT_HPP
#define ADHOC3_INIT_HPP

#include "adhoc.hpp"

#include <tuple>
#include <utility>

namespace adhoc3 {

namespace detail {

template <std::size_t... Idxs>
auto inline Init_aux(std::index_sequence<Idxs...> const & /* in */) {
    return std::tuple<double_t<Idxs>...>{};
}

} // namespace detail

template <std::size_t N> auto inline Init() {
    return detail::Init_aux(std::make_index_sequence<N>{});
}

} // namespace adhoc3

#endif // ADHOC3_INIT_HPP
