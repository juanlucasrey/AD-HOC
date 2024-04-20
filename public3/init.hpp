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
