/*
 * AD-HOC, Automatic Differentiation for High Order Calculations
 *
 * This file is part of the AD-HOC distribution
 * (https://github.com/juanlucasrey/AD-HOC).
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

#ifndef ADHOC4_UTILS_INDEX_SEQUENCE_HPP
#define ADHOC4_UTILS_INDEX_SEQUENCE_HPP

#include <array>
#include <cmath>

namespace adhoc4 {

template <std::size_t Idx, std::size_t... I>
constexpr auto get(std::index_sequence<I...> /* i */) {
    return std::array<std::size_t, sizeof...(I)>{I...}[Idx];
}

} // namespace adhoc4

#endif // ADHOC4_UTILS_INDEX_SEQUENCE_HPP
