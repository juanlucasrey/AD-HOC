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

#ifndef ADHOC3_PARTITION_COMBINATIONS_HPP
#define ADHOC3_PARTITION_COMBINATIONS_HPP

#include "binomial_coefficient.hpp"

#include <cstddef>

namespace adhoc3 {

constexpr inline auto combinations(std::size_t bins, std::size_t balls) noexcept
    -> std::size_t {
    return binomial_coefficient(bins + balls - 1, bins - 1);
}

} // namespace adhoc3

#endif // ADHOC3_PARTITION_COMBINATIONS_HPP
