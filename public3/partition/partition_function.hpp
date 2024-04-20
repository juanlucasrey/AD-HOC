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

#ifndef ADHOC3_PARTITION_PARTITION_FUNCTION_HPP
#define ADHOC3_PARTITION_PARTITION_FUNCTION_HPP

namespace adhoc3 {

namespace detail {
constexpr auto partition_function_aux(int n, int k) -> int {
    if (k == 0) {
        return 0;
    }
    if (n == 0) {
        return 1;
    }
    if (n < 0) {
        return 0;
    }

    return partition_function_aux(n, k - 1) + partition_function_aux(n - k, k);
}
} // namespace detail

constexpr auto partition_function(int n) -> int {
    // convention
    if (n == 0) {
        return 1;
    }
    return detail::partition_function_aux(n, n);
}

} // namespace adhoc3

#endif // ADHOC3_PARTITION_PARTITION_FUNCTION_HPP
