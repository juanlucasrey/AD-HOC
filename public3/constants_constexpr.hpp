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

#ifndef ADHOC3_CONSTANTS_CONSTEXPR_HPP
#define ADHOC3_CONSTANTS_CONSTEXPR_HPP

#include <limits>

namespace adhoc3::constexpression {

namespace detail {

template <class T>
auto constexpr bailey_borwein_plouffe(T curr, T prev, T N, T power) -> T {
    return curr == prev ? curr
                        : bailey_borwein_plouffe(
                              (1. / power *
                               (4. / (8. * N + 1.) - 2. / (8. * N + 4.) -
                                1. / (8. * N + 5.) - 1. / (8. * N + 6.))) +
                                  curr,
                              curr, N + 1., power * 16.);
}

template <class T> auto constexpr sqrtNewtonRaphson(T x, T curr, T prev) -> T {
    return curr == prev ? curr
                        : sqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
}

} // namespace detail

template <class T> auto constexpr pi() -> T {
    return detail::bailey_borwein_plouffe<T>(0., 1., 0., 1.);
}

template <class T> auto constexpr sqrt(T x) -> T {
    return x >= 0 && x < std::numeric_limits<T>::infinity()
               ? detail::sqrtNewtonRaphson<T>(x, x, 0)
               : std::numeric_limits<T>::quiet_NaN();
}

} // namespace adhoc3::constexpression

#endif // ADHOC3_CONSTANTS_CONSTEXPR_HPP
