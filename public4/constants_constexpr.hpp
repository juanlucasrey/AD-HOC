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

#ifndef ADHOC4_CONSTANTS_CONSTEXPR_HPP
#define ADHOC4_CONSTANTS_CONSTEXPR_HPP

#include <limits>

namespace adhoc4::constexpression {

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

} // namespace adhoc4::constexpression

#endif // ADHOC4_CONSTANTS_CONSTEXPR_HPP
