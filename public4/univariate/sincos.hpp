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

#ifndef ADHOC4_UNIVARIATE_SINCOS_HPP
#define ADHOC4_UNIVARIATE_SINCOS_HPP

#include "../base.hpp"

#include <array>
#include <cmath>

namespace adhoc4 {

namespace detail {

template <std::size_t N, std::size_t Order, std::size_t Output>
inline void cossin(std::array<double, Output> &res) {
    res[N] = -res[N - 2];
    if constexpr ((N + 1) < Order) {
        cossin<N + 1, Order>(res);
    }
}

} // namespace detail

template <class Input> struct sin_t : public Base<sin_t<Input>> {
    static inline auto v(double in) -> double { return std::sin(in); }
    template <std::size_t Order, std::size_t Output>
    static inline void d(double thisv, double in,
                         std::array<double, Output> &res) {
        // we use f''(x) + f(x) = 0
        static_assert(Order <= Output);

        if constexpr (Order >= 1) {
            res[0] = std::cos(in);
        }

        if constexpr (Order >= 2) {
            res[1] = -thisv;
        }

        if constexpr (Order >= 3) {
            detail::cossin<2, Order>(res);
        }
    }
};

template <class Derived> constexpr auto sin(Base<Derived> /* in */) {
    return sin_t<Derived>{};
}

template <class Input> struct cos_t : public Base<cos_t<Input>> {
    static inline auto v(double in) -> double { return std::cos(in); }
    template <std::size_t Order, std::size_t Output>
    static inline void d(double thisv, double in,
                         std::array<double, Output> &res) {
        // we use f''(x) + f(x) = 0
        static_assert(Order <= Output);

        if constexpr (Order >= 1) {
            res[0] = -std::sin(in);
        }

        if constexpr (Order >= 2) {
            res[1] = -thisv;
        }

        if constexpr (Order >= 3) {
            detail::cossin<2, Order>(res);
        }
    }
};

template <class Derived> constexpr auto cos(Base<Derived> /* in */) {
    return cos_t<Derived>{};
}

} // namespace adhoc4

#endif // ADHOC4_UNIVARIATE_SINCOS_HPP
