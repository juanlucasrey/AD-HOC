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

#ifndef ADHOC3_ADHOC_HPP
#define ADHOC3_ADHOC_HPP

#include "base.hpp"
#include "constants_constexpr.hpp"

#include <array>
#include <cmath>

namespace adhoc3 {

template <class Input> struct exp_t : public Base<exp_t<Input>> {
    static inline auto v(double in) -> double { return std::exp(in); }
    static inline auto d(double thisv, double /* in */) -> double {
        return thisv;
    }

    template <std::size_t Order>
    static inline auto d2(double thisv, double /* in */)
        -> std::array<double, Order> {
        // we use f'(x) - f(x) = 0
        std::array<double, Order> res;
        res.fill(thisv);
        return res;
    }
};

template <class Derived> auto exp(Base<Derived> /* in */) {
    return exp_t<Derived>{};
}

namespace detail {

template <std::size_t N, std::size_t Order>
inline void cossin(std::array<double, Order> &res) {
    res[N] = -res[N - 2];
    if constexpr ((N + 1) < Order) {
        cossin<N + 1>(res);
    }
}

} // namespace detail

template <class Input> struct cos_t : public Base<cos_t<Input>> {
    static inline auto v(double in) -> double { return std::cos(in); }
    static inline auto d(double /* thisv */, double in) -> double {
        return -std::sin(in);
    }

    template <std::size_t Order>
    static inline auto d2(double thisv, double in)
        -> std::array<double, Order> {
        // we use f''(x) + f(x) = 0
        std::array<double, Order> res;

        if constexpr (Order >= 1) {
            res[0] = -std::sin(in);
        }

        if constexpr (Order >= 2) {
            res[1] = -thisv;
        }

        if constexpr (Order >= 3) {
            detail::cossin<2>(res);
        }

        return res;
    }
};

template <class Derived> auto cos(Base<Derived> /* in */) {
    return cos_t<Derived>{};
}

template <class Input> struct sin_t : public Base<sin_t<Input>> {
    static inline auto v(double in) -> double { return std::sin(in); }
    static inline auto d(double /* thisv */, double in) -> double {
        return std::cos(in);
    }

    template <std::size_t Order>
    static inline auto d2(double thisv, double in)
        -> std::array<double, Order> {
        // we use f''(x) + f(x) = 0
        std::array<double, Order> res;

        if constexpr (Order >= 1) {
            res[0] = std::cos(in);
        }

        if constexpr (Order >= 2) {
            res[1] = -thisv;
        }

        if constexpr (Order >= 3) {
            detail::cossin<2>(res);
        }

        return res;
    }
};

template <class Derived> auto sin(Base<Derived> /* in */) {
    return sin_t<Derived>{};
}

namespace detail {

template <std::size_t N, std::size_t Order>
inline void sqrt_aux(std::array<double, Order> &res, double one_over_in) {
    constexpr double coeff = -static_cast<double>(N * 2 - 1) / 2.;
    res[N] = res[N - 1] * one_over_in * coeff;
    if constexpr ((N + 1) < Order) {
        sqrt_aux<N + 1>(res, one_over_in);
    }
}

} // namespace detail

template <class Input> struct sqrt_t : public Base<sqrt_t<Input>> {
    static inline auto v(double in) -> double { return std::sqrt(in); }
    // option 1
    // static inline auto d(double thisv, double in) -> double {
    //     return 0.5 * thisv / in;
    // }
    // option 2
    static inline auto d(double thisv, double /* in */) -> double {
        return 0.5 / thisv;
    }

    template <std::size_t Order>
    static inline auto d2(double thisv, double in)
        -> std::array<double, Order> {
        // we use 2 * x * f'(x) - f(x) = 0
        // => 2 * (x + e) * f'(x + e) - f(x + e) = 0
        // => 2 * (x + e) * sum(e^n * f(n+1)(x)/n!) - sum(e^n * f(n)(x)/n!) = 0
        // => 2*x*sum(e^n * f(n+1)(x)/n!) + 2*sum(e^n * f(n)(x)/(n-1)!) -
        // sum(e^n * f(n)(x)/n!) = 0
        // => 2*x*f(n+1)(x) + 2*n*f(n)(x) - f(n)(x) = 0
        // => f(n+1)(x) = (1/2 - n)*f(n)(x)/x = ((1 - 2n)/2)*f(n)(x)/x
        std::array<double, Order> res;

        double const one_over_in = 1.0 / in;

        if constexpr (Order >= 1) {
            res[0] = 0.5 * thisv * one_over_in;
        }

        if constexpr (Order >= 2) {
            detail::sqrt_aux<1>(res, one_over_in);
        }

        return res;
    }
};

template <class Derived> auto sqrt(Base<Derived> /* in */) {
    return sqrt_t<Derived>{};
}

namespace detail {

template <std::size_t N, std::size_t Order>
inline void lnders(std::array<double, Order> &res) {
    res[N] = -static_cast<double>(N) * res[N - 1] * res[0];
    if constexpr ((N + 1) < Order) {
        lnders<N + 1>(res);
    }
}

} // namespace detail

template <class Input> struct log_t : public Base<log_t<Input>> {
    static inline auto v(double in) -> double { return std::log(in); }
    static inline auto d(double /* thisv */, double in) -> double {
        return 1.0 / in;
    }

    template <std::size_t Order>
    static inline auto d2(double /* thisv */, double in)
        -> std::array<double, Order> {
        // we use x * f'(x) - 1 = 0
        std::array<double, Order> res;

        if constexpr (Order >= 1) {
            res[0] = 1.0 / in;
        }

        if constexpr (Order >= 2) {
            res[1] = -res[0] * res[0];
        }

        if constexpr (Order >= 3) {
            detail::lnders<2>(res);
        }

        return res;
    }
};

template <class Derived> auto log(Base<Derived> /* in */) {
    return log_t<Derived>{};
}

namespace detail {

template <std::size_t N, std::size_t Order>
inline void erfc_ders(std::array<double, Order> &res, double in) {
    static_assert(N > 0);
    res[N] =
        -2 * in * res[N - 1] - static_cast<double>(2 * (N - 1)) * res[N - 2];

    if constexpr ((N + 1) < Order) {
        erfc_ders<N + 1>(res, in);
    }
}

} // namespace detail

template <class Input> struct erfc_t : public Base<erfc_t<Input>> {
    static inline auto v(double in) -> double { return std::erfc(in); }
    static inline auto d(double /* thisv */, double in) -> double {
        constexpr double two_over_root_pi =
            2.0 / constexpression::sqrt(constexpression::pi<double>());
        return -std::exp(-in * in) * two_over_root_pi;
    }

    template <std::size_t Order>
    static inline auto d2(double /* thisv */, double in)
        -> std::array<double, Order> {
        // we use f''(x) + 2 * x * f'(x) = 0
        std::array<double, Order> res;

        if constexpr (Order >= 1) {
            constexpr double minus_two_over_root_pi =
                2. / constexpression::sqrt(constexpression::pi<double>());
            res[0] = -std::exp(-in * in) * minus_two_over_root_pi;
        }

        if constexpr (Order >= 2) {
            res[1] = -2.0 * in * res[0];
        }

        if constexpr (Order >= 3) {
            detail::erfc_ders<2>(res, in);
        }

        return res;
    }
};

template <class Derived> auto erfc(Base<Derived> /* in */) {
    return erfc_t<Derived>{};
}

template <class Input> struct minus_t : public Base<minus_t<Input>> {
    static inline auto v(double in) -> double { return -in; }
    static inline auto d(double /* thisv */, double /* in */) -> double {
        return -1.0;
    }
};

template <class Derived> auto operator-(Base<Derived> /* in */) {
    return minus_t<Derived>{};
}

template <class Input1, class Input2>
struct add_t : public Base<add_t<Input1, Input2>> {
    static inline auto v(double in1, double in2) -> double { return in1 + in2; }
    static inline constexpr auto d1(double /* thisv */, double /* in1 */,
                                    double /* in2 */) -> double {
        return 1.0;
    }
    static inline constexpr auto d2(double /* thisv */, double /* in1 */,
                                    double /* in2 */) -> double {
        return 1.0;
    }

    static inline constexpr auto d(double /* thisv */, double /* in1 */,
                                   double /* in2 */) {
        return std::array<double, 2>{1.0, 1.0};
    }
};

template <class Input1, class Input2>
struct sub_t : public Base<sub_t<Input1, Input2>> {
    static inline auto v(double in1, double in2) -> double { return in1 - in2; }
    static inline constexpr auto d1(double /* thisv */, double /* in1 */,
                                    double /* in2 */) -> double {
        return 1.0;
    }
    static inline constexpr auto d2(double /* thisv */, double /* in1 */,
                                    double /* in2 */) -> double {
        return -1.0;
    }
    static inline constexpr auto d(double /* thisv */, double /* in1 */,
                                   double /* in2 */) {
        return std::array<double, 2>{1.0, -1.0};
    }
};

template <class Input1, class Input2>
struct mul_t : public Base<mul_t<Input1, Input2>> {
    static inline auto v(double in1, double in2) -> double { return in1 * in2; }
    static inline auto d1(double /* thisv */, double /* in1 */, double in2)
        -> double {
        return in2;
    }
    static inline auto d2(double /* thisv */, double in1, double /* in2 */)
        -> double {
        return in1;
    }
    static inline auto d(double /* thisv */, double in1, double in2) {
        return std::array<double, 2>{in2, in1};
    }
};

template <class Input1, class Input2>
struct div_t : public Base<div_t<Input1, Input2>> {
    static inline auto v(double in1, double in2) -> double { return in1 / in2; }
    static inline auto d1(double /* thisv */, double /* in1 */, double in2)
        -> double {
        return 1.0 / in2;
    }

    // option 1 avoids a multiplication but requires value (only
    // bivariate derivative to have this requirement!)
    static inline auto d2(double thisv, double /* in1 */, double in2)
        -> double {
        return -thisv / in2;
    }

    // option 2 does not require value but is more ineficcient
    // static inline auto d2(double /* thisv */, double in1, double in2)
    //     -> double {
    //     return -in1 / (in2 * in2);
    // }

    static inline auto d(double thisv, double /* in1 */, double in2) {
        std::array<double, 2> res{1.0 / in2, 0};
        res[1] = -thisv * res[0];
        return res;
    }
};

template <std::size_t N> struct double_t : public Base<double_t<N>> {};

} // namespace adhoc3

#endif // ADHOC3_ADHOC_HPP
