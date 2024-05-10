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

#ifndef ADHOC4_ADHOC_HPP
#define ADHOC4_ADHOC_HPP

#include "base.hpp"
#include "constants_constexpr.hpp"
#include "univariate/tan.hpp"

#include <array>
#include <cmath>

namespace adhoc4 {

template <class Input> struct exp_t : public Base<exp_t<Input>> {
    static inline auto v(double in) -> double { return std::exp(in); }
    template <std::size_t Order, std::size_t Output>
    static inline void d(double thisv, double /* in */,
                         std::array<double, Output> &res) {
        static_assert(Order <= Output);
        // we use f'(x) - f(x) = 0
        res.fill(thisv);
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
            detail::cossin<2>(res);
        }
    }
};

template <class Derived> auto sin(Base<Derived> /* in */) {
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
            detail::cossin<2>(res);
        }
    }
};

template <class Derived> auto cos(Base<Derived> /* in */) {
    return cos_t<Derived>{};
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

    template <std::size_t Order, std::size_t Output>
    static inline void d(double thisv, double in,
                         std::array<double, Output> &res) {
        // we use 2 * x * f'(x) - f(x) = 0
        // => 2 * (x + e) * f'(x + e) - f(x + e) = 0
        // => 2 * (x + e) * sum(e^n * f(n+1)(x)/n!) - sum(e^n * f(n)(x)/n!) = 0
        // => 2*x*sum(e^n * f(n+1)(x)/n!) + 2*sum(e^n * f(n)(x)/(n-1)!) -
        // sum(e^n * f(n)(x)/n!) = 0
        // => 2*x*f(n+1)(x) + 2*n*f(n)(x) - f(n)(x) = 0
        // => f(n+1)(x) = (1/2 - n)*f(n)(x)/x = ((1 - 2n)/2)*f(n)(x)/x
        static_assert(Order <= Output);

        double const one_over_in = 1.0 / in;

        if constexpr (Order >= 1) {
            res[0] = 0.5 * thisv * one_over_in;
        }

        if constexpr (Order >= 2) {
            detail::sqrt_aux<1>(res, one_over_in);
        }
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
    template <std::size_t Order, std::size_t Output>
    static inline void d(double /* thisv */, double in,
                         std::array<double, Output> &res) {
        // we use x * f'(x) - 1 = 0
        static_assert(Order <= Output);

        if constexpr (Order >= 1) {
            res[0] = 1.0 / in;
        }

        if constexpr (Order >= 2) {
            res[1] = -res[0] * res[0];
        }

        if constexpr (Order >= 3) {
            detail::lnders<2>(res);
        }
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

    template <std::size_t Order, std::size_t Output>
    static inline void d(double /* thisv */, double in,
                         std::array<double, Output> &res) {
        // we use f''(x) + 2 * x * f'(x) = 0
        static_assert(Order <= Output);

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
    }
};

template <class Derived> auto erfc(Base<Derived> /* in */) {
    return erfc_t<Derived>{};
}

template <class Input> struct minus_t : public Base<minus_t<Input>> {
    static inline auto v(double in) -> double { return -in; }
};

template <class Derived> auto operator-(Base<Derived> /* in */) {
    return minus_t<Derived>{};
}

template <class Input1, class Input2>
struct add_t : public Base<add_t<Input1, Input2>> {
    static inline auto v(double in1, double in2) -> double { return in1 + in2; }
};

template <class Input1, class Input2>
struct sub_t : public Base<sub_t<Input1, Input2>> {
    static inline auto v(double in1, double in2) -> double { return in1 - in2; }
};

template <class Input1, class Input2>
struct mul_t : public Base<mul_t<Input1, Input2>> {
    static inline auto v(double in1, double in2) -> double { return in1 * in2; }
};

template <class Input1, class Input2>
struct div_t : public Base<div_t<Input1, Input2>> {
    static inline auto v(double in1, double in2) -> double { return in1 / in2; }
};

template <std::size_t N> struct double_t : public Base<double_t<N>> {};

} // namespace adhoc4

#endif // ADHOC4_ADHOC_HPP
