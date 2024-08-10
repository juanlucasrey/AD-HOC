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

#ifndef ADHOC4_ADHOC_HPP
#define ADHOC4_ADHOC_HPP

#include "base.hpp"
#include "univariate/gamma.hpp"
#include "univariate/sincos.hpp"
#include "univariate/tan.hpp"

#ifndef __clang__
#include "univariate/zeta.hpp"
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>

namespace adhoc4 {

template <class Input> struct exp_t : public Base<exp_t<Input>> {
    static inline auto v(double in) -> double { return std::exp(in); }
    template <std::size_t Order, std::size_t Output>
    static inline void d(double thisv, double /* in */,
                         std::array<double, Output> &res) {
        static_assert(Order <= Output);
        // we use f'(x) - f(x) = 0
        std::fill_n(res.begin(), Order, thisv);
    }
};

template <class Derived> auto exp(Base<Derived> /* in */) {
    return exp_t<Derived>{};
}

namespace detail {

template <std::size_t N, std::size_t Order, std::size_t Output>
inline void asin_aux(std::array<double, Output> &res, double in,
                     double denominator) {
    auto constexpr coeff1 = static_cast<double>(2 * (N - 1) + 1);
    auto constexpr coeff2 = static_cast<double>((N - 1) * (N - 1));
    res[N] = (coeff1 * in * res[N - 1] + coeff2 * res[N - 2]) * denominator;
    if constexpr ((N + 1) < Order) {
        asin_aux<N + 1, Order>(res, in, denominator);
    }
}

} // namespace detail

template <class Input> struct asin_t : public Base<asin_t<Input>> {
    static inline auto v(double in) -> double { return std::asin(in); }
    template <std::size_t Order, std::size_t Output>
    static inline void d(double /* thisv */, double in,
                         std::array<double, Output> &res) {
        // we use (x^2 - 1)f''(x) + x * f'(x) = 0
        static_assert(Order <= Output);

        const double denominator = 1. / (1. - in * in);
        if constexpr (Order >= 1) {
            res[0] = std::sqrt(denominator);
        }

        if constexpr (Order >= 2) {
            res[1] = in * res[0] * denominator;
        }

        if constexpr (Order >= 3) {
            res[2] = (3. * in * res[1] + res[0]) * denominator;
        }

        if constexpr (Order >= 4) {
            detail::asin_aux<3, Order>(res, in, denominator);
        }
    }
};

template <class Derived> auto asin(Base<Derived> /* in */) {
    return asin_t<Derived>{};
}

template <class Input> struct acos_t : public Base<acos_t<Input>> {
    static inline auto v(double in) -> double { return std::acos(in); }
    template <std::size_t Order, std::size_t Output>
    static inline void d(double /* thisv */, double in,
                         std::array<double, Output> &res) {
        // we use (x^2 - 1)f''(x) + x * f'(x) = 0
        static_assert(Order <= Output);

        const double denominator = 1. / (1. - in * in);
        if constexpr (Order >= 1) {
            res[0] = -std::sqrt(denominator);
        }

        if constexpr (Order >= 2) {
            res[1] = in * res[0] * denominator;
        }

        if constexpr (Order >= 3) {
            res[2] = (3. * in * res[1] + res[0]) * denominator;
        }

        if constexpr (Order >= 4) {
            detail::asin_aux<3, Order>(res, in, denominator);
        }
    }
};

template <class Derived> auto acos(Base<Derived> /* in */) {
    return acos_t<Derived>{};
}

namespace detail {

template <std::size_t N, std::size_t Order, std::size_t Output>
inline void atan_aux(std::array<double, Output> &res, double in,
                     double denominator) {
    auto constexpr coeff1 = static_cast<double>(2 * N);
    auto constexpr coeff2 = static_cast<double>(N * (N - 1));
    res[N] = -(coeff1 * in * res[N - 1] + coeff2 * res[N - 2]) * denominator;
    if constexpr ((N + 1) < Order) {
        atan_aux<N + 1, Order>(res, in, denominator);
    }
}

} // namespace detail

template <class Input> struct atan_t : public Base<atan_t<Input>> {
    static inline auto v(double in) -> double { return std::atan(in); }
    template <std::size_t Order, std::size_t Output>
    static inline void d(double /* thisv */, double in,
                         std::array<double, Output> &res) {
        // we use (x^2 - 1)f''(x) + x * f'(x) = 0
        static_assert(Order <= Output);

        const double denominator = 1. / (1. + in * in);
        if constexpr (Order >= 1) {
            res[0] = denominator;
        }

        if constexpr (Order >= 2) {
            res[1] = -2. * in * res[0] * denominator;
        }

        if constexpr (Order >= 3) {
            detail::atan_aux<2, Order>(res, in, denominator);
        }
    }
};

template <class Derived> auto atan(Base<Derived> /* in */) {
    return atan_t<Derived>{};
}

namespace detail {

template <std::size_t N, std::size_t Order, std::size_t Output>
inline void cossinh(std::array<double, Output> &res) {
    res[N] = res[N - 2];
    if constexpr ((N + 1) < Order) {
        cossinh<N + 1, Order>(res);
    }
}

} // namespace detail

template <class Input> struct sinh_t : public Base<sinh_t<Input>> {
    static inline auto v(double in) -> double { return std::sinh(in); }
    template <std::size_t Order, std::size_t Output>
    static inline void d(double thisv, double in,
                         std::array<double, Output> &res) {
        // we use f''(x) + f(x) = 0
        static_assert(Order <= Output);

        if constexpr (Order >= 1) {
            res[0] = std::cosh(in);
        }

        if constexpr (Order >= 2) {
            res[1] = thisv;
        }

        if constexpr (Order >= 3) {
            detail::cossinh<2, Order>(res);
        }
    }
};

template <class Derived> auto sinh(Base<Derived> /* in */) {
    return sinh_t<Derived>{};
}

template <class Input> struct cosh_t : public Base<cosh_t<Input>> {
    static inline auto v(double in) -> double { return std::cosh(in); }
    template <std::size_t Order, std::size_t Output>
    static inline void d(double thisv, double in,
                         std::array<double, Output> &res) {
        // we use f''(x) + f(x) = 0
        static_assert(Order <= Output);

        if constexpr (Order >= 1) {
            res[0] = std::sinh(in);
        }

        if constexpr (Order >= 2) {
            res[1] = thisv;
        }

        if constexpr (Order >= 3) {
            detail::cossinh<2, Order>(res);
        }
    }
};

template <class Derived> auto cosh(Base<Derived> /* in */) {
    return cosh_t<Derived>{};
}

namespace detail {

template <std::size_t N, std::size_t Order, std::size_t Output>
inline void asinh_aux(std::array<double, Output> &res, double in,
                      double denominator) {
    auto constexpr coeff1 = static_cast<double>(2 * (N - 1) + 1);
    auto constexpr coeff2 = static_cast<double>((N - 1) * (N - 1));
    res[N] = -(coeff1 * in * res[N - 1] + coeff2 * res[N - 2]) * denominator;
    if constexpr ((N + 1) < Order) {
        asinh_aux<N + 1, Order>(res, in, denominator);
    }
}

} // namespace detail

template <class Input> struct asinh_t : public Base<asinh_t<Input>> {
    static inline auto v(double in) -> double { return std::asinh(in); }
    template <std::size_t Order, std::size_t Output>
    static inline void d(double /* thisv */, double in,
                         std::array<double, Output> &res) {
        // we use (x^2 - 1)f''(x) + x * f'(x) = 0
        static_assert(Order <= Output);

        const double denominator = 1. / (1. + in * in);
        if constexpr (Order >= 1) {
            res[0] = std::sqrt(denominator);
        }

        if constexpr (Order >= 2) {
            res[1] = -in * res[0] * denominator;
        }

        if constexpr (Order >= 3) {
            res[2] = -(3. * in * res[1] + res[0]) * denominator;
        }

        if constexpr (Order >= 4) {
            detail::asinh_aux<3, Order>(res, in, denominator);
        }
    }
};

template <class Derived> auto asinh(Base<Derived> /* in */) {
    return asinh_t<Derived>{};
}

template <class Input> struct acosh_t : public Base<acosh_t<Input>> {
    static inline auto v(double in) -> double { return std::acosh(in); }
    template <std::size_t Order, std::size_t Output>
    static inline void d(double /* thisv */, double in,
                         std::array<double, Output> &res) {
        // we use (x^2 - 1)f''(x) + x * f'(x) = 0
        static_assert(Order <= Output);

        const double denominator = 1. / (in * in - 1.);
        if constexpr (Order >= 1) {
            res[0] = std::sqrt(denominator);
        }

        if constexpr (Order >= 2) {
            res[1] = -in * res[0] * denominator;
        }

        if constexpr (Order >= 3) {
            res[2] = -(3. * in * res[1] + res[0]) * denominator;
        }

        if constexpr (Order >= 4) {
            detail::asinh_aux<3, Order>(res, in, denominator);
        }
    }
};

template <class Derived> auto acosh(Base<Derived> /* in */) {
    return acosh_t<Derived>{};
}

template <class Input> struct atanh_t : public Base<atanh_t<Input>> {
    static inline auto v(double in) -> double { return std::atanh(in); }
    template <std::size_t Order, std::size_t Output>
    static inline void d(double /* thisv */, double in,
                         std::array<double, Output> &res) {
        // we use (x^2 - 1)f''(x) + x * f'(x) = 0
        static_assert(Order <= Output);

        const double denominator = 1. / (in * in - 1.);
        if constexpr (Order >= 1) {
            res[0] = -denominator;
        }

        if constexpr (Order >= 2) {
            res[1] = -2. * in * res[0] * denominator;
        }

        if constexpr (Order >= 3) {
            detail::atan_aux<2, Order>(res, in, denominator);
        }
    }
};

template <class Derived> auto atanh(Base<Derived> /* in */) {
    return atanh_t<Derived>{};
}

namespace detail {

template <std::size_t N, std::size_t Order, std::size_t Output>
inline void sqrt_aux(std::array<double, Output> &res, double one_over_in) {
    constexpr double coeff = -static_cast<double>(N * 2 - 1) / 2.;
    res[N] = res[N - 1] * one_over_in * coeff;
    if constexpr ((N + 1) < Order) {
        sqrt_aux<N + 1, Order>(res, one_over_in);
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
        // => 2 * (x + e) * sum(e^n * f(n+1)(x)/n!) - sum(e^n * f(n)(x)/n!)
        // = 0
        // => 2*x*sum(e^n * f(n+1)(x)/n!) + 2*sum(e^n * f(n)(x)/(n-1)!) -
        // sum(e^n * f(n)(x)/n!) = 0
        // => 2*x*f(n+1)(x) + 2*n*f(n)(x) - f(n)(x) = 0
        // => f(n+1)(x) = (1/2 - n)*f(n)(x)/x = ((1 - 2n)/2)*f(n)(x)/x
        static_assert(Order <= Output);

        double const one_over_in = 1. / in;

        if constexpr (Order >= 1) {
            res[0] = 0.5 * thisv * one_over_in;
        }

        if constexpr (Order >= 2) {
            detail::sqrt_aux<1, Order>(res, one_over_in);
        }
    }
};

template <class Derived> auto sqrt(Base<Derived> /* in */) {
    return sqrt_t<Derived>{};
}

namespace detail {

template <std::size_t N, std::size_t Order, std::size_t Output>
inline void lnders(std::array<double, Output> &res) {
    res[N] = -static_cast<double>(N) * res[N - 1] * res[0];
    if constexpr ((N + 1) < Order) {
        lnders<N + 1, Order>(res);
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
            res[0] = 1. / in;
        }

        if constexpr (Order >= 2) {
            res[1] = -res[0] * res[0];
        }

        if constexpr (Order >= 3) {
            detail::lnders<2, Order>(res);
        }
    }
};

template <class Derived> auto log(Base<Derived> /* in */) {
    return log_t<Derived>{};
}

namespace detail {

template <std::size_t N, std::size_t Order, std::size_t Output>
inline void erfc_ders(std::array<double, Output> &res, double in) {
    static_assert(N > 0);
    res[N] =
        -2. * in * res[N - 1] - static_cast<double>(2 * (N - 1)) * res[N - 2];

    if constexpr ((N + 1) < Order) {
        erfc_ders<N + 1, Order>(res, in);
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
                2. * std::numbers::inv_sqrtpi_v<double>;
            res[0] = -std::exp(-in * in) * minus_two_over_root_pi;
        }

        if constexpr (Order >= 2) {
            res[1] = -2. * in * res[0];
        }

        if constexpr (Order >= 3) {
            detail::erfc_ders<2, Order>(res, in);
        }
    }
};

template <class Derived> auto erfc(Base<Derived> /* in */) {
    return erfc_t<Derived>{};
}

#ifndef __clang__
namespace detail {

template <std::size_t N, std::size_t Order, std::size_t Output>
inline void comp_ellint_1_ders(std::array<double, Output> &res, double in,
                               double coeff, double denominator) {
    static_assert(N > 0);
    auto constexpr coeff1 = static_cast<double>(N);
    auto constexpr coeff2 =
        static_cast<double>(3 * (N - 1) * (N - 1) + 3 * (N - 1) + 1);
    auto constexpr coeff3 = static_cast<double>((N - 1) * (N - 1) * (N - 1));
    res[N] = (coeff3 * res[N - 3] + coeff2 * in * res[N - 2] -
              coeff1 * coeff * res[N - 1]) *
             denominator;

    if constexpr ((N + 1) < Order) {
        comp_ellint_1_ders<N + 1, Order>(res, in, coeff, denominator);
    }
}

} // namespace detail

template <class Input>
struct comp_ellint_1_t : public Base<comp_ellint_1_t<Input>> {
    static inline auto v(double in) -> double { return std::comp_ellint_1(in); }

    template <std::size_t Order, std::size_t Output>
    static inline void d(double thisv, double in,
                         std::array<double, Output> &res) {
        // we use x * (1 - x^2) f''(x) + (1 - 3 * x^2) * f'(x) - x * f(x) =
        // 0
        static_assert(Order <= Output);

        const double in_sq = in * in;
        const double temp1 = (1. - in_sq);
        const double denominator = 1. / (in * temp1);
        if constexpr (Order >= 1) {
            res[0] = (std::comp_ellint_2(in) - temp1 * thisv) * denominator;
        }

        if constexpr (Order >= 2) {
            const double coeff = (1. - 3 * in_sq);
            res[1] = (in * thisv - coeff * res[0]) * denominator;
            if constexpr (Order >= 3) {
                res[2] = (thisv + 7. * in * res[0] - 2. * coeff * res[1]) *
                         denominator;
            }

            if constexpr (Order >= 4) {
                detail::comp_ellint_1_ders<3, Order>(res, in, coeff,
                                                     denominator);
            }
        }
    }
};

template <class Derived> auto comp_ellint_1(Base<Derived> /* in */) {
    return comp_ellint_1_t<Derived>{};
}

namespace detail {

template <std::size_t N, std::size_t Order, std::size_t Output>
inline void comp_ellint_2_ders(std::array<double, Output> &res, double in,
                               double in_sq, double denominator) {
    static_assert(N > 0);
    auto constexpr coeff1 = static_cast<double>((N - 1) * (N - 1) * (N - 1) -
                                                2 * (N - 1) * (N - 1));
    auto constexpr coeff2 = static_cast<double>(3 * (N - 1) * (N - 1) - N);
    auto constexpr coeff3 = static_cast<double>(3 * N - 2);
    res[N] = (coeff1 * res[N - 3] + coeff2 * in * res[N - 2] -
              (static_cast<double>(N) - coeff3 * in_sq) * res[N - 1]) *
             denominator;

    if constexpr ((N + 1) < Order) {
        comp_ellint_2_ders<N + 1, Order>(res, in, in_sq, denominator);
    }
}

} // namespace detail

template <class Input>
struct comp_ellint_2_t : public Base<comp_ellint_2_t<Input>> {
    static inline auto v(double in) -> double { return std::comp_ellint_2(in); }

    template <std::size_t Order, std::size_t Output>
    static inline void d(double thisv, double in,
                         std::array<double, Output> &res) {
        // we use x * (1 - x^2) f''(x) + (1 - x^2) * f'(x) + x * f(x) = 0
        static_assert(Order <= Output);

        const double in_sq = in * in;
        const double temp1 = (1. - in_sq);
        const double denominator = 1. / (in * temp1);
        if constexpr (Order >= 1) {
            res[0] = (thisv - std::comp_ellint_1(in)) / in;
        }

        if constexpr (Order >= 2) {
            res[1] = -(in * thisv + temp1 * res[0]) * denominator;
            if constexpr (Order >= 3) {
                res[2] = (-thisv + in * res[0] - (2. - 4. * in_sq) * res[1]) *
                         denominator;
            }

            if constexpr (Order >= 4) {
                detail::comp_ellint_2_ders<3, Order>(res, in, in_sq,
                                                     denominator);
            }
        }
    }
};

template <class Derived> auto comp_ellint_2(Base<Derived> /* in */) {
    return comp_ellint_2_t<Derived>{};
}
#endif

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

namespace detail {
template <size_t N> struct StringLiteral {
    constexpr StringLiteral(const char (&str)[N]) {
        std::copy_n(str, N, value);
    }

    char value[N]{};
};

} // namespace detail

template <detail::StringLiteral literal>
struct double_t : public Base<double_t<literal>> {};

#define ADHOC(x) double_t<#x> x

} // namespace adhoc4

#endif // ADHOC4_ADHOC_HPP
