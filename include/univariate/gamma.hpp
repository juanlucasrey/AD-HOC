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

#ifndef UNIVARIATE_GAMMA_HPP
#define UNIVARIATE_GAMMA_HPP

#include "../base.hpp"
#include "../combinatorics/combinations.hpp"
#include "../combinatorics/factorial.hpp"
#include "../combinatorics/pow.hpp"
#include "tan.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <utility>

namespace adhoc4 {

namespace detail {

template <typename T> struct DoubleAsTemplateArg {
    T buffer{};
    constexpr DoubleAsTemplateArg(T t) : buffer(t) {}
    constexpr operator T() const { return this->buffer; }
};

template <DoubleAsTemplateArg<double> F> struct Double {
    constexpr operator double() const { return F; }
    constexpr auto static v() -> double { return F; }
};

template <class NProd, std::size_t N, std::size_t End, std::size_t Idx = 0>
constexpr auto Poly_Zeta(double x_inv_2, double b, double result) -> double {

    constexpr double convergence = std::numeric_limits<double>::epsilon() * 0.5;
    // Expansion coefficients
    // for Euler-Maclaurin summation formula
    // (2k)! / B2k
    // where B2k are Bernoulli numbers
    constexpr std::array<double, 12> BCoeffs{
        1.0 / 12.0,
        1.0 / -720.0,
        1.0 / 30240.0,
        1.0 / -1209600.0,
        1.0 / 47900160.0,
        1.0 / -1.8924375803183791606e9, /*1.307674368e12/691*/
        1.0 / 7.47242496e10,
        1.0 / -2.950130727918164224e12,  /*1.067062284288e16/3617*/
        1.0 / 1.1646782814350067249e14,  /*5.109094217170944e18/43867*/
        1.0 / -4.5979787224074726105e15, /*8.028576626982912e20/174611*/
        1.0 / 1.8152105401943546773e17,  /*1.5511210043330985984e23/854513*/
        1.0 /
            -7.1661652561756670113e18 /*1.6938241367317436694528e27/236364091*/};

    constexpr double newNProd =
        NProd{}.v() * static_cast<double>(N + 1) * static_cast<double>(N + 2);

    constexpr double coeff = newNProd * BCoeffs[Idx];
    b *= x_inv_2;

    double t = coeff * b;
    result = result + t;
    if (fabs(t / result) < convergence) {
        return result;
    }

    if constexpr ((Idx + 1) < End) {
        return Poly_Zeta<Double<newNProd>, N + 2, End, Idx + 1>(x_inv_2, b,
                                                                result);
    } else {
        return result;
    }
}

template <std::size_t N> inline auto zeta(double x) -> double {

    constexpr double convergence = std::numeric_limits<double>::epsilon() * 0.5;

    // Euler-Maclaurin summation formula
    double result = std::pow(x, -static_cast<double>(N));
    double b = 0.0;
    std::size_t i = 0;
    while ((i < 9) || (x <= 9.0)) {
        i++;
        x += 1.0;
        b = std::pow(x, -static_cast<double>(N));
        result += b;
        if (fabs(b / result) < convergence) {
            return result;
        }
    }

    constexpr double denom = 1.0 / (static_cast<double>(N) - 1.0);
    result += b * (x * denom - 0.5);

    const double x_inv = 1.0 / x;
    const double x_inv_2 = x_inv * x_inv;
    b *= x_inv;
    constexpr double t_coeff = static_cast<double>(N) * (1.0 / 12.0);
    double t = b * t_coeff;
    result = result + t;
    if (fabs(t / result) < convergence) {
        return result;
    }

    return Poly_Zeta<Double<static_cast<double>(N)>, N, 12, 1>(x_inv_2, b,
                                                               result);
}

template <std::size_t Order, std::size_t N, std::size_t Output>
inline void polygamma_aux_positive(double x, std::array<double, Output> &res) {

    constexpr double coeff =
        factorial(static_cast<unsigned int>(N)) * (N & 1 ? 1. : -1.);

    res[N] = coeff * zeta<N + 1>(x);

    if constexpr ((N + 1) < Order) {
        polygamma_aux_positive<Order, N + 1>(x, res);
    }
}

template <std::size_t Order, std::size_t N, std::size_t Output>
inline void polygamma_aux_negative(double x, std::array<double, Output> &res,
                                   std::array<double, Order> const &cotd) {
    constexpr double coeff1 = -factorial(static_cast<unsigned int>(N));
    constexpr double coeff2 = -pow<N + 1>(std::numbers::pi_v<double>);
    res[N] = coeff1 * zeta<N + 1>(1 - x) + coeff2 * cotd[N - 1];

    if constexpr ((N + 1) < Order) {
        polygamma_aux_negative<Order, N + 1>(x, res, cotd);
    }
}

template <std::size_t Order, std::size_t Output>
inline void cotd(double thisv, double /* in */,
                 std::array<double, Output> &res) {
    // we use f'(x) + (f(x))^2 + 1 = 0
    static_assert(Order <= Output);

    const double cot2 = thisv * thisv;
    tan_aux<false, false, Order>(thisv, cot2, res,
                                 std::integer_sequence<int, 0, 1>{});
}

template <std::size_t End, std::size_t Idx = 0>
inline void Poly_Aux(double &result, double in) {
    constexpr std::array<double, 7> coeffs = {
        8.33333333333333333333E-2, -2.10927960927960927961E-2,
        7.57575757575757575758E-3, -4.16666666666666666667E-3,
        3.96825396825396825397E-3, -8.33333333333333333333E-3,
        8.33333333333333333333E-2};

    result = result * in + coeffs[Idx];

    if constexpr (Idx < End) {
        Poly_Aux<End, Idx + 1>(result, in);
    }
}

template <std::size_t End> constexpr auto Poly(double in) -> double {
    double result = 0.;
    Poly_Aux<End>(result, in);
    return result;
}

inline auto digamma(double x) -> double {
    double w = 0.0;
    while (x < 10.0) {
        w += 1.0 / x;
        x += 1.0;
    }

    const double x_inv = 1.0 / x;
    const double z = x_inv * x_inv;
    return std::log(x) - (0.5 * x_inv) - z * Poly<6>(z) - w;
}

} // namespace detail

template <class Input> struct lgamma_t : public Base<lgamma_t<Input>> {
    static inline auto v(double in) -> double { return std::lgamma(in); }
    template <std::size_t Order, std::size_t Output>
    static inline void d(double /* thisv */, double in,
                         std::array<double, Output> &res) {
        static_assert(Order <= Output);

        if (in < 0) {
            double cotarg = std::numbers::pi_v<double> * std::remainder(in, 1.);
            double thisv = std::cos(cotarg) / std::sin(cotarg);

            if constexpr (Order >= 1) {
                res[0] = detail::digamma(1.0 - in) -
                         std::numbers::pi_v<double> * thisv;
            }

            if constexpr (Order >= 2) {
                std::array<double, Order> cot_d;
                cot_d.fill(0);
                detail::cotd<Order>(thisv, in, cot_d);
                detail::polygamma_aux_negative<Order, 1>(in, res, cot_d);
            }
        } else {
            if constexpr (Order >= 1) {
                res[0] = detail::digamma(in);
            }

            if constexpr (Order >= 2) {
                detail::polygamma_aux_positive<Order, 1>(in, res);
            }
        }
    }
};

template <class Derived> constexpr auto lgamma(Base<Derived> /* in */) {
    return lgamma_t<Derived>{};
}

namespace detail {

template <std::size_t Current, std::size_t CoeffIdx = 0, std::size_t Order,
          std::size_t Output>
inline void tgamma_calc(double thisv, std::array<double, Order> const &lgamma,
                        std::array<double, Output> &res) {

    if constexpr (CoeffIdx == 0) {
        res[Current] = thisv * lgamma[Current];

    } else {
        constexpr auto coeff = binomial_coefficient(Current, CoeffIdx);
        if constexpr (coeff == 1) {
            res[Current] += res[CoeffIdx - 1] * lgamma[Current - CoeffIdx];

        } else {
            res[Current] +=
                coeff * res[CoeffIdx - 1] * lgamma[Current - CoeffIdx];
        }
    }

    if constexpr (CoeffIdx < Current) {
        tgamma_calc<Current, CoeffIdx + 1>(thisv, lgamma, res);
    }
}

template <std::size_t Idx = 0, std::size_t Order, std::size_t Output>
inline void tgamma_aux(double thisv, std::array<double, Order> const &lgamma,
                       std::array<double, Output> &res) {

    tgamma_calc<Idx>(thisv, lgamma, res);

    if constexpr ((Idx + 1) < Order) {
        tgamma_aux<Idx + 1>(thisv, lgamma, res);
    }
}

} // namespace detail

template <class Input> struct tgamma_t : public Base<tgamma_t<Input>> {
    static inline auto v(double in) -> double { return std::tgamma(in); }
    template <std::size_t Order, std::size_t Output>
    static inline void d(double thisv, double in,
                         std::array<double, Output> &res) {
        static_assert(Order <= Output);

        std::array<double, Order> lgamma;
        lgamma_t<double>::d<Order>(0., in, lgamma);
        detail::tgamma_aux(thisv, lgamma, res);
    }
};

template <class Derived> constexpr auto tgamma(Base<Derived> /* in */) {
    return tgamma_t<Derived>{};
}

} // namespace adhoc4

#endif // UNIVARIATE_GAMMA_HPP
