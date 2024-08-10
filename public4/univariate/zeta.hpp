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

#ifndef ADHOC4_UNIVARIATE_ZETA_HPP
#define ADHOC4_UNIVARIATE_ZETA_HPP

#include "../base.hpp"
#include "../combinatorics/combinations.hpp"
#include "../combinatorics/integer_partition_index_sequence.hpp"
#include "../combinatorics/partition_function.hpp"
#include "../utils/index_sequence.hpp"
#include "sincos.hpp"

#include <array>
#include <cmath>
#include <numbers>

#include <vector>

namespace adhoc4 {

namespace detail {

template <std::size_t n> constexpr inline auto zeta_borwein_coeffs() -> auto {
    std::array<double, n + 1> d;
    double prev = 0;
    for (std::size_t i = 0; i < (n + 1); i++) {
        auto combinations = binomial_coefficient(n + i, 2 * i);
        double temp =
            static_cast<double>(combinations) / static_cast<double>(n + i);

        double coeff2 = std::pow(4, i) * temp;
        d[i] = prev + coeff2;
        prev = d[i];
    }

    std::array<double, n> res;
    for (std::size_t i = 0; i < n; i++) {
        double const_coeff = std::pow(-1., i) * (d[i] - d[n]);
        double const_coeff2 = const_coeff / d[n];
        res[i] = const_coeff2;
    }

    return res;
}

template <std::size_t n> inline auto zeta_borwein_positive(double x) -> double {
    constexpr auto coeffs = zeta_borwein_coeffs<n>();

    double result = 0;
    for (std::size_t i = 0; i < n; i++) {
        double coeff = coeffs[i] * std::pow(i + 1, -x);
        result += coeff;
    }

    result *= -1. / (1 - std::pow(2, 1 - x));

    return result;
}

template <std::size_t n> inline auto zeta_borwein(double x) -> double {
    if (x < 0.) {
        double const oneminuss = 1. - x;

        double sinpi = std::sin(0.5 * std::numbers::pi_v<double> * x);
        double pow = std::pow(2. * std::numbers::pi_v<double>, -oneminuss);
        double tgamm = std::tgamma(oneminuss);
        double zetatemp = zeta_borwein_positive<n>(oneminuss);
        return sinpi * 2.0 * pow * tgamm * zetatemp;
    }

    return zeta_borwein_positive<n>(x);
}

template <std::size_t End, std::size_t CoeffIdx = 0, std::size_t Order,
          std::size_t... I>
void zeta4dd_aux3(std::index_sequence<I...> current, double &res,
                  std::array<double, Order> const &buff1,
                  std::array<double, Order> const &buff2) {
    constexpr auto sum_v = sum(current);
    constexpr auto coeff = BellCoeff(current);
    if constexpr (coeff == 1) {
        res += buff2[sum_v - 1] * inner_product_power(current, buff1);

    } else {
        res += buff2[sum_v - 1] * BellCoeff(current) *
               inner_product_power(current, buff1);
    }

    if constexpr ((CoeffIdx + 1) < End) {
        auto constexpr next = NextPartitionIS(current);
        zeta4dd_aux3<End, CoeffIdx + 1>(next, res, buff1, buff2);
    }
}

template <std::size_t Idx = 0, std::size_t Size>
void zeta4dd_aux2(double &res, std::array<double, Size> const &buff1,
                  std::array<double, Size> const &buff2) {
    res = 0;
    constexpr auto first = FirstPartitionIS<Idx + 1>();
    zeta4dd_aux3<partition_function(Idx + 1)>(first, res, buff1, buff2);
}

template <std::size_t Idx = 0, std::size_t Order, std::size_t Output>
void zeta4dd_aux(std::array<double, Output> &res,
                 std::array<double, Order> const &buff1,
                 std::array<double, Order> const &buff2) {
    zeta4dd_aux2<Idx>(res[Idx], buff1, buff2);
    if constexpr ((Idx + 1) < Order) {
        zeta4dd_aux<Idx + 1>(res, buff1, buff2);
    }
}

template <std::size_t End, std::size_t CoeffIdx = 0, std::size_t Order>
void mult_aux(double &res, double rightval,
              std::array<double, Order> const &right, double leftval,
              std::array<double, Order> const &left) {
    if constexpr (CoeffIdx == 0) {
        res += rightval * left[End - 1];
    } else if constexpr (CoeffIdx == End) {
        res += leftval * right[End - 1];
    } else {
        constexpr auto binomial = binomial_coefficient(End, CoeffIdx);
        res += binomial * left[CoeffIdx - 1] * right[End - 1 - CoeffIdx];
    }
    if constexpr (CoeffIdx < End) {
        mult_aux<End, CoeffIdx + 1>(res, rightval, right, leftval, left);
    }
}

template <std::size_t Idx = 0, std::size_t Order, std::size_t Output>
void mult(std::array<double, Output> &res, double rightval,
          std::array<double, Order> const &right, double leftval,
          std::array<double, Order> const &left) {
    static_assert(Order <= Output);
    res[Idx] = 0;
    mult_aux<Idx + 1>(res[Idx], rightval, right, leftval, left);
    if constexpr ((Idx + 1) < Order) {
        mult<Idx + 1>(res, rightval, right, leftval, left);
    }
}

template <std::size_t Order, std::size_t n, std::size_t Output>
inline auto zeta_borwein_positive_der(double x, std::array<double, Output> &res)
    -> double {

    constexpr auto coeffs = zeta_borwein_coeffs<n>();
    std::array<double, Order> res1;
    res1.fill(0);
    double result = 0;
    for (std::size_t i = 0; i < n; i++) {
        double coeff = coeffs[i] * std::pow(i + 1, -x);
        result += coeff;

        if (i != 0) {
            double temp = coeff;
            for (std::size_t j = 0; j < Order; j++) {
                temp *= -std::log(i + 1);
                res1[j] += temp;
            }
        }
    }

    std::array<double, Order> buffer;
    double const powval = std::pow(2.0, (1.0 - x));

    buffer[0] = powval * std::numbers::ln2_v<double>;
    for (std::size_t k = 1; k < Order; k++) {
        buffer[k] = buffer[k - 1] * -std::numbers::ln2_v<double>;
    }

    double const frac = 1.0 / (1.0 - powval);

    std::array<double, Order> buffer2;
    inv_t<double>::d<Order>(frac, frac, buffer2);

    std::array<double, Order> res2;
    zeta4dd_aux(res2, buffer, buffer2);

    for (std::size_t j = 0; j < Order; j++) {
        res2[j] *= -1.;
    }
    double result2 = -1. / (1 - std::pow(2, 1 - x));

    std::array<double, Order> resfinal;
    mult(resfinal, result2, res2, result, res1);
    result *= result2;

    for (std::size_t j = 0; j < Order; j++) {
        res[j] = resfinal[j];
    }

    return result;
}

template <std::size_t Order, std::size_t n, std::size_t Output>
inline void zeta_borwein_der(double x, std::array<double, Output> &res) {
    if (x < 0.) {
        double const oneminuss = 1. - x;

        constexpr double coeff = 0.5 * std::numbers::pi_v<double>;
        double const val = coeff * std::remainder(x, 4);
        double const sinpi = std::sin(val);

        std::array<double, Order> sinders;
        sin_t<double>::d<Order>(sinpi, val, sinders);
        double multcoeff = coeff;
        for (std::size_t j = 0; j < Order; j++) {
            sinders[j] *= multcoeff;
            multcoeff *= coeff;
        }

        double const pow =
            std::pow(2. * std::numbers::pi_v<double>, -oneminuss);

        std::array<double, Order> powders;
#ifndef __clang__
        // in GCC, std::log can be called as a constexpr function!
        constexpr double ln2pi = std::log(2.0 * std::numbers::pi_v<double>);
#else
        constexpr double ln2pi = 1.8378770664093453;
#endif

        double multpow = ln2pi * pow;
        for (std::size_t j = 0; j < Order; j++) {
            powders[j] = multpow;
            multpow *= ln2pi;
        }

        std::array<double, Order> gammaders;
        double const tgamma_val = std::tgamma(oneminuss);
        tgamma_t<double>::d<Order>(tgamma_val, oneminuss, gammaders);
        for (std::size_t j = 0; j < Order; j++) {
            gammaders[j] *= std::pow(-1., j + 1);
        }

        std::array<double, Order> zetaders;
        double const zetatemp =
            zeta_borwein_positive_der<Order, n>(oneminuss, zetaders);
        for (std::size_t j = 0; j < Order; j++) {
            zetaders[j] *= std::pow(-1., j + 1);
        }

        std::array<double, Order> temp1;
        double val1 = tgamma_val * zetatemp;
        mult(temp1, zetatemp, zetaders, tgamma_val, gammaders);

        std::array<double, Order> temp2;
        double val2 = pow * val1;
        mult(temp2, val1, temp1, pow, powders);

        std::array<double, Order> temp3;
        double val3 = sinpi * val2;
        mult(temp3, val2, temp2, sinpi, sinders);

        for (std::size_t j = 0; j < Order; j++) {
            temp3[j] *= 2.0;
        }

        for (std::size_t j = 0; j < Order; j++) {
            res[j] = temp3[j];
        }
    } else {
        zeta_borwein_positive_der<Order, n>(x, res);
    }
}

} // namespace detail

template <class Input>
struct riemann_zeta_t : public Base<riemann_zeta_t<Input>> {
    static inline auto v(double in) -> double { return std::riemann_zeta(in); }
    template <std::size_t Order, std::size_t Output>
    static inline void d(double /* thisv */, double in,
                         std::array<double, Output> &res) {
        static_assert(Order <= Output);

        detail::zeta_borwein_der<Order, 20>(in, res);
    }
};

template <class Derived> auto riemann_zeta(Base<Derived> /* in */) {
    return riemann_zeta_t<Derived>{};
}

} // namespace adhoc4

#endif // ADHOC4_UNIVARIATE_ZETA_HPP
