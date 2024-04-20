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

#ifndef ADHOC_ADHOC_HPP
#define ADHOC_ADHOC_HPP

#include "base.hpp"
#include "constants_constexpr.hpp"

#include <array>
#include <cmath>

namespace adhoc {

template <class Input> struct exp_t : public Base<exp_t<Input>> {
    static inline auto v(double in) -> double { return std::exp(in); }
    static inline auto d(double thisv, double /* in */) -> double {
        return thisv;
    }
};

template <class Derived> auto exp(Base<Derived> /* in */) {
    return exp_t<Derived>{};
}

template <class Input> struct cos_t : public Base<cos_t<Input>> {
    static inline auto v(double in) -> double { return std::cos(in); }
    static inline auto d(double /* thisv */, double in) -> double {
        return -std::sin(in);
    }
};

template <class Derived> auto cos(Base<Derived> /* in */) {
    return cos_t<Derived>{};
}

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
};

template <class Derived> auto sqrt(Base<Derived> /* in */) {
    return sqrt_t<Derived>{};
}

template <class Input> struct log_t : public Base<log_t<Input>> {
    static inline auto v(double in) -> double { return std::log(in); }
    static inline auto d(double /* thisv */, double in) -> double {
        return 1.0 / in;
    }
};

template <class Derived> auto log(Base<Derived> /* in */) {
    return log_t<Derived>{};
}

template <class Input> struct erfc_t : public Base<erfc_t<Input>> {
    static inline auto v(double in) -> double { return std::erfc(in); }
    static inline auto d(double /* thisv */, double in) -> double {
        constexpr double two_over_root_pi =
            2.0 / constexpression::sqrt(constexpression::pi<double>());
        return -std::exp(-in * in) * two_over_root_pi;
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

} // namespace adhoc

#endif // ADHOC_ADHOC_HPP
