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

#ifndef ADHOC_CONSTANTS_TYPE_HPP
#define ADHOC_CONSTANTS_TYPE_HPP

#include "base.hpp"

#include <limits>

#if !(__cplusplus >= 202002L)
#include <array>
#include <cstdint>
#endif

namespace adhoc::constants {

template <class Derived> struct ConstBase {
    operator double() const { return Derived::v(); }
};

#if __cplusplus >= 202002L
namespace detail {
template <typename T> struct AsTemplateArg {
    T buffer{};
    constexpr AsTemplateArg(T t) : buffer(t) {}
    constexpr operator T() const { return this->buffer; }
};
} // namespace detail

template <typename T, detail::AsTemplateArg<T> F>
struct C : public Base<C<T, F>> {
    constexpr operator T() const { return F; }
    constexpr auto static v() -> T { return F; }
};

template <detail::AsTemplateArg<float> F> using CF = C<float, F>;
template <detail::AsTemplateArg<double> F> using CD = C<double, F>;
template <detail::AsTemplateArg<long double> F> using CLD = C<long double, F>;

constexpr auto encode(double x) -> double {
    // do nothing, for backward compatibility
    return x;
}

using ArgType = constants::detail::AsTemplateArg<double>;

#else

namespace detail {

template <bool positive, unsigned int End, unsigned int Count = 1>
auto constexpr uint64_to_double_aux2(double x, std::uint32_t e) -> double {
    if constexpr (Count == End) {
        return x;
    } else {
        bool apply = (e >> (Count - 1)) & 1U;

        constexpr std::array<double, 10> pow2{
            2.,                     // 2^0b1
            4.,                     // 2^0b10
            16.,                    // 2^0b100
            256.,                   // 2^0b1000
            65536.,                 // 2^0b10000
            4294967296.,            // 2^0b100000
            1.8446744073709552e+19, // 2^0b1000000
            3.4028236692093846e+38, // 2^0b10000000
            1.157920892373162e+77,  // 2^0b10000000
            1.3407807929942597e+154 // 2^0b100000000
        };

        if constexpr (positive) {
            x = apply ? x * pow2[Count - 1] : x;
        } else {
            x = apply ? x / pow2[Count - 1] : x;
        }

        return uint64_to_double_aux2<positive, End, Count + 1>(x, e);
    }
}

auto constexpr uint64_to_double_aux(double x, std::int32_t e) -> double {
    // std::abs breaks constexpr!!
    // auto eu = static_cast<std::uint32_t>(std::abs(e));
    std::uint32_t eu =
        e < 0 ? static_cast<std::uint32_t>(-e) : static_cast<std::uint32_t>(e);

    if (e == 0U) {
        return x;
    }

    if (e < 0) {
        // 11 is the size of the exponent. with this configuration we
        // are guaranteed 11 recursions only
        return uint64_to_double_aux2<false, 11U>(x, eu);
    }

    return uint64_to_double_aux2<true, 11U>(x, eu);
}

auto constexpr uint64_to_double(std::uint64_t x) -> double {
    // based on
    // https://en.wikipedia.org/wiki/Double-precision_floating-point_format#Exponent_encoding

    bool sign = static_cast<bool>(x & 0x8000000000000000);
    std::uint64_t exp = (x >> 52U) & 0x7ffU;
    std::uint64_t val = (x & 0xfffffffffffffU);

    // exceptions inf and nan
    if (exp == 2047) {
        if (val == 0) {
            return sign ? -std::numeric_limits<double>::infinity()
                        : std::numeric_limits<double>::infinity();
        }

        // we don't care about different types of nan
        return std::numeric_limits<double>::quiet_NaN();

        // bool msb = static_cast<bool>(val & 0x8000000000000U);
        // if (msb) {
        //     return sign ? -std::numeric_limits<double>::quiet_NaN()
        //                 : std::numeric_limits<double>::quiet_NaN();
        // }
        // return sign ? -std::numeric_limits<double>::signaling_NaN()
        //             : std::numeric_limits<double>::signaling_NaN();
    }

    // subnormal case
    if (exp == 0) {
        if (val == 0) {
            // return sign ? -0. : 0.;
            return 0.; // deliberate choice to not differentiate + and -
        }
        auto vald = static_cast<double>(val);
        double unsigned_result = vald * 4.9406564584124654e-324;
        return sign ? -unsigned_result : unsigned_result;
    }

    val |= 0x10000000000000U;

    auto exp_trunc = static_cast<std::int32_t>(exp) - 0x3ff;
    auto denom = static_cast<double>(0x10000000000000U);
    auto vald = static_cast<double>(val);

    double unsigned_result = uint64_to_double_aux(vald / denom, exp_trunc);
    return sign ? -unsigned_result : unsigned_result;
}

template <bool positive, unsigned int End, unsigned int Count = End>
void constexpr double_to_uint64_aux(double &x, std::uint64_t &exp) {
    if constexpr (Count == 0) {
        return;
    } else {
        // bool apply = (e >> (Count - 1)) & 1U;

        constexpr std::array<double, 10> pow2{
            2.,                     // 2^0b1
            4.,                     // 2^0b10
            16.,                    // 2^0b100
            256.,                   // 2^0b1000
            65536.,                 // 2^0b10000
            4294967296.,            // 2^0b100000
            1.8446744073709552e+19, // 2^0b1000000
            3.4028236692093846e+38, // 2^0b10000000
            1.157920892373162e+77,  // 2^0b10000000
            1.3407807929942597e+154 // 2^0b100000000
        };

        if constexpr (positive) {
            if (x >= pow2[Count - 1]) {
                exp |= 0x1UL << (Count - 1);
                x /= pow2[Count - 1];
            }
        } else {
            if (x < (2. / pow2[Count - 1])) {
                exp |= 0x1UL << (Count - 1);
                x *= pow2[Count - 1];
            }
        }

        double_to_uint64_aux<positive, End, Count - 1>(x, exp);
    }
}

constexpr auto double_to_uint64(double x) -> std::uint64_t {
    // we treat 0. and -0. the same. we could have tried to treat separately
    // with
    // __builtin_copysign(1.0, x) == 1.0 ?
    // instead of
    // x >= 0 ?
    // but really, what for? 0. and -0. will produce the same results.
    // besides, they should never be used as a constant.
    // besides, __builtin_copysign is not portable.
    // to treat 0. and -0. separately, we could use std::signbit in the future,
    // since it will be constexpr in C++23 but by c++20 we won't need this
    // conversion, so, it's useless.
    if (x == 0.) {
        return 0;
    }

    // we treat all nans the same. again, no point differentiating between nans
    // for this encoding usecase.
    if (x != x) {
        std::uint64_t exp = 0x7ffUL << 52U;
        std::uint64_t val = 0x1UL << 51U;
        return exp | val;
    }

    // after 0 and nan have been treated, the sign is meaningful
    uint64_t sign = x >= 0 ? 0x0 : 0x8000000000000000;
    x = x < 0.0 ? -x : x;

    if (x == std::numeric_limits<double>::infinity()) {
        std::uint64_t exp = 0x7ffUL << 52U;
        return sign | exp;
    }

    // limit between subnormal and normal
    if (x < std::numeric_limits<double>::min()) {
        auto val = static_cast<uint64_t>(x / 4.9406564584124654e-324);
        return sign | val;
    }

    std::uint64_t exp = 0;
    if (x >= 1.) {
        double_to_uint64_aux<true, 10>(x, exp);
        exp += 0x3ff;
    } else {
        double_to_uint64_aux<false, 10>(x, exp);
        auto exps = static_cast<std::int64_t>(exp);
        exps = 0x3ff - exps;
        exp = static_cast<std::uint64_t>(exps);
    }

    auto val =
        static_cast<uint64_t>(x * static_cast<double>(0x10000000000000U));
    return sign | (exp << 52U) | (val & 0xfffffffffffffU);
}

} // namespace detail

constexpr auto encode(double x) -> std::uint64_t {
    return detail::double_to_uint64(x);
}

template <std::uint64_t D>
struct CD : public Base<CD<D>>, public ConstBase<CD<D>> {
    constexpr auto static v() -> double { return detail::uint64_to_double(D); }
};

using ArgType = std::uint64_t;
#endif

} // namespace adhoc::constants

#endif // ADHOC_CONSTANTS_TYPE_HPP
