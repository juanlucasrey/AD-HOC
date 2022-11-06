#ifndef CONSTANTS_DOUBLE_HPP
#define CONSTANTS_DOUBLE_HPP

#include <base.hpp>
#include <constants_constexpr.hpp>
#include <constants_type.hpp>

#include <array>
#include <cstdint>

namespace adhoc2::constants {

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
        double temp = x / 4.9406564584124654e-324;
        auto val = static_cast<uint64_t>(x / 4.9406564584124654e-324);
        return sign | val;
    }

    std::uint64_t exp = 0x3ff;

    if (x < 2. && x >= 1.) {
        auto val =
            static_cast<uint64_t>(x * static_cast<double>(0x10000000000000U));
        return sign | (exp << 52U) | (val & 0xfffffffffffffU);
    }

    exp = 0;
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
    if (x >= 2.) {
        if (x >= pow2[9]) {
            exp |= 0x1UL << 9U;
            x /= pow2[9];
        }

        if (x >= pow2[8]) {
            exp |= 0x1UL << 8U;
            x /= pow2[8];
        }

        if (x >= pow2[7]) {
            exp |= 0x1UL << 7U;
            x /= pow2[7];
        }

        if (x >= pow2[6]) {
            exp |= 0x1UL << 6U;
            x /= pow2[6];
        }

        if (x >= pow2[5]) {
            exp |= 0x1UL << 5U;
            x /= pow2[5];
        }

        if (x >= pow2[4]) {
            exp |= 0x1UL << 4U;
            x /= pow2[4];
        }

        if (x >= pow2[3]) {
            exp |= 0x1UL << 3U;
            x /= pow2[3];
        }

        if (x >= pow2[2]) {
            exp |= 0x1UL << 2U;
            x /= pow2[2];
        }

        if (x >= pow2[1]) {
            exp |= 0x1UL << 1U;
            x /= pow2[1];
        }

        if (x >= pow2[0]) {
            exp |= 0x1UL << 0U;
            x /= pow2[0];
        }

        exp += 0x3ff;
        auto val =
            static_cast<uint64_t>(x * static_cast<double>(0x10000000000000U));
        return sign | (exp << 52U) | (val & 0xfffffffffffffU);
    }

    // if (x < 1.) {
    if (x < (2. / pow2[9])) {
        exp |= 0x1UL << 9U;
        x *= pow2[9];
    }

    if (x < (2. / pow2[8])) {
        exp |= 0x1UL << 8U;
        x *= pow2[8];
    }

    if (x < (2. / pow2[7])) {
        exp |= 0x1UL << 7U;
        x *= pow2[7];
    }

    if (x < (2. / pow2[6])) {
        exp |= 0x1UL << 6U;
        x *= pow2[6];
    }

    if (x < (2. / pow2[5])) {
        exp |= 0x1UL << 5U;
        x *= pow2[5];
    }

    if (x < (2. / pow2[4])) {
        exp |= 0x1UL << 4U;
        x *= pow2[4];
    }

    if (x < (2. / pow2[3])) {
        exp |= 0x1UL << 3U;
        x *= pow2[3];
    }

    if (x < (2. / pow2[2])) {
        exp |= 0x1UL << 2U;
        x *= pow2[2];
    }

    if (x < (2. / pow2[1])) {
        exp |= 0x1UL << 1U;
        x *= pow2[1];
    }

    if (x < (2. / pow2[0])) {
        exp |= 0x1UL << 0U;
        x *= pow2[0];
    }

    auto exps = static_cast<std::int64_t>(exp);
    exps = 0x3ff - exps;
    exp = static_cast<std::uint64_t>(exps);
    auto val =
        static_cast<uint64_t>(x * static_cast<double>(0x10000000000000U));
    return sign | (exp << 52U) | (val & 0xfffffffffffffU);
}

} // namespace detail

template <std::uint64_t D>
struct Double : public Base<Double<D>>, public ConstBase<Double<D>> {
    constexpr auto static v() -> double { return detail::uint64_to_double(D); }
};

} // namespace adhoc2::constants

#endif // CONSTANTS_DOUBLE_HPP
