#include <constants_double.hpp>

#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <cstring>
#include <iostream>

#include <gtest/gtest.h>

namespace adhoc2::constants {

uint64_t double_to_uint64_runtime(double x) {
    union {
        double d;
        uint64_t i;
    } y = {x};
    return y.i;
}

double uint64_to_double_runtime(uint64_t x) {
    union {
        uint64_t i;
        double d;
    } y = {x};
    return y.d;
}

uint64_t constexpr double_to_uint64(double x) {

    // std::signbit is constexpr only in C++23, so until then we have
    // __builtin_copysign
    uint64_t sign =
        __builtin_copysign(1.0, x) == 1.0 ? 0x0 : 0x8000000000000000;

    if (x == 0.) {
        return sign;
    }

    x = x < 0.0 ? -x : x;
    if (x == std::numeric_limits<double>::infinity() || x != x) {
        // std::uint64_t exp = 0x7ff0000000000000U;
        std::uint64_t exp = 0x7ffUL << 52U;

        // nan check
        if (x != x) {
            std::uint64_t val = 0x1UL << 51U;
            return sign | exp | val;
        }
        return sign | exp;
    }

    // limit between subnormal and normal
    if (x < std::numeric_limits<double>::epsilon()) {
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

    if (x < 1.) {
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
    return sign;
}

TEST(constants, double_to_uint_internet) {
    auto check = [](double temp) {
        auto tempui = double_to_uint64_runtime(temp);
        auto tempui_mine = double_to_uint64(temp);

        auto bitset1 = std::bitset<64>(tempui);
        // std::cout << bitset1 << std::endl;

        auto tempui_mine_bitset = std::bitset<64>(tempui_mine);
        // std::cout << tempui_mine_bitset << std::endl;

        EXPECT_EQ(tempui, tempui_mine);
    };

    // check2(std::numeric_limits<double>::quiet_NaN());
    // check2(std::numeric_limits<double>::infinity());
    // auto constexpr pos0 = double_to_uint64(0.);
    // auto constexpr neg0 = double_to_uint64(-0.);
    // check2(0.);
    // check2(-0.);

    // check2(2.2250738585072009e-308); // max subnormal
    // check2(4.9406564584124654e-324); // min subnormal
    // auto constexpr maxsub = double_to_uint64(2.2250738585072009e-308);
    // auto constexpr minsub = double_to_uint64(4.9406564584124654e-324);
    // check2(1.234);
    // check2(1.3407807929942597e+154);
    // check2(1.3245e+77);
    // check2(1.15792222e+89);
    // check2(1.85593e+145);
    // check2(1.999999999);
    // check2(2.0);

    check(0.5);
    check(0.25001);
    check(0.501);
    check(0.125);
    check(0.1250001);
    check(0.124999);
    check(0.25002);
    check(0.25);
    check(0.125);

    check(-0.5);
    check(-0.25001);
    check(-0.501);
    check(-0.125);
    check(-0.1250001);
    check(-0.124999);
    check(-0.25002);
    check(-0.25);
    check(-0.125);

    check(1.);
    check(2.);
    check(4.);
    check(8.);
    check(16.);
    check(0.5);
    check(0.46293);
    check(0.46293e3);
    check(0.4629334e300);
    check(-1.);
    check(-2.);
    check(-0.5);
    check(-0.46293);
    check(-0.46293e3);
    check(-0.4629334e300);

    // check(std::numeric_limits<double>::epsilon());
    // check(2.2250738585072009e-308);  // max subnormal
    // check(4.9406564584124654e-324);  // min subnormal
    // check(-2.2250738585072009e-308); // max subnormal negative
    // check(-4.9406564584124654e-324); // min subnormal negative
    // check(4.9406564584124654e-314);  // some subnormal
    // check(0);
    // check(-0);
    // check(std::numeric_limits<double>::infinity());
    // check(-std::numeric_limits<double>::infinity());
    // check(std::numeric_limits<double>::quiet_NaN());
    // check(-std::numeric_limits<double>::quiet_NaN());
    // check(std::numeric_limits<double>::max()); // 9218868437227405311
    // check(std::numeric_limits<double>::min()); // 4503599627370496
}

} // namespace adhoc2::constants
