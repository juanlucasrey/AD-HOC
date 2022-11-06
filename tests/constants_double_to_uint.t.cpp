#include <constants_double.hpp>

#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <cstring>
#include <iostream>

#include <gtest/gtest.h>

template <typename T> void Ref(T &&);

typedef std::array<char, sizeof(float)> float_buf;

namespace std {

// In June 16, 2019 gcc, header <bit> exists but does not have bit_cast in it.
// stolen from https://en.cppreference.com/w/cpp/numeric/bit_cast
template <class To, class From>
typename std::enable_if<
    (sizeof(To) == sizeof(From)) && std::is_trivially_copyable<From>::value &&
        std::is_trivial<To>::value && true,
    // this implementation requires that To is trivially default constructible
    To>::type constexpr bit_cast(const From &src) noexcept {
    To dst = {};
    std::memcpy(&dst, &src, sizeof(To));
    return dst;
}

// In June 16, 2019 gcc, bit_cast above is not really constexpr, so we have to
// define a custom bit_cast that works for 32-bit float *and* is constexpr.
template <>
constexpr float_buf bit_cast<float_buf, float>(const float &src) noexcept {
    // float_buf bit_cast<float_buf, float>(const float &src) noexcept {
    if (src == 0)
        return {0, 0, 0, 0};

    bool sign = src < 0.0f;
    float f = sign ? -src : src;
    uint8_t exponent = 127;
    while (f >= 2 && exponent < 255) {
        exponent += 1, f *= 0.5f;
    }
    while (f < 1 && exponent > 1)
        exponent -= 1, f *= 2.0f;
    if (f < 1)
        exponent -= 1; // denormal
    uint32_t mantissa =
        exponent == 255 ? 0 : 0x00800000 * static_cast<std::uint32_t>(f);
    uint32_t rep = (sign ? 0x80000000 : 0U) +
                   static_cast<std::uint32_t>(exponent << 23U) +
                   (mantissa & 0x7FFFFFU);
    return {(char)(rep), (char)(rep >> 8), (char)(rep >> 16),
            (char)(rep >> 24)};
}

// void bit_cast2(const double src) noexcept {
//     if (src == 0)
//         return;
//     // return {0, 0, 0, 0};

//     bool sign = src < 0.0;
//     double f = sign ? -src : src;
//     uint8_t exponent = 127;
//     while (f >= 2 && exponent < 255) {
//         exponent += 1, f *= 0.5;
//     }
//     while (f < 1 && exponent > 1)
//         exponent -= 1, f *= 2.0;
//     if (f < 1)
//         exponent -= 1; // denormal
//     uint32_t mantissa =
//         exponent == 255 ? 0 : 0x00800000 * static_cast<std::uint32_t>(f);
//     uint32_t rep = (sign ? 0x80000000 : 0U) +
//                    static_cast<std::uint32_t>(exponent << 23U) +
//                    (mantissa & 0x7FFFFFU);

//     char val1 = (char)(rep);
//     char val2 = (char)(rep >> 8);
//     char val3 = (char)(rep >> 16);
//     char val4 = (char)(rep >> 24);
//     // return {(char)(rep), (char)(rep >> 8), (char)(rep >> 16),
//     //         (char)(rep >> 24)};
//     // return {(char)(rep), (char)(rep >> 8), (char)(rep >> 16),
//     //         (char)(rep >> 24)};
// }

} // namespace std

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

uint32_t float_to_uint32_runtime(float x) {
    union {
        float d;
        uint32_t i;
    } y = {x};
    return y.i;
}

float uint32_to_float_runtime(uint32_t x) {
    union {
        uint32_t i;
        float d;
    } y = {x};
    return y.d;
}

// TEST(constants, double_to_uint) {
//     double temp = 1.0;
//     auto tempui = double_to_uint64_runtime(temp);
//     auto temp2 = uint_to_double_runtime(tempui);
//     EXPECT_EQ(temp, temp2);
// }

using namespace std;

template <size_t N1, size_t N2>
bitset<N1 + N2> concat(const bitset<N1> &b1, const bitset<N2> &b2) {
    string s1 = b1.to_string();
    string s2 = b2.to_string();
    return bitset<N1 + N2>(s1 + s2);
}

constexpr bool isPosZero(double d) {
    return d == 0.0 && __builtin_copysign(1.0, d) == 1.0;
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
        std::cout << std::bitset<64>(exp) << std::endl;
        // exp += 0x3ff;
        auto val =
            static_cast<uint64_t>(x * static_cast<double>(0x10000000000000U));
        return sign /* | (1UL << 52U) */ | (exp << 52U) |
               (val & 0xfffffffffffffU);
    }
    return sign;
}

namespace nonc {
template <bool positive, unsigned int End, unsigned int Count = 1>
auto uint64_to_double_aux2(double x, std::uint32_t e) -> double {
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

auto uint64_to_double_aux(double x, std::int32_t e) -> double {
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

auto uint64_to_double(std::uint64_t x) -> double {
    // based on
    // https://en.wikipedia.org/wiki/Double-precision_floating-point_format#Exponent_encoding

    bool sign = static_cast<bool>(x & 0x8000000000000000);
    std::uint64_t exp = (x >> 52U) & 0x7ffU;
    std::cout << std::bitset<11>(exp) << std::endl;
    std::uint64_t val = (x & 0xfffffffffffffU);

    // exceptions inf and nan
    if (exp == 2047) {
        if (val == 0) {
            return sign ? -std::numeric_limits<double>::infinity()
                        : std::numeric_limits<double>::infinity();
        }
        return sign ? -std::numeric_limits<double>::quiet_NaN()
                    : std::numeric_limits<double>::quiet_NaN();
    }

    // subnormal case
    if (exp == 0) {
        if (val == 0) {
            return sign ? -0. : 0.;
        }
        auto vald = static_cast<double>(val);
        double unsigned_result = vald * 4.9406564584124654e-324;
        return sign ? -unsigned_result : unsigned_result;
    }

    val |= 0x10000000000000U;

    std::cout << std::bitset<32>(static_cast<std::uint32_t>(exp)) << std::endl;
    auto exp_trunc = static_cast<std::int32_t>(exp) - 0x3ff;
    std::cout << std::bitset<32>(static_cast<std::uint32_t>(exp_trunc))
              << std::endl;
    auto denom = static_cast<double>(0x10000000000000U);
    auto vald = static_cast<double>(val);

    double unsigned_result = uint64_to_double_aux(vald / denom, exp_trunc);
    return sign ? -unsigned_result : unsigned_result;
}
} // namespace nonc

TEST(constants, double_to_uint_internet) {
    // auto check = [](double temp) {
    //     auto tempui = double_to_uint64_runtime(temp);
    //     auto tempui_mine = double_to_uint64(temp);

    //     auto tempuibitcheck = tempui & 0x8000000000000000;
    //     auto tempui_minebitcheck = tempui & 0x8000000000000000;
    //     EXPECT_EQ(tempuibitcheck, tempui_minebitcheck);

    //     auto bitset1 = std::bitset<64>(tempui);
    //     std::cout << bitset1 << std::endl;

    //     std::array<char, sizeof(double)> tempui2 =
    //         std::bit_cast<std::array<char, sizeof(double)>>(temp);

    //     auto print =
    //         concat(std::bitset<8>(static_cast<unsigned char>(tempui2[7])),
    //                std::bitset<8>(static_cast<unsigned char>(tempui2[6])));
    //     auto print2 = concat(
    //         print, std::bitset<8>(static_cast<unsigned char>(tempui2[5])));
    //     auto print3 = concat(
    //         print2, std::bitset<8>(static_cast<unsigned char>(tempui2[4])));
    //     auto print4 = concat(
    //         print3, std::bitset<8>(static_cast<unsigned char>(tempui2[3])));
    //     auto print5 = concat(
    //         print4, std::bitset<8>(static_cast<unsigned char>(tempui2[2])));
    //     auto print6 = concat(
    //         print5, std::bitset<8>(static_cast<unsigned char>(tempui2[1])));
    //     auto print7 = concat(
    //         print6, std::bitset<8>(static_cast<unsigned char>(tempui2[0])));

    //     EXPECT_EQ(bitset1, print7);
    // };

    auto check = [](double temp) {
        auto tempui = double_to_uint64_runtime(temp);
        auto tempui_mine = double_to_uint64(temp);

        auto bitset1 = std::bitset<64>(tempui);
        std::cout << bitset1 << std::endl;

        auto tempui_mine_bitset = std::bitset<64>(tempui_mine);
        std::cout << tempui_mine_bitset << std::endl;

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

TEST(constants, float_to_uint_internet) {
    auto check = [](float temp) {
        auto tempui = float_to_uint32_runtime(temp);
        auto bitset1 = std::bitset<32>(tempui);
        std::cout << bitset1 << std::endl;

        auto tempui2 = std::bit_cast<std::array<char, sizeof(float)>>(temp);
        // std::cout << tempui2.size() << std::endl;

        auto print =
            concat(std::bitset<8>(static_cast<unsigned char>(tempui2[3])),
                   std::bitset<8>(static_cast<unsigned char>(tempui2[2])));
        // auto print2 = concat(
        //     print, std::bitset<8>(static_cast<unsigned char>(tempui2[5])));
        // auto print3 = concat(
        //     print2, std::bitset<8>(static_cast<unsigned char>(tempui2[4])));
        // auto print4 = concat(
        //     print3, std::bitset<8>(static_cast<unsigned char>(tempui2[3])));
        // auto print5 = concat(
        //     print, std::bitset<8>(static_cast<unsigned char>(tempui2[2])));
        auto print6 = concat(
            print, std::bitset<8>(static_cast<unsigned char>(tempui2[1])));
        auto print7 = concat(
            print6, std::bitset<8>(static_cast<unsigned char>(tempui2[0])));

        std::cout << print7 << std::endl;

        EXPECT_EQ(bitset1, print7);
    };

    // check(1.f);
    // check(2.f);
    // check(4.f);
    // check(8.f);
    // check(16.f);
    // check(0.5f);

    check(0.46293f);
    // check(0.46293e3f);
    // check(-1.f);
    // check(-2.f);
    // check(-0.5f);
    // check(-0.46293f);
    // check(-0.46293e3f);
    // check(std::numeric_limits<float>::epsilon());
    // check(0.f);
    // check(-0.f);
    // check(std::numeric_limits<float>::infinity());
    // check(-std::numeric_limits<float>::infinity());
    // check(std::numeric_limits<float>::quiet_NaN());
    // check(-std::numeric_limits<float>::quiet_NaN());
    // check(std::numeric_limits<float>::max()); // 9218868437227405311
    // check(std::numeric_limits<float>::min()); // 4503599627370496
}

} // namespace adhoc2::constants
