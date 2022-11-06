#include <constants_double.hpp>
#include <constants_type.hpp>

#include <adhoc2.hpp>

#include <bitset>
#include <random>

#include <gtest/gtest.h>

namespace adhoc2::constants {

auto double_to_uint64_runtime(double x) -> uint64_t {
    union {
        double d;
        uint64_t i;
    } y = {x};
    return y.i;
}

auto uint64_to_double_runtime(uint64_t x) -> double {
    union {
        uint64_t i;
        double d;
    } y = {x};
    return y.d;
}

TEST(constants, frac) {
    std::cout.precision(std::numeric_limits<double>::max_digits10);
    constexpr double one_over_root_two = Frac<Const<1>, Sqrt<Const<2>>>::v();
    static_assert(one_over_root_two == 0.70710678118654757);
}

TEST(constants, fracop) {
    std::cout.precision(std::numeric_limits<double>::max_digits10);
    constexpr double one_over_root_two = Frac<Const<1>, Sqrt<Const<2>>>::v();
    static_assert(one_over_root_two == 0.70710678118654757);

    double S = 100;
    adhoc<ID> aS(S);
    auto result = aS * Frac<Const<1>, Sqrt<Const<2>>>();
}

TEST(constants, muldouble) {
    std::cout.precision(std::numeric_limits<double>::max_digits10);
    constexpr double one_over_root_two = Frac<Const<1>, Sqrt<Const<2>>>::v();
    static_assert(one_over_root_two == 0.70710678118654757);

    double S = 100;
    double result = S * Frac<Const<1>, Sqrt<Const<2>>>();
}

TEST(constants, uint_to_double) {

    auto check = [](double value) {
        union converter {
            double dValue;
            uint64_t iValue;
        };
        converter conv;
        conv.dValue = value;
        std::uint64_t exp = conv.iValue >> 52 & 0x7ff;
        std::uint64_t val = conv.iValue & 0xfffffffffffff;
        bool sign = conv.iValue >> 63;

        std::uint64_t valbin2 =
            val | (exp << 52U) | (static_cast<std::uint64_t>(sign) << 63U);
        double back = detail::uint64_to_double(valbin2);
        if (value != value) {
            // we treat nan separately
            converter conv2;
            conv2.dValue = back;
            EXPECT_EQ(conv.iValue, conv2.iValue);

        } else {
            EXPECT_EQ(value, back);
        }
    };

    check(1.);

    // constexpr works for normal values
    std::uint64_t constexpr oneeq = 4607182418800017408;
    double constexpr back = detail::uint64_to_double(oneeq);
    static_assert(back == 1.0);
    static_assert(Double<4607182418800017408>::v() == 1.0);

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
    check(std::numeric_limits<double>::epsilon());
    check(2.2250738585072009e-308);  // max subnormal
    check(4.9406564584124654e-324);  // min subnormal
    check(-2.2250738585072009e-308); // max subnormal negative
    check(-4.9406564584124654e-324); // min subnormal negative
    check(4.9406564584124654e-314);  // some subnormal

    // constexpr works for subnormal values
    std::uint64_t constexpr subnormaleq = 10000000000;
    double constexpr backsub = detail::uint64_to_double(subnormaleq);
    static_assert(backsub == 4.9406564584124654e-314);
    static_assert(Double<10000000000>::v() == 4.9406564584124654e-314);

    check(0.);
    check(-0.);
    check(std::numeric_limits<double>::infinity());
    check(-std::numeric_limits<double>::infinity());
    check(std::numeric_limits<double>::quiet_NaN());
    // check(-std::numeric_limits<double>::quiet_NaN());
    // check(std::numeric_limits<double>::signaling_NaN());
    // check(-std::numeric_limits<double>::signaling_NaN());

    // these 2 calls to uint64_to_double are the heaviest in terms of recursion
    // we modify uint64_to_double to make sure if can provide a constexpr on
    // these cases too
    check(std::numeric_limits<double>::max()); // 9218868437227405311
    check(std::numeric_limits<double>::min()); // 4503599627370496
    std::uint64_t constexpr max = 9218868437227405311;
    double constexpr backmax = detail::uint64_to_double(max);
    static_assert(backmax == std::numeric_limits<double>::max());
    static_assert(Double<9218868437227405311>::v() ==
                  std::numeric_limits<double>::max());

    std::uint64_t constexpr min = 4503599627370496;
    double constexpr backmin = detail::uint64_to_double(min);
    static_assert(backmin == std::numeric_limits<double>::min());

    static_assert(Double<4503599627370496>::v() ==
                  std::numeric_limits<double>::min());
}

TEST(constants, double_to_uint) {
    auto check = [](double temp) {
        auto tempui = double_to_uint64_runtime(temp);
        auto tempui_mine = detail::double_to_uint64(temp);

        auto bitset1 = std::bitset<64>(tempui);
        // std::cout << bitset1 << std::endl;

        auto tempui_mine_bitset = std::bitset<64>(tempui_mine);
        // std::cout << tempui_mine_bitset << std::endl;

        EXPECT_EQ(tempui, tempui_mine);
    };

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
    check(std::numeric_limits<double>::epsilon());
    check(2.2250738585072009e-308);  // max subnormal
    check(4.9406564584124654e-324);  // min subnormal
    check(-2.2250738585072009e-308); // max subnormal negative
    check(-4.9406564584124654e-324); // min subnormal negative
    check(4.9406564584124654e-314);  // some subnormal
    check(0.);
    // check(-0.);
    check(std::numeric_limits<double>::infinity());
    check(-std::numeric_limits<double>::infinity());
    check(std::numeric_limits<double>::quiet_NaN());
    // check(-std::numeric_limits<double>::quiet_NaN());
    // check(std::numeric_limits<double>::signaling_NaN());
    // check(-std::numeric_limits<double>::signaling_NaN());
    check(std::numeric_limits<double>::max()); // 9218868437227405311
    check(std::numeric_limits<double>::min()); // 4503599627370496

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

    check(std::numeric_limits<double>::epsilon());
    check(2.2250738585072009e-308);  // max subnormal
    check(4.9406564584124654e-324);  // min subnormal
    check(-2.2250738585072009e-308); // max subnormal negative
    check(-4.9406564584124654e-324); // min subnormal negative
    check(4.9406564584124654e-314);  // some subnormal
    check(0);
    // check(-0);
    check(std::numeric_limits<double>::infinity());
    check(-std::numeric_limits<double>::infinity());
    check(std::numeric_limits<double>::quiet_NaN());
    // check(-std::numeric_limits<double>::quiet_NaN());
    check(std::numeric_limits<double>::max()); // 9218868437227405311
    check(std::numeric_limits<double>::min()); // 4503599627370496
}

void check_encode_and_decode(double in) {
    // this gets rid of -0.
    // at the end, 0. and -0. have the same results on operations, so
    // there's no point separating them in the encoding
    // bool iszero = false;
    if (in == 0.) {
        // iszero = true;
        in = 0.;
    }

    // this gets rid of -quiet_NaN, +signaling_NaN, and -signaling_NaN
    // no point separating those on encoding.
    // bool isnan = false;
    if (in != in) {
        // isnan = true;
        in = std::numeric_limits<double>::quiet_NaN();
    }

    // the previous exceptions are accepted because..
    // why would anyone want to enode a NAN or a 0 on an operation?
    // adding 0 or multiplying by 0 means you need to refactor the code..
    // 0. and NaN are still supported but differentiating between 0. and -0.
    // or nans is useless
    auto in_ui = double_to_uint64_runtime(in);
    auto in_ui_to_check = detail::double_to_uint64(in);
    EXPECT_EQ(in_ui, in_ui_to_check);

    auto in_to_check = detail::uint64_to_double(in_ui_to_check);
    if (in != in) {
        EXPECT_EQ(in_to_check != in_to_check, true);
    } else {
        EXPECT_EQ(in, in_to_check);
    }
}

TEST(constants, check_encode_and_decode) {
    check_encode_and_decode(0.);
    check_encode_and_decode(-0.);
    check_encode_and_decode(std::numeric_limits<double>::quiet_NaN());
    check_encode_and_decode(-std::numeric_limits<double>::quiet_NaN());
    check_encode_and_decode(std::numeric_limits<double>::signaling_NaN());
    check_encode_and_decode(-std::numeric_limits<double>::signaling_NaN());
    check_encode_and_decode(std::numeric_limits<double>::infinity());
    check_encode_and_decode(-std::numeric_limits<double>::infinity());

    check_encode_and_decode(1.);
    check_encode_and_decode(2.);
    check_encode_and_decode(4.);
    check_encode_and_decode(8.);
    check_encode_and_decode(16.);
    check_encode_and_decode(0.5);
    check_encode_and_decode(0.46293);
    check_encode_and_decode(0.46293e3);
    check_encode_and_decode(0.4629334e300);
    check_encode_and_decode(-1.);
    check_encode_and_decode(-2.);
    check_encode_and_decode(-0.5);
    check_encode_and_decode(-0.46293);
    check_encode_and_decode(-0.46293e3);
    check_encode_and_decode(-0.4629334e300);
    check_encode_and_decode(std::numeric_limits<double>::epsilon());
    check_encode_and_decode(2.2250738585072014e-308);  // min normal
    check_encode_and_decode(2.2250738585072009e-308);  // max subnormal
    check_encode_and_decode(4.9406564584124654e-324);  // min subnormal
    check_encode_and_decode(-2.2250738585072009e-308); // max subnormal negative
    check_encode_and_decode(-4.9406564584124654e-324); // min subnormal negative
    check_encode_and_decode(4.9406564584124654e-314);  // some subnormal
    check_encode_and_decode(
        std::numeric_limits<double>::max()); // 9218868437227405311
    check_encode_and_decode(
        std::numeric_limits<double>::min()); // 4503599627370496
}

auto check_decode_and_encode(std::uint64_t in) -> bool {
    auto in_double = uint64_to_double_runtime(in);
    auto in_double_to_check = detail::uint64_to_double(in);

    if (in_double != in_double) {
        EXPECT_EQ(in_double_to_check != in_double_to_check, true);
        if ((in_double_to_check != in_double_to_check) == false) {
            return false;
        }
    } else {
        EXPECT_EQ(in_double, in_double_to_check);
        if (in_double != in_double_to_check) {
            return false;
        }
    }

    auto in_ui_to_check = detail::double_to_uint64(in_double_to_check);

    // treat 0. separately to take care of -0.
    if (in_double == 0.) {
        in = 0;
    }

    // treat nan separately, map all values to +quiet_nan
    if (in_double != in_double) {
        std::uint64_t exp = 0x7ffUL << 52U;
        std::uint64_t val = 0x1UL << 51U;
        in = exp | val;
    }
    EXPECT_EQ(in_ui_to_check, in);
    if (in_ui_to_check != in) {
        return false;
    }
    return true;
}

TEST(constants, check_decode_and_encode) {
    check_decode_and_encode(0U);                   // 0.
    check_decode_and_encode(9223372036854775808U); // -0.
    check_decode_and_encode(
        9221120237041090560U); // std::numeric_limits<double>::quiet_NaN()
    check_decode_and_encode(
        18444492273895866368U); // -std::numeric_limits<double>::quiet_NaN()
    check_decode_and_encode(
        9219994337134247936U); // std::numeric_limits<double>::signaling_NaN()
    check_decode_and_encode(
        18443366373989023744U); // -std::numeric_limits<double>::signaling_NaN()
}

TEST(constants, random_uint64) {
    unsigned int seed = 213;
    std::mt19937 gen(seed);
    std::uniform_int_distribution<std::uint64_t> dis(
        std::numeric_limits<std::uint64_t>::min(),
        std::numeric_limits<std::uint64_t>::max());

    unsigned long sims = 10000;
    for (unsigned long i = 0; i < sims; ++i) {
        auto val = dis(gen);
        bool report = check_decode_and_encode(val);
        if (report == false) {
            std::cout << "error on path: " << i << std::endl;
        }
    }
}

} // namespace adhoc2::constants
