#include <constants_type.hpp>

#include <bitset>
#include <random>

#include <gtest/gtest.h>

namespace adhoc::constants {

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

#if __cplusplus >= 202002L
#else
TEST(Constants, uint_to_double) {

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
    static_assert(CD<4607182418800017408>::v() == 1.0);

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
    static_assert(CD<10000000000>::v() == 4.9406564584124654e-314);

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
    static_assert(CD<9218868437227405311>::v() ==
                  std::numeric_limits<double>::max());

    std::uint64_t constexpr min = 4503599627370496;
    double constexpr backmin = detail::uint64_to_double(min);
    static_assert(backmin == std::numeric_limits<double>::min());

    static_assert(CD<4503599627370496>::v() ==
                  std::numeric_limits<double>::min());
}

TEST(Constants, double_to_uint) {
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

TEST(Constants, check_encode_and_decode) {
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

TEST(Constants, check_decode_and_encode) {
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

TEST(Constants, random_uint64) {
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

TEST(Constants, static_constants) {
    // this is what we care about:
    // constexpr conversion of all kinds of doubles to and from uint64
    constexpr double one =
        detail::uint64_to_double(detail::double_to_uint64(1.0));
    static_assert(one == 1.0);

    constexpr double nan = detail::uint64_to_double(
        detail::double_to_uint64(std::numeric_limits<double>::quiet_NaN()));
    static_assert(nan != std::numeric_limits<double>::quiet_NaN());

    constexpr double zero =
        detail::uint64_to_double(detail::double_to_uint64(0.));
    static_assert(zero == 0.);

    constexpr double mzero =
        detail::uint64_to_double(detail::double_to_uint64(-0.));
    static_assert(mzero == -0.);

    constexpr double pinf = detail::uint64_to_double(
        detail::double_to_uint64(std::numeric_limits<double>::infinity()));
    static_assert(pinf == std::numeric_limits<double>::infinity());

    constexpr double minf = detail::uint64_to_double(
        detail::double_to_uint64(-std::numeric_limits<double>::infinity()));
    static_assert(minf == -std::numeric_limits<double>::infinity());

    constexpr double max = detail::uint64_to_double(
        detail::double_to_uint64(std::numeric_limits<double>::max()));
    static_assert(max == std::numeric_limits<double>::max());

    constexpr double min = detail::uint64_to_double(
        detail::double_to_uint64(std::numeric_limits<double>::min()));
    static_assert(min == std::numeric_limits<double>::min());

    constexpr double epsilon = detail::uint64_to_double(
        detail::double_to_uint64(std::numeric_limits<double>::epsilon()));
    static_assert(epsilon == std::numeric_limits<double>::epsilon());

    constexpr double maxsubnormal = detail::uint64_to_double(
        detail::double_to_uint64(2.2250738585072009e-308));
    static_assert(maxsubnormal == 2.2250738585072009e-308);

    constexpr double minsubnormal = detail::uint64_to_double(
        detail::double_to_uint64(4.9406564584124654e-324));
    static_assert(minsubnormal == 4.9406564584124654e-324);

    constexpr double maxsubnormaln = detail::uint64_to_double(
        detail::double_to_uint64(-2.2250738585072009e-308));
    static_assert(maxsubnormaln == -2.2250738585072009e-308);

    constexpr double minsubnormaln = detail::uint64_to_double(
        detail::double_to_uint64(-4.9406564584124654e-324));
    static_assert(minsubnormaln == -4.9406564584124654e-324);

    constexpr double somesubnormal = detail::uint64_to_double(
        detail::double_to_uint64(4.9406564584124654e-314));
    static_assert(somesubnormal == 4.9406564584124654e-314);

    constexpr double somesubnormaln = detail::uint64_to_double(
        detail::double_to_uint64(-4.9406564584124654e-314));
    static_assert(somesubnormaln == -4.9406564584124654e-314);

    constexpr double two =
        detail::uint64_to_double(detail::double_to_uint64(2.0));
    static_assert(two == 2.0);

    constexpr double four =
        detail::uint64_to_double(detail::double_to_uint64(4.0));
    static_assert(four == 4.0);

    constexpr double eight =
        detail::uint64_to_double(detail::double_to_uint64(8.0));
    static_assert(eight == 8.0);

    constexpr double sixteen =
        detail::uint64_to_double(detail::double_to_uint64(16.0));
    static_assert(sixteen == 16.0);

    constexpr double half =
        detail::uint64_to_double(detail::double_to_uint64(0.5));
    static_assert(half == 0.5);

    constexpr double v1 =
        detail::uint64_to_double(detail::double_to_uint64(0.46293));
    static_assert(v1 == 0.46293);

    constexpr double v2 =
        detail::uint64_to_double(detail::double_to_uint64(0.46293e3));
    static_assert(v2 == 0.46293e3);

    constexpr double v3 =
        detail::uint64_to_double(detail::double_to_uint64(0.4629334e300));
    static_assert(v3 == 0.4629334e300);

    constexpr double minusone =
        detail::uint64_to_double(detail::double_to_uint64(-1.0));
    static_assert(minusone == -1.0);

    constexpr double minustwo =
        detail::uint64_to_double(detail::double_to_uint64(-2.0));
    static_assert(minustwo == -2.0);

    constexpr double minushalf =
        detail::uint64_to_double(detail::double_to_uint64(-0.5));
    static_assert(minushalf == -0.5);

    constexpr double v1n =
        detail::uint64_to_double(detail::double_to_uint64(-0.46293));
    static_assert(v1n == -0.46293);

    constexpr double v2n =
        detail::uint64_to_double(detail::double_to_uint64(-0.46293e3));
    static_assert(v2n == -0.46293e3);

    constexpr double v3n =
        detail::uint64_to_double(detail::double_to_uint64(-0.4629334e300));
    static_assert(v3n == -0.4629334e300);
}
#endif

TEST(Constants, type_value_static_check) {
    // this is what we care about:
    // constexpr conversion of all kinds of doubles to and from uint64

    static_assert(1.0 == CD<encode(1.0)>::v());

    static_assert(std::numeric_limits<double>::quiet_NaN() !=
                  CD<encode(std::numeric_limits<double>::quiet_NaN())>::v());

    static_assert(0. == CD<encode(0.)>::v());
    static_assert(-0. == CD<encode(-0.)>::v());
    static_assert(std::numeric_limits<double>::infinity() ==
                  CD<encode(std::numeric_limits<double>::infinity())>::v());

    static_assert(-std::numeric_limits<double>::infinity() ==
                  CD<encode(-std::numeric_limits<double>::infinity())>::v());

    static_assert(std::numeric_limits<double>::max() ==
                  CD<encode(std::numeric_limits<double>::max())>::v());

    static_assert(std::numeric_limits<double>::min() ==
                  CD<encode(std::numeric_limits<double>::min())>::v());

    static_assert(std::numeric_limits<double>::epsilon() ==
                  CD<encode(std::numeric_limits<double>::epsilon())>::v());

    static_assert(2.2250738585072009e-308 ==
                  CD<encode(2.2250738585072009e-308)>::v());
    static_assert(4.9406564584124654e-324 ==
                  CD<encode(4.9406564584124654e-324)>::v());
    static_assert(-2.2250738585072009e-308 ==
                  CD<encode(-2.2250738585072009e-308)>::v());
    static_assert(-4.9406564584124654e-324 ==
                  CD<encode(-4.9406564584124654e-324)>::v());
    static_assert(4.9406564584124654e-314 ==
                  CD<encode(4.9406564584124654e-314)>::v());
    static_assert(-4.9406564584124654e-314 ==
                  CD<encode(-4.9406564584124654e-314)>::v());

    static_assert(2.0 == CD<encode(2.0)>::v());
    static_assert(4.0 == CD<encode(4.0)>::v());
    static_assert(8.0 == CD<encode(8.0)>::v());
    static_assert(16.0 == CD<encode(16.0)>::v());
    static_assert(0.5 == CD<encode(0.5)>::v());
    static_assert(0.46293 == CD<encode(0.46293)>::v());
    static_assert(0.46293e3 == CD<encode(0.46293e3)>::v());
    static_assert(0.4629334e300 == CD<encode(0.4629334e300)>::v());
    static_assert(-1.0 == CD<encode(-1.0)>::v());
    static_assert(1 == CD<encode(1)>::v());
    static_assert(-2.0 == CD<encode(-2.0)>::v());
    static_assert(-0.5 == CD<encode(-0.5)>::v());
    static_assert(-0.46293 == CD<encode(-0.46293)>::v());
    static_assert(-0.46293e3 == CD<encode(-0.46293e3)>::v());
    static_assert(-0.4629334e300 == CD<encode(-0.4629334e300)>::v());
}

#if __cplusplus >= 202002L
TEST(Constants, sizes) {
    float constexpr onef = 1.0f;
    C<float, 1.0f> onef_type;
    // different sizes, type is empty
    static_assert(sizeof(onef) == 4);
    static_assert(sizeof(onef_type) == 1);
    // same values
    static_assert(onef == onef_type);

    double constexpr oned = 1.0;
    C<double, 1.0> oned_type;
    // different sizes, type is empty
    static_assert(sizeof(oned) == 8);
    static_assert(sizeof(oned_type) == 1);
    // same values
    static_assert(oned == oned_type);

    long double constexpr oneld = 1.0l;
    C<long double, 1.0l> oneld_type;
    // different sizes, type is empty
    static_assert(sizeof(oneld) == 16);
    static_assert(sizeof(oneld_type) == 1);
    // same values
    static_assert(oneld == oneld_type);

    double constexpr one = 1.0;
    CD<one> one_type;
    CD<one> one_type2;
    // same type for same values
    static_assert(std::is_same_v<decltype(one_type), decltype(one_type)>);
    static_assert(one == one_type);
    static_assert(one == one_type2);

    // 1.000000000000000222 = std::nextafter(1)
    double constexpr oneeps = 1.000000000000000222;
    CD<oneeps> oneeps_type;
    // different types for epsilon
    static_assert(!std::is_same_v<decltype(one_type), decltype(oneeps_type)>);
    static_assert(oneeps == oneeps_type);
}
#else
TEST(Constants, sizes) {
    double constexpr oned = 1.0;
    CD<encode(1.0)> oned_type;
    // different sizes, type is empty
    static_assert(sizeof(oned) == 8);
    static_assert(sizeof(oned_type) == 1);
    // same values
    static_assert(oned == oned_type.v());

    double constexpr one = 1.0;
    CD<encode(one)> one_type;
    CD<encode(one)> one_type2;
    // same type for same values
    static_assert(std::is_same_v<decltype(one_type), decltype(one_type)>);
    static_assert(one == one_type.v());
    static_assert(one == one_type2.v());

    // 1.000000000000000222 = std::nextafter(1)
    double constexpr oneeps = 1.000000000000000222;
    CD<encode(oneeps)> oneeps_type;
    // different types for epsilon
    static_assert(!std::is_same_v<decltype(one_type), decltype(oneeps_type)>);
    static_assert(oneeps == oneeps_type.v());
}
#endif

} // namespace adhoc::constants
