#include <constants_double.hpp>
#include <constants_type.hpp>

#include <adhoc2.hpp>

#include <gtest/gtest.h>

namespace adhoc2::constants {

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

TEST(constants, double) {
    // constexpr double half = 0.5;
    // Double<half> temp;
    // std::cout << temp << std::endl;

    constexpr std::size_t N = convert_to_type(0.1);
    constexpr std::size_t M = convert_to_type(0.8);

    static_assert(N == 2);
    static_assert(M == 1);
}

double frexp2(double x, int *e) {
    union {
        double d;
        uint64_t i;
    } y = {x};
    int ee = y.i >> 52 & 0x7ff;

    if (!ee) {
        if (static_cast<bool>(x)) {
            x = frexp2(x * 0x1p64, e);
            *e -= 64;
        } else
            *e = 0;
        return x;
    } else if (ee == 0x7ff) {
        return x;
    }

    *e = ee - 0x3fe;
    y.i &= 0x800fffffffffffffull;
    y.i |= 0x3fe0000000000000ull;
    return y.d;
}

double constexpr mantissa_exp_to_value(double x, int e) {
    if (e == 0) {
        return x;
    } else if (e < 0) {
        return mantissa_exp_to_value(x / 2, e + 1);
    } else {
        return mantissa_exp_to_value(x * 2, e - 1);
    }
}

TEST(constants, double2) {
    double constexpr value = 0.15625;
    union {
        double d;
        uint64_t i;
    } constexpr y = {value};
    int ee = y.i >> 52 & 0x7ff;

    //  Separate value into significand in [.5, 1) and exponent.
    int exponent;
    double significand = std::frexp(value, &exponent);

    int exponent2;
    double significand2 = frexp2(value, &exponent2);
    double significand2_scaled = significand2 * 9007199254740992;

    //  Scale significand by number of digits in it, to produce an integer.
    significand = scalb(significand, std::numeric_limits<double>::digits);

    //  Adjust exponent to compensate for scaling.
    exponent -= std::numeric_limits<double>::digits;

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    std::cout << 5629499534213120.0 / 0.625 << std::endl;

    //  Set stream to print significand in full.
    std::cout << std::fixed << std::setprecision(0);

    //  Output triple with significand, base, and exponent.
    std::cout << "(" << significand << ", "
              << std::numeric_limits<double>::radix << ", " << exponent
              << ")\n";

    std::cout << std::numeric_limits<double>::digits << std::endl;
    std::cout << std::numeric_limits<double>::digits10 << std::endl;
    std::cout << std::numeric_limits<double>::max_digits10 << std::endl;
    std::cout << std::numeric_limits<double>::radix << std::endl;

    constexpr double val1 = 5629499534213120;
    constexpr int val2 = -55;
    constexpr double result = mantissa_exp_to_value(val1, val2);
    std::cout.precision(std::numeric_limits<double>::max_digits10);
    std::cout << result << std::endl;

    // constexpr double result = val1 * std::pow(2, val2);
}

TEST(constants, double3) {
    double constexpr value = 0.15625;
    constexpr double val1 = 5629499534213120;
    constexpr int val2 = -55;

    // double constexpr value = std::numeric_limits<double>::max();
    // constexpr double val1 = 9007199254740991;
    // constexpr int val2 = 971;

    int exponent;
    double significand = frexp2(value, &exponent);
    double significand_scaled = significand * 9007199254740992;
    int exponent_unbiased = exponent - std::numeric_limits<double>::digits;

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    std::cout << significand_scaled << std::endl;
    std::cout << exponent_unbiased << std::endl;
    std::int64_t final =
        static_cast<std::int64_t>(significand_scaled) +
        static_cast<std::int64_t>(exponent_unbiased) * 9007199254740992;

    // back
    double result = mantissa_exp_to_value(val1, val2);
    std::cout << result << std::endl;
    std::cout << value << std::endl;

    union myUnion {
        double dValue;
        uint64_t iValue;
    };

    myUnion myValue;
    myValue.dValue = result;
    std::cout << myValue.iValue << std::endl;

    myUnion myValue2;
    myValue2.dValue = value;
    std::cout << myValue2.iValue << std::endl;
    std::cout << final << std::endl;

    std::bitset<64> correct(myValue2.iValue);

    // static_cast<std::uint64_t>(exponent)
    // std::uint64_t final2 =
    //     static_cast<std::uint64_t>(significand_scaled) +
    //     static_cast<std::uint64_t>(exponent) * 9007199254740992;
    // std::cout << final2 << std::endl;

    // std::bitset<64> incorrect(final2);
    std::cout << correct << std::endl;
    // std::cout << incorrect << std::endl;

    union {
        double d;
        uint64_t i;
    } constexpr y = {value};
    std::uint64_t ee = y.i >> 52 & 0x7ff;
    std::bitset<64> expuint(ee);
    std::cout << expuint << std::endl;

    std::uint64_t val = (y.i << 12) >> 12;
    std::bitset<64> valuint(ee);
    std::cout << valuint << std::endl;
}

template <bool positive, unsigned int End, unsigned int Count = 1>
auto constexpr uint64_to_double_aux2(double x, std::uint32_t e,
                                     double multiplier = 2.0) -> double {
    if constexpr (Count == End) {
        return x;
    } else {
        if constexpr (positive) {
            bool apply = (e >> (Count - 1)) & 1U;
            return uint64_to_double_aux2<positive, End, Count + 1>(
                apply ? x * multiplier : x, e, multiplier * multiplier);
        } else {
            bool apply = (e >> (Count - 1)) & 1U;
            return uint64_to_double_aux2<positive, End, Count + 1>(
                apply ? x / multiplier : x, e, multiplier * multiplier);
        }
    }
}

auto constexpr uint64_to_double_aux(double x, std::int32_t e) -> double {
    // std::abs breaks constexpr!!
    // auto eu = static_cast<std::uint32_t>(std::abs(e));
    std::uint32_t eu =
        e < 0 ? static_cast<std::uint32_t>(-e) : static_cast<std::uint32_t>(e);

    if (e == 0U) {
        return x;
    } else if (e < 0) {
        // 11 is the size of the exponent. with this configuration we
        // are guaranteed 11 recursions only
        return uint64_to_double_aux2<false, 11U>(x, eu);
    } else {
        return uint64_to_double_aux2<true, 11U>(x, eu);
    }
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
        return sign ? -std::numeric_limits<double>::quiet_NaN()
                    : std::numeric_limits<double>::quiet_NaN();
    }

    // subnormal case
    if (exp == 0) {
        if (val == 0) {
            return sign ? -0 : 0;
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

template <std::uint64_t N> class Double64 {};

TEST(constants, exp_and_val_to_uint64) {
    // double constexpr value = 0.15625;
    // std::uint64_t constexpr exp2 = 1169;
    // std::uint64_t constexpr val2 = 3385410959988115;

    // double constexpr value = 0.156257953223e45;
    // std::uint64_t constexpr exp2 = 1023;
    // std::uint64_t constexpr val2 = 0;

    double constexpr value = 1.;
    union converter {
        double dValue;
        uint64_t iValue;
    };

    converter conv;
    conv.dValue = value;
    std::cout << "original double uint64" << std::endl;
    std::cout << conv.iValue << std::endl;
    std::bitset<64> valbin(conv.iValue);
    std::cout << "original double bitset" << std::endl;
    std::cout << valbin << std::endl;

    std::uint64_t exp = conv.iValue >> 52 & 0x7ff;
    std::bitset<64> expbitset(exp);
    std::cout << "exponent int" << std::endl;
    std::cout << exp << std::endl;
    std::cout << "exponent bitset" << std::endl;
    std::cout << expbitset << std::endl;

    std::uint64_t val = conv.iValue & 0xfffffffffffff;
    std::bitset<64> valbitset(val);
    std::cout << "significand int" << std::endl;
    std::cout << val << std::endl;
    std::cout << "significand bitset" << std::endl;
    std::cout << valbitset << std::endl;

    std::uint64_t exp2 = exp;
    std::uint64_t val2 = val;

    std::uint64_t valbin2 = val2 | (exp2 << 52);

    std::cout << "original double uint64 from exponent and significand"
              << std::endl;
    std::cout << valbin2 << std::endl;
    std::bitset<64> valbin2bitset(valbin2);
    std::cout << "original double bitset from exponent and significand"
              << std::endl;
    std::cout << valbin2bitset << std::endl;

    std::cout << "debug" << std::endl;
    std::cout << std::bitset<64>(valbin2) << std::endl;
    std::cout << std::bitset<64>(valbin2 >> 52U) << std::endl;
    std::cout << std::bitset<64>((valbin2 >> 52U) & 0x7ffU) << std::endl;
    double back = uint64_to_double(valbin2);

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    std::cout << "original double" << std::endl;
    std::cout << value << std::endl;
    std::cout << "original double reproduced" << std::endl;
    std::cout << back << std::endl;

    converter conv2;
    conv2.dValue = back;
    std::cout << "back double uint64" << std::endl;
    std::cout << conv2.iValue << std::endl;
    std::bitset<64> backbin(conv2.iValue);
    std::cout << "back double bitset" << std::endl;
    std::cout << backbin << std::endl;

    // int exponent2;
    // double significand2 = frexp2(value, &exponent2);
    // std::cout << exponent2 << std::endl;
    // std::cout << significand2 << std::endl;
}

TEST(constants, check_double_to_uint64) {

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
        double back = uint64_to_double(valbin2);
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
    std::uint64_t constexpr oneeq = 4607182418800017408;
    double constexpr back = uint64_to_double(oneeq);
    static_assert(back == 1.0);
    // constexpr works for normal values

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
    std::uint64_t constexpr subnormaleq = 10000000000;
    double constexpr backsub = uint64_to_double(subnormaleq);
    static_assert(backsub == 4.9406564584124654e-314);
    // constexpr works for subnormal values

    check(0);
    check(-0);
    check(std::numeric_limits<double>::infinity());
    check(-std::numeric_limits<double>::infinity());
    check(std::numeric_limits<double>::quiet_NaN());
    check(-std::numeric_limits<double>::quiet_NaN());

    // these 2 calls to uint64_to_double are the heaviest in terms of recursion
    // we modify uint64_to_double to make sure if can provide a constexpr on
    // these cases too
    check(std::numeric_limits<double>::max()); // 9218868437227405311
    check(std::numeric_limits<double>::min()); // 4503599627370496
    std::uint64_t constexpr max = 9218868437227405311;
    double constexpr backmax = uint64_to_double(max);
    static_assert(backmax == std::numeric_limits<double>::max());

    std::uint64_t constexpr min = 4503599627370496;
    double constexpr backmin = uint64_to_double(min);
    static_assert(backmin == std::numeric_limits<double>::min());
}

} // namespace adhoc2::constants
