#include <adhoc2.hpp>
#include <evaluate.hpp>
#include <tape_size.hpp>
#include <utils.hpp>

#include <gtest/gtest.h>

namespace adhoc2 {

TEST(adhoc2, Initial) {
    adouble val(1.);
    auto valplus = val + adouble(1.0);
    auto const valmult = valplus * adouble(2.0);
    auto valmult2 = (valmult + adouble(2.0)) * adouble(3.0);
}

TEST(adhoc2, depends2) {
    adouble val1(1.);
    adouble val2(2.);
    constexpr auto res =
        decltype(val1)::depends3<decltype(val1), decltype(val1),
                                 decltype(val2)>();
    EXPECT_EQ(res, 2);
}

TEST(adhoc2, complexdepends2) {
    adouble val1(1.);
    adouble val2(2.);
    adouble val3(2.);
    auto valsum = val1 + val2;
    auto valprod = val1 * val2;
    auto t1 = valsum * valprod;
    auto t2 = valprod + t1;

    constexpr auto res = decltype(t2)::depends3<decltype(val1), decltype(val2),
                                                decltype(val3)>();
    EXPECT_EQ(res, 2);
}

TEST(adhoc2, complexdepends3) {
    adouble val1(1.);
    adouble val2(2.);
    adouble val3(2.);
    auto valsum = val1 + val2;
    auto valprod = val1 * val2;
    auto t1 = valsum * valprod;
    auto t2 = valprod + t1;

    constexpr auto res = decltype(t2)::depends<decltype(t1)>();
    EXPECT_EQ(res, true);
}

template <int N> constexpr auto first_n_fibs() -> std::array<double, N> {
    std::array<double, N> ret{};
    for (std::size_t i = 0; i < N; i++)
        ret[i] = 1.0;
    return ret;
}

TEST(adhoc2, complexd2) {
    adouble val1(1.);
    adouble val2(2.);
    adouble val3(2.);
    auto valsum = val1 + val2;
    auto valprod = val1 * val2;
    auto t1 = valsum * valprod;
    auto t2 = valprod + t1;

    constexpr auto a = first_n_fibs<3>();
    // constexpr auto res = t2.d2(val1, val2);
    EXPECT_EQ(a.size(), 3);
}

template <std::size_t N> class sometype {};

TEST(adhoc2, position) {
    double temp = 0;
    constexpr auto res1 = has_type<int, float, double, int, unsigned int>(temp);
    constexpr auto res2 = has_type<int, float, int, unsigned int>(temp);
    static_assert(res1);
    static_assert(!res2);
    constexpr auto res3 = idx_type<int, float, double, int, unsigned int>(temp);
    static_assert(res3 == 2);

    constexpr auto res4 = has_type2<double>();
    static_assert(!res4);
}

TEST(adhoc2, derivcomplexdexp2) {
    adouble val1(1.);
    adouble val2(2.);
    adouble val3(2.);
    auto exp2 = exp(val2);
    auto tmp1 = val1 + exp2;
    auto valsum = tmp1 + val3;
    auto val12 = val1 * val2;
    auto valprod = val12 + adouble(3);
    auto val22 = val2 * val2;
    auto valprod2 = val22 * adouble(4);
    auto valprod3 = valprod * valprod2;
    auto res = valsum * valprod3;

    constexpr std::size_t size = tape_size(res, val1, val2, val3);
    static_assert(size == 5);
    std::array<double, size> vals = {};
    vals[0] = res.d1(); // valsum
    vals[1] = res.d2(); // valprod3

    vals[2] = vals[0];
    vals[0] *= valsum.d1(); // tmp1
    vals[2] *= valsum.d2(); // val3

    vals[3] = vals[0];
    vals[0] *= tmp1.d1(); // val1
    vals[3] *= tmp1.d2(); // exp2

    vals[3] *= exp2.d(); // val2

    vals[4] = vals[1];
    vals[1] *= valprod3.d1(); // valprod
    vals[4] *= valprod3.d2(); // valprod2

    vals[1] *= valprod.d1(); // val12

    vals[0] += vals[1] * val12.d1(); // val1
    vals[3] += vals[1] * val12.d2(); // val2

    vals[4] *= valprod2.d1(); // val22

    vals[3] += vals[1] * val22.d1(); // val2
    vals[3] += vals[1] * val22.d2(); // val2

    EXPECT_NEAR(vals[0], 412.449795166, 1e-9);

    std::cout << vals[0] << ", " << vals[3] << ", " << vals[2] << std::endl;
    auto res2 = valprod3 * valsum;
    constexpr std::size_t size2 = tape_size(res2, val1, val2, val3);
    static_assert(size2 == 4);
}

TEST(adhoc2, derivcomplexdexpevaluate) {
    adouble val1(1.);
    adouble val2(2.);
    adouble val3(2.);
    auto exp2 = exp(val2);
    auto tmp1 = val1 + exp2;
    auto valsum = tmp1 + val3;
    auto val12 = val1 * val2;
    auto valprod = val12 + adouble(3);
    auto val22 = val2 * val2;
    auto valprod2 = val22 * adouble(4);
    auto valprod3 = valprod * valprod2;
    auto res = valsum * valprod3;

    auto derivatives = evaluate(res, val1, val2, val3);
}

TEST(adhoc2, derivativezero) {
    adouble val1(1.);
    adouble val2(1.);
    auto derivatives = evaluate(val1, val2);
    static_assert(derivatives.size() == 0);
}

TEST(adhoc2, derivativeexp) {
    adouble val1(1.);
    auto temp = exp(val1);
    auto derivatives = evaluate(temp, val1);
    static_assert(derivatives.size() == 1);
}

TEST(adhoc2, derivativecos) {
    adouble val1(1.);
    auto temp = cos(val1);
    auto derivatives = evaluate(temp, val1);
    static_assert(derivatives.size() == 1);
}

TEST(adhoc2, derivativeadd) {
    adouble val1(1.);
    adouble val2(1.);
    auto temp = val1 + val2;
    auto derivatives = evaluate(temp, val1, val2);
    static_assert(derivatives.size() == 2);
}

TEST(adhoc2, derivativemul) {
    adouble val1(1.);
    adouble val2(1.);
    auto temp = val1 * val2;
    auto derivatives = evaluate(temp, val1, val2);
    static_assert(derivatives.size() == 2);
}

TEST(adhoc2, derivative) {
    adouble val1(1.);
    auto valexp = exp(val1);
    auto valcos = cos(valexp);
    constexpr std::size_t size = tape_size(valcos, val1);
    static_assert(size == 1);

    double derivative = valcos.d();
    derivative *= valexp.d();
    // www.wolframalpha.com/input?i=d%2Fdx+cos%28exp%28x%29%29+%7C+x%3D1
    EXPECT_NEAR(derivative, -1.116619317445013, 1e-10);
}

TEST(adhoc2, derivativeuniv) {
    adouble val1(1.);
    auto valexp = exp(val1);
    auto valcos = cos(valexp);
    auto result = evaluate(valcos, val1);

    static_assert(result.size() == 1);
    // www.wolframalpha.com/input?i=d%2Fdx+cos%28exp%28x%29%29+%7C+x%3D1
    EXPECT_NEAR(result[0], -1.116619317445013, 1e-10);
}

} // namespace adhoc2