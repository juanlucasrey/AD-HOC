#include <adhoc2.hpp>
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

TEST(adhoc2, derivative) {
    adouble val1(1.);
    auto valexp = exp(val1);
    auto valcos = cos(valexp);
    double derivative = 0.;
    valcos.d(derivative);
    double derivative2 = 0.;
    valexp.d(derivative2);
    derivative *= derivative2;
    // https://www.wolframalpha.com/input?i=d%2Fdx+cos%28exp%28x%29%29+%7C+x%3D1
    EXPECT_NEAR(derivative, -1.116619317445013, 1e-10);
}

} // namespace adhoc2