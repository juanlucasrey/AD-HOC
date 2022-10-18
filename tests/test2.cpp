#include <adhoc2.hpp>

#include <gtest/gtest.h>

namespace adhoc2 {

TEST(adhoc2, Initial) {
    adouble val(1.);
    auto valplus = val + 1.0;
    auto const valmult = valplus * 2.0;
    auto valmult2 = (valmult + 2.0) * 3.0;
    // double result = valmult2;
    // double deriv = d(result) / d(val);
    // double deriv = d(result, val);
    // EXPECT_EQ(result, 18.);
    // EXPECT_EQ(deriv, 0.);
}

TEST(adhoc2, doperator) {
    adouble val1(1.);
    adouble val2(2.);
    auto valplus = val1 + 1.0;
    double res1 = val1.d(val1);
    double res2 = val1.d(val2);
    double res3 = val1.d(valplus);
    EXPECT_EQ(res1, 1.);
    EXPECT_EQ(res2, 0.);
    EXPECT_EQ(res3, 0.);
}

TEST(adhoc2, doperatorplus) {
    adouble val1(1.);
    adouble val2(2.);
    auto valplus = val1 + 1.0;
    double res1 = valplus.d(val1);
    double res2 = valplus.d(val2);
    EXPECT_EQ(res1, 1.);
    EXPECT_EQ(res2, 0.);
}

TEST(adhoc2, doperatormult) {
    adouble val1(1.);
    adouble val2(2.);
    auto valmul = val1 * 2.0;
    double res1 = valmul.d(val1);
    double res2 = valmul.d(val2);
    EXPECT_EQ(res1, 2.);
    EXPECT_EQ(res2, 0.);
}

TEST(adhoc2, activeadd) {
    adouble val1(1.);
    adouble val2(2.);
    auto valsum = (val1 * 3.0) + val2;

    double res1 = valsum.d(val1);
    double res2 = valsum.d(val2);
    EXPECT_EQ(res1, 3.);
    EXPECT_EQ(res2, 1.);

    auto valsum2 = (val1 * 3.0) + val1;

    double res1_2 = valsum2.d(val1);
    double res2_2 = valsum2.d(val2);
    EXPECT_EQ(res1_2, 4.);
    EXPECT_EQ(res2_2, 0.);
}

TEST(adhoc2, activemul) {
    adouble val1(1.);
    adouble val2(2.);
    auto valsum = val1 * val2;

    double res1 = valsum.d(val1);
    double res2 = valsum.d(val2);
    EXPECT_EQ(res1, 2.);
    EXPECT_EQ(res2, 1.);
}

TEST(adhoc2, complex) {
    adouble val1(1.);
    adouble val2(2.);
    auto valsum = val1 + val2;
    auto valprod = val1 * val2;
    auto t1 = valsum * valprod;
    auto t2 = valprod + t1;

    double res1 = t2.d(val1);
    double res2 = t2.d(val2);
    EXPECT_EQ(res1, 10.);
    EXPECT_EQ(res2, 6.);
}

TEST(adhoc2, depends2) {
    adouble val1(1.);
    adouble val2(2.);
    constexpr auto res =
        decltype(val1)::depends2<decltype(val1), decltype(val1),
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

    constexpr auto res = decltype(t2)::depends2<decltype(val1), decltype(val2),
                                                decltype(val3)>();
    EXPECT_EQ(res, 2);
}

template <int N> constexpr std::array<double, N> first_n_fibs() {
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

TEST(adhoc2, Initial2) {
    adouble val1(1.);
    adouble val2(2.);
    constexpr auto res = val1.d2(val1, val2);
    EXPECT_EQ(res.size(), 1);
    EXPECT_EQ(res[0], 1.);
}

TEST(adhoc2, InitialMany) {
    adouble val1(1.);
    adouble val2(2.);
    adouble val3(3.);
    constexpr auto res = val1.dmany(val1, val2, val3);
    EXPECT_EQ(res.size(), 1);
    EXPECT_EQ(res[0], 1.);
}

TEST(adhoc2, InitialManyInterface) {
    adouble val1(1.);
    adouble val2(2.);
    adouble val3(3.);
    constexpr auto res = val1.dinterface(val1, val2, val3);
    EXPECT_EQ(res.size(), 3);
    EXPECT_EQ(res[0], 1.);
    EXPECT_EQ(res[1], 0.);
    EXPECT_EQ(res[2], 0.);
}

} // namespace adhoc2