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
    int id1 = val1.id();
    int id2 = val2.id();
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

} // namespace adhoc2