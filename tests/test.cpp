#include <adhoc.hpp>

#include <gtest/gtest.h>

namespace adhoc {

TEST(Test, Initial) {
    adouble val(1.);
    auto valplus = val + 1.0;
    auto valmult = valplus * 2.0;
    auto valmult2 = (valmult + 2.0) * 3.0;
    // double result = valmult2;
    // double deriv = d(result) / d(val);
    // double deriv = d(result, val);
    // EXPECT_EQ(result, 18.);
    // EXPECT_EQ(deriv, 0.);
}

TEST(Test, doperator) {
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

TEST(Test, doperatorplus) {
    adouble val1(1.);
    adouble val2(2.);
    auto valplus = val1 + 1.0;
    double res1 = valplus.d(val1);
    double res2 = valplus.d(val2);
    EXPECT_EQ(res1, 1.);
    EXPECT_EQ(res2, 0.);
}

TEST(Test, doperatormult) {
    adouble val1(1.);
    adouble val2(2.);
    auto valmul = val1 * 2.0;
    double res1 = valmul.d(val1);
    double res2 = valmul.d(val2);
    EXPECT_EQ(res1, 2.);
    EXPECT_EQ(res2, 0.);
}

TEST(Test, activeadd) {
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

TEST(Test, activemul) {
    adouble val1(1.);
    adouble val2(2.);
    auto valsum = val1 * val2;

    double res1 = valsum.d(val1);
    double res2 = valsum.d(val2);
    EXPECT_EQ(res1, 2.);
    EXPECT_EQ(res2, 1.);

    // auto valsum2 = (val1 * 3.0) + val1;

    // double res1_2 = valsum2.d(val1);
    // double res2_2 = valsum2.d(val2);
    // EXPECT_EQ(res1_2, 4.);
    // EXPECT_EQ(res2_2, 0.);
}

} // namespace adhoc