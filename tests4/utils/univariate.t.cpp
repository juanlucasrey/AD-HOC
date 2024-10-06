#include <utils/univariate.hpp>

#include <combinatorics/factorial.hpp>

#include <gtest/gtest.h>

namespace adhoc4 {

namespace {
template <class Array>
void compare_array(Array const &in1, Array const &in2, double tolerance) {
    for (std::size_t i = 0; i < in1.size(); i++) {
        EXPECT_NEAR(in1[i], in2[i], tolerance);
    }
}
} // namespace

TEST(UtilsUnivariate, ElevateUnivariate) {
    constexpr std::size_t MaxOrder = 5;
    std::array<double, MaxOrder> ua = {0.1, 0.2, 0.3, 0.5, 0.7};
    std::array<double, MaxOrder> ua_elevated{0};
    elevate_univariate<0, 1>(ua, ua_elevated);
    EXPECT_EQ(ua_elevated, ua);
    elevate_univariate<1, 2>(ua, ua_elevated);
    {
        std::array<double, MaxOrder> res = {0., 0.01, 0.04, 0.1, 0.22};
        compare_array(ua_elevated, res, 1e-14);

        std::array<double, MaxOrder> ua_elevated2{0};
        elevate_univariate<0, 2>(ua, ua_elevated2);
        EXPECT_EQ(ua_elevated, ua_elevated2);
    }
    elevate_univariate<2, 3>(ua, ua_elevated);
    {
        std::array<double, MaxOrder> res = {0., 0., 0.001, 0.006, 0.021};
        compare_array(ua_elevated, res, 1e-14);

        std::array<double, MaxOrder> ua_elevated2{0};
        elevate_univariate<0, 3>(ua, ua_elevated2);
        EXPECT_EQ(ua_elevated, ua_elevated2);

        std::array<double, MaxOrder> ua_elevated3{0};
        elevate_univariate<0, 1>(ua, ua_elevated3);
        elevate_univariate<1, 3>(ua, ua_elevated3);
        EXPECT_EQ(ua_elevated, ua_elevated3);
    }
    elevate_univariate<3, 4>(ua, ua_elevated);
    {
        std::array<double, MaxOrder> res = {0., 0., 0., 0.0001, 0.0008};
        compare_array(ua_elevated, res, 1e-14);

        std::array<double, MaxOrder> ua_elevated2{0};
        elevate_univariate<0, 4>(ua, ua_elevated2);
        EXPECT_EQ(ua_elevated, ua_elevated2);

        std::array<double, MaxOrder> ua_elevated3{0};
        elevate_univariate<0, 1>(ua, ua_elevated3);
        elevate_univariate<1, 4>(ua, ua_elevated3);
        EXPECT_EQ(ua_elevated, ua_elevated3);

        std::array<double, MaxOrder> ua_elevated4{0};
        elevate_univariate<0, 2>(ua, ua_elevated4);
        elevate_univariate<2, 4>(ua, ua_elevated4);
        EXPECT_EQ(ua_elevated, ua_elevated4);
    }
    elevate_univariate<4, 5>(ua, ua_elevated);
    {
        std::array<double, MaxOrder> res = {0., 0., 0., 0., 0.00001};
        compare_array(ua_elevated, res, 1e-14);

        std::array<double, MaxOrder> ua_elevated2{0};
        elevate_univariate<0, 5>(ua, ua_elevated2);
        EXPECT_EQ(ua_elevated, ua_elevated2);

        std::array<double, MaxOrder> ua_elevated3{0};
        elevate_univariate<0, 1>(ua, ua_elevated3);
        elevate_univariate<1, 5>(ua, ua_elevated3);
        EXPECT_EQ(ua_elevated, ua_elevated3);

        std::array<double, MaxOrder> ua_elevated4{0};
        elevate_univariate<0, 2>(ua, ua_elevated4);
        elevate_univariate<2, 5>(ua, ua_elevated4);
        EXPECT_EQ(ua_elevated, ua_elevated4);

        std::array<double, MaxOrder> ua_elevated5{0};
        elevate_univariate<0, 3>(ua, ua_elevated5);
        elevate_univariate<3, 5>(ua, ua_elevated5);
        EXPECT_EQ(ua_elevated, ua_elevated5);
    }
    elevate_univariate<5, 6>(ua, ua_elevated);
    {
        std::array<double, MaxOrder> res = {0., 0., 0., 0., 0.};
        compare_array(ua_elevated, res, 1e-14);

        std::array<double, MaxOrder> ua_elevated2{0};
        elevate_univariate<0, 6>(ua, ua_elevated2);
        EXPECT_EQ(ua_elevated, ua_elevated2);
    }
}

TEST(UtilsUnivariate, ToTaylorCoeff) {
    constexpr std::size_t MaxOrder = 5;
    std::array<double, MaxOrder> ua = {0.1, 0.2, 0.3, 0.5, 0.7};
    derivative_to_taylor(ua);

    std::array<double, MaxOrder> res = {0.1, 0.4, 1.8, 12., 84.};
    compare_array(ua, res, 1e-14);
}

} // namespace adhoc4
