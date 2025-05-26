#include <utils/univariate.hpp>

#include "../../test_simple/test_simple_include.hpp"

#include <array>

int main() {
    using namespace adhoc;
    // ElevateUnivariate
    {
        constexpr std::size_t MaxOrder = 5;
        std::array<double, MaxOrder> ua = {0.1, 0.2, 0.3, 0.5, 0.7};
        std::array<double, MaxOrder> ua_elevated{0};
        elevate_univariate<0, 1>(ua, ua_elevated);
        EXPECT_EQUAL_ARRAY(ua_elevated, ua);
        elevate_univariate<1, 2>(ua, ua_elevated);
        {
            std::array<double, MaxOrder> res = {0., 0.01, 0.04, 0.1, 0.22};
            EXPECT_NEAR_REL_ARRAY(ua_elevated, res, 1e-14);

            std::array<double, MaxOrder> ua_elevated2{0};
            elevate_univariate<0, 2>(ua, ua_elevated2);
            EXPECT_EQUAL_ARRAY(ua_elevated, ua_elevated2);
        }
        elevate_univariate<2, 3>(ua, ua_elevated);
        {
            std::array<double, MaxOrder> res = {0., 0., 0.001, 0.006, 0.021};
            EXPECT_NEAR_REL_ARRAY(ua_elevated, res, 1e-14);

            std::array<double, MaxOrder> ua_elevated2{0};
            elevate_univariate<0, 3>(ua, ua_elevated2);
            EXPECT_EQUAL_ARRAY(ua_elevated, ua_elevated2);

            std::array<double, MaxOrder> ua_elevated3{0};
            elevate_univariate<0, 1>(ua, ua_elevated3);
            elevate_univariate<1, 3>(ua, ua_elevated3);
            EXPECT_EQUAL_ARRAY(ua_elevated, ua_elevated3);
        }
        elevate_univariate<3, 4>(ua, ua_elevated);
        {
            std::array<double, MaxOrder> res = {0., 0., 0., 0.0001, 0.0008};
            EXPECT_NEAR_REL_ARRAY(ua_elevated, res, 1e-14);

            std::array<double, MaxOrder> ua_elevated2{0};
            elevate_univariate<0, 4>(ua, ua_elevated2);
            EXPECT_EQUAL_ARRAY(ua_elevated, ua_elevated2);

            std::array<double, MaxOrder> ua_elevated3{0};
            elevate_univariate<0, 1>(ua, ua_elevated3);
            elevate_univariate<1, 4>(ua, ua_elevated3);
            EXPECT_EQUAL_ARRAY(ua_elevated, ua_elevated3);

            std::array<double, MaxOrder> ua_elevated4{0};
            elevate_univariate<0, 2>(ua, ua_elevated4);
            elevate_univariate<2, 4>(ua, ua_elevated4);
            EXPECT_EQUAL_ARRAY(ua_elevated, ua_elevated4);
        }
        elevate_univariate<4, 5>(ua, ua_elevated);
        {
            std::array<double, MaxOrder> res = {0., 0., 0., 0., 0.00001};
            EXPECT_NEAR_REL_ARRAY(ua_elevated, res, 1e-14);

            std::array<double, MaxOrder> ua_elevated2{0};
            elevate_univariate<0, 5>(ua, ua_elevated2);
            EXPECT_EQUAL_ARRAY(ua_elevated, ua_elevated2);

            std::array<double, MaxOrder> ua_elevated3{0};
            elevate_univariate<0, 1>(ua, ua_elevated3);
            elevate_univariate<1, 5>(ua, ua_elevated3);
            EXPECT_EQUAL_ARRAY(ua_elevated, ua_elevated3);

            std::array<double, MaxOrder> ua_elevated4{0};
            elevate_univariate<0, 2>(ua, ua_elevated4);
            elevate_univariate<2, 5>(ua, ua_elevated4);
            EXPECT_EQUAL_ARRAY(ua_elevated, ua_elevated4);

            std::array<double, MaxOrder> ua_elevated5{0};
            elevate_univariate<0, 3>(ua, ua_elevated5);
            elevate_univariate<3, 5>(ua, ua_elevated5);
            EXPECT_EQUAL_ARRAY(ua_elevated, ua_elevated5);
        }
        elevate_univariate<5, 6>(ua, ua_elevated);
        {
            std::array<double, MaxOrder> res = {0., 0., 0., 0., 0.};
            EXPECT_NEAR_REL_ARRAY(ua_elevated, res, 1e-14);

            std::array<double, MaxOrder> ua_elevated2{0};
            elevate_univariate<0, 6>(ua, ua_elevated2);
            EXPECT_EQUAL_ARRAY(ua_elevated, ua_elevated2);
        }
    }

    // ToTaylorCoeff
    {
        constexpr std::size_t MaxOrder = 5;
        std::array<double, MaxOrder> ua = {0.1, 0.2, 0.3, 0.5, 0.7};
        derivative_to_taylor(ua);

        std::array<double, MaxOrder> res = {0.1, 0.2 / 2., 0.3 / 6., 0.5 / 24.,
                                            0.7 / 120.};
        EXPECT_NEAR_REL_ARRAY(ua, res, 1e-14);
    }

    TEST_END;
}
