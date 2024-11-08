#include <adhoc.hpp>
#include <base.hpp>

#include "finite_differences.hpp"
#include "test_tools.hpp"

#include <cmath>

int main() {
    TEST_START;
    using namespace adhoc4;

    // CompEllint1
    {
        std::array<double, 10> results1;
        results1.fill(0.);
        double val = 0.32;
        double res = comp_ellint_1_t<double>::v(val);
        comp_ellint_1_t<double>::d<6>(res, val, results1);

        std::function<double(double)> lambdainput = [](double d) {
            return std::comp_ellint_1(d);
        };
        double epsilon = 0.01;
        auto results2 = finite_differences(val, epsilon, lambdainput);
        EXPECT_NEAR_ABS(results1[0], results2[0], 1e-12);
        EXPECT_NEAR_ABS(results1[1], results2[1], 1e-11);
        EXPECT_NEAR_ABS(results1[2], results2[2], 1e-8);
        EXPECT_NEAR_ABS(results1[3], results2[3], 1e-6);
        EXPECT_NEAR_ABS(results1[4], results2[4], 1e-4);
        EXPECT_NEAR_ABS(results1[5], results2[5], 1e-1);
        EXPECT_EQUAL(results1[6], 0.);
        EXPECT_EQUAL(results1[7], 0.);
        EXPECT_EQUAL(results1[8], 0.);
        EXPECT_EQUAL(results1[9], 0.);
    }

    // CompEllint2
    {
        std::array<double, 10> results1;
        results1.fill(0.);
        double val = 0.32;
        double res = comp_ellint_2_t<double>::v(val);
        comp_ellint_2_t<double>::d<6>(res, val, results1);

        std::function<double(double)> lambdainput = [](double d) {
            return std::comp_ellint_2(d);
        };
        double epsilon = 0.01;
        auto results2 = finite_differences(val, epsilon, lambdainput);
        EXPECT_NEAR_ABS(results1[0], results2[0], 1e-14);
        EXPECT_NEAR_ABS(results1[1], results2[1], 1e-10);
        EXPECT_NEAR_ABS(results1[2], results2[2], 1e-9);
        EXPECT_NEAR_ABS(results1[3], results2[3], 1e-6);
        EXPECT_NEAR_ABS(results1[4], results2[4], 1e-4);
        EXPECT_NEAR_ABS(results1[5], results2[5], 1e-1);
        EXPECT_EQUAL(results1[6], 0.);
        EXPECT_EQUAL(results1[7], 0.);
        EXPECT_EQUAL(results1[8], 0.);
        EXPECT_EQUAL(results1[9], 0.);
    }

    // RiemannZeta
    {
        {
            std::array<double, 10> results1;
            results1.fill(0.);
            double val = 0.32;
            double res = riemann_zeta_t<double>::v(val);
            riemann_zeta_t<double>::d<6>(res, val, results1);

            std::function<double(double)> lambdainput = [](double d) {
                return std::riemann_zeta(d);
            };
            double epsilon = 0.01;
            auto results2 = finite_differences(val, epsilon, lambdainput);
            EXPECT_LESS_THAN(
                std::abs((results2[0] - results1[0]) / results2[0]), 1e-11);
            EXPECT_LESS_THAN(
                std::abs((results2[1] - results1[1]) / results2[1]), 1e-11);
            EXPECT_LESS_THAN(
                std::abs((results2[2] - results1[2]) / results2[2]), 1e-8);
            EXPECT_LESS_THAN(
                std::abs((results2[3] - results1[3]) / results2[3]), 1e-7);
            EXPECT_LESS_THAN(
                std::abs((results2[4] - results1[4]) / results2[4]), 1e-7);
            EXPECT_LESS_THAN(
                std::abs((results2[5] - results1[5]) / results2[5]), 1e-5);
            EXPECT_EQUAL(results1[6], 0.);
            EXPECT_EQUAL(results1[7], 0.);
            EXPECT_EQUAL(results1[8], 0.);
            EXPECT_EQUAL(results1[9], 0.);
        }

        {
            std::array<double, 10> results1;
            results1.fill(0.);
            double val = -0.32;
            double res = riemann_zeta_t<double>::v(val);
            riemann_zeta_t<double>::d<6>(res, val, results1);

            std::function<double(double)> lambdainput = [](double d) {
                return std::riemann_zeta(d);
            };
            double epsilon = 0.01;
            auto results2 = finite_differences(val, epsilon, lambdainput);
            EXPECT_LESS_THAN(
                std::abs((results2[0] - results1[0]) / results2[0]), 1e-13);
            EXPECT_LESS_THAN(
                std::abs((results2[1] - results1[1]) / results2[1]), 1e-11);
            EXPECT_LESS_THAN(
                std::abs((results2[2] - results1[2]) / results2[2]), 1e-9);
            EXPECT_LESS_THAN(
                std::abs((results2[3] - results1[3]) / results2[3]), 1e-7);
            EXPECT_LESS_THAN(
                std::abs((results2[4] - results1[4]) / results2[4]), 1e-6);
            EXPECT_LESS_THAN(
                std::abs((results2[5] - results1[5]) / results2[5]), 1e-4);
            EXPECT_EQUAL(results1[6], 0.);
            EXPECT_EQUAL(results1[7], 0.);
            EXPECT_EQUAL(results1[8], 0.);
            EXPECT_EQUAL(results1[9], 0.);
        }
    }

    TEST_END;
}
