#include <adhoc.hpp>
#include <base.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <numeric>

namespace adhoc4 {

/*double der1_1(double x, double epsilon) {
    double f1 = std::comp_ellint_1(x + epsilon);
    double f2 = std::comp_ellint_1(x - epsilon);
    return 0.5 * (f1 - f2) / epsilon;
}

double der1_2(double x, double epsilon) {
    double f1 = std::comp_ellint_2(x + epsilon);
    double f2 = std::comp_ellint_2(x - epsilon);
    return 0.5 * (f1 - f2) / epsilon;
}

double der2_1(double x, double epsilon) {
    double f1 = std::comp_ellint_1(x + epsilon);
    double f = std::comp_ellint_1(x);
    double f2 = std::comp_ellint_1(x - epsilon);
    return 0.5 * (f1 + f2 - 2 * f) / (epsilon * epsilon);
}

double der1c_1(double x) {
    double e = std::comp_ellint_2(x);
    double f = std::comp_ellint_1(x);
    double num = (x * x - 1) * f + e;
    double den = x * (1 - x * x);
    return num / den;
}

double der1c_2(double x) {
    double num = std::comp_ellint_2(x) - std::comp_ellint_1(x);
    double den = x;
    return num / den;
}

TEST(Base, CompEllint1) {
    double x1 = 0.02;
    double epsilon = 0.00001;
    double val = std::comp_ellint_1(x1);
    std::cout.precision(std::numeric_limits<double>::max_digits10);
    std::cout << val << std::endl;

    std::cout << der1_1(x1, epsilon) << std::endl;
    std::cout << der1c_1(x1) << std::endl;

    std::cout << der1_2(x1, epsilon) << std::endl;
    std::cout << der1c_2(x1) << std::endl;
    // std::cout << der2(x1, epsilon) << std::endl;
    // std::cout << std::comp_ellint_1(0) << std::endl;
    // constexpr double π{std::numbers::pi};

    // std::cout << "K(0) ≈ " << std::comp_ellint_1(0) << '\n'
    //           << "π/2 ≈ " << π / 2 << '\n'
    //           << "K(0.5) ≈ " << std::comp_ellint_1(0.5) << '\n'
    //           << "F(0.5, π/2) ≈ " << std::ellint_1(0.5, π / 2) << '\n'
    //           << "The period of a pendulum length 1m at 10° initial angle ≈ "
    //           << 4 * std::sqrt(1 / 9.80665) *
    //                  std::comp_ellint_1(std::sin(π / 18 / 2))
    //           << "s,\n"
    //              "whereas the linear approximation gives ≈ "
    //           << 2 * π * std::sqrt(1 / 9.80665) << '\n';
}*/

std::array<double, 6> finite_differences(double x, double epsilon,
                                         std::function<double(double)> func) {
    std::vector<double> fvals{
        func(x - 5. * epsilon), func(x - 4. * epsilon), func(x - 3. * epsilon),
        func(x - 2. * epsilon), func(x - 1. * epsilon), func(x),
        func(x + 1. * epsilon), func(x + 2. * epsilon), func(x + 3. * epsilon),
        func(x + 4. * epsilon), func(x + 5. * epsilon)};

    std::vector<double> weights1 = {0.,        1. / 280., -4. / 105., 1. / 5.,
                                    -4. / 5.,  0.,        4. / 5.,    -1. / 5.,
                                    4. / 105., -1. / 280, 0.};

    std::vector<double> weights2 = {0.,        -1. / 560.,  8. / 315., -1. / 5.,
                                    8. / 5.,   -205. / 72., 8. / 5.,   -1. / 5.,
                                    8. / 315., -1. / 560.,  0.};

    std::vector<double> weights3 = {
        0.,         -7. / 240.,  3. / 10.,  -169. / 120., 61. / 30., 0.,
        -61. / 30., 169. / 120., -3. / 10., 7. / 240.,    0.};

    std::vector<double> weights4 = {
        0.,          7. / 240.,  -2. / 5., 169. / 60., -122. / 15., 91. / 8.,
        -122. / 15., 169. / 60., -2. / 5., 7. / 240.,  0.};

    std::vector<double> weights5 = {
        -13. / 288., 19. / 36., -87. / 32., 13. / 2.,   -323. / 48., 0.,
        323. / 48.,  -13. / 2., 87. / 32.,  -19. / 36., 13. / 288.};

    std::vector<double> weights6 = {
        13. / 240., -19. / 24., 87. / 16., -39. / 2.,  323. / 8., -1023. / 20.,
        323. / 8.,  -39. / 2.,  87. / 16., -19. / 24., 13. / 240.};

    std::array<double, 6> result;
    result.fill(0);
    result[0] = std::transform_reduce(fvals.begin(), fvals.end(),
                                      weights1.begin(), 0.) /
                epsilon;

    result[1] = std::transform_reduce(fvals.begin(), fvals.end(),
                                      weights2.begin(), 0.) /
                std::pow(epsilon, 2);

    result[2] = std::transform_reduce(fvals.begin(), fvals.end(),
                                      weights3.begin(), 0.) /
                std::pow(epsilon, 3);

    result[3] = std::transform_reduce(fvals.begin(), fvals.end(),
                                      weights4.begin(), 0.) /
                std::pow(epsilon, 4);

    result[4] = std::transform_reduce(fvals.begin(), fvals.end(),
                                      weights5.begin(), 0.) /
                std::pow(epsilon, 5);

    result[5] = std::transform_reduce(fvals.begin(), fvals.end(),
                                      weights6.begin(), 0.) /
                std::pow(epsilon, 6);

    return result;
}

TEST(UnivariateFunctions, Inverse) {

    std::array<double, 6> results1;
    double val = 100.0;
    double res = inv_t<double>::v(val);
    inv_t<double>::d<6>(res, val, results1);

    std::function<double(double)> lambdainput = [](double d) { return 1. / d; };
    double epsilon = 0.1;
    auto results2 = finite_differences(val, epsilon, lambdainput);
    EXPECT_NEAR(results1[0], results2[0], 1e-15);
    EXPECT_NEAR(results1[1], results2[1], 1e-15);
    EXPECT_NEAR(results1[2], results2[2], 3e-15);
    EXPECT_NEAR(results1[3], results2[3], 2e-13);
    EXPECT_NEAR(results1[4], results2[4], 1e-12);
    EXPECT_NEAR(results1[5], results2[5], 1e-10);
}

TEST(UnivariateFunctions, Exponential) {

    std::array<double, 6> results1;
    double val = .32;
    double res = exp_t<double>::v(val);
    exp_t<double>::d<6>(res, val, results1);

    std::function<double(double)> lambdainput = [](double d) {
        return std::exp(d);
    };
    double epsilon = 0.1;
    auto results2 = finite_differences(val, epsilon, lambdainput);
    EXPECT_NEAR(results1[0], results2[0], 1e-10);
    EXPECT_NEAR(results1[1], results2[1], 1e-11);
    EXPECT_NEAR(results1[2], results2[2], 1e-7);
    EXPECT_NEAR(results1[3], results2[3], 1e-8);
    EXPECT_NEAR(results1[4], results2[4], 1e-7);
    EXPECT_NEAR(results1[5], results2[5], 1e-8);
}

TEST(UnivariateFunctions, Logarithm) {

    std::array<double, 6> results1;
    double val = 1.32;
    double res = log_t<double>::v(val);
    log_t<double>::d<6>(res, val, results1);

    std::function<double(double)> lambdainput = [](double d) {
        return std::log(d);
    };
    double epsilon = 0.01;
    auto results2 = finite_differences(val, epsilon, lambdainput);
    EXPECT_NEAR(results1[0], results2[0], 1e-14);
    EXPECT_NEAR(results1[1], results2[1], 1e-12);
    EXPECT_NEAR(results1[2], results2[2], 1e-10);
    EXPECT_NEAR(results1[3], results2[3], 1e-7);
    EXPECT_NEAR(results1[4], results2[4], 1e-5);
    EXPECT_NEAR(results1[5], results2[5], 1e-2);
}

TEST(UnivariateFunctions, SquareRoot) {

    std::array<double, 6> results1;
    double val = 1.32;
    double res = sqrt_t<double>::v(val);
    sqrt_t<double>::d<6>(res, val, results1);

    std::function<double(double)> lambdainput = [](double d) {
        return std::sqrt(d);
    };
    double epsilon = 0.01;
    auto results2 = finite_differences(val, epsilon, lambdainput);
    EXPECT_NEAR(results1[0], results2[0], 1e-14);
    EXPECT_NEAR(results1[1], results2[1], 1e-11);
    EXPECT_NEAR(results1[2], results2[2], 1e-9);
    EXPECT_NEAR(results1[3], results2[3], 1e-6);
    EXPECT_NEAR(results1[4], results2[4], 1e-4);
    EXPECT_NEAR(results1[5], results2[5], 1e-2);
}

TEST(UnivariateFunctions, Sine) {

    std::array<double, 6> results1;
    double val = 1.32;
    double res = sin_t<double>::v(val);
    sin_t<double>::d<6>(res, val, results1);

    std::function<double(double)> lambdainput = [](double d) {
        return std::sin(d);
    };
    double epsilon = 0.01;
    auto results2 = finite_differences(val, epsilon, lambdainput);
    EXPECT_NEAR(results1[0], results2[0], 1e-14);
    EXPECT_NEAR(results1[1], results2[1], 1e-11);
    EXPECT_NEAR(results1[2], results2[2], 1e-9);
    EXPECT_NEAR(results1[3], results2[3], 1e-6);
    EXPECT_NEAR(results1[4], results2[4], 1e-5);
    EXPECT_NEAR(results1[5], results2[5], 1e-2);
}

TEST(UnivariateFunctions, Cosine) {

    std::array<double, 6> results1;
    double val = 1.32;
    double res = cos_t<double>::v(val);
    cos_t<double>::d<6>(res, val, results1);

    std::function<double(double)> lambdainput = [](double d) {
        return std::cos(d);
    };
    double epsilon = 0.01;
    auto results2 = finite_differences(val, epsilon, lambdainput);
    EXPECT_NEAR(results1[0], results2[0], 1e-14);
    EXPECT_NEAR(results1[1], results2[1], 1e-12);
    EXPECT_NEAR(results1[2], results2[2], 1e-10);
    EXPECT_NEAR(results1[3], results2[3], 1e-7);
    EXPECT_NEAR(results1[4], results2[4], 1e-6);
    EXPECT_NEAR(results1[5], results2[5], 1e-2);
}

TEST(UnivariateFunctions, Tangent) {

    std::array<double, 6> results1;
    double val = 1.32;
    double res = tan_t<double>::v(val);
    tan_t<double>::d<6>(res, val, results1);

    std::function<double(double)> lambdainput = [](double d) {
        return std::tan(d);
    };
    double epsilon = 0.01;
    auto results2 = finite_differences(val, epsilon, lambdainput);
    EXPECT_NEAR(results1[0], results2[0], 1e-7);
    EXPECT_NEAR(results1[1], results2[1], 1e-6);
    EXPECT_NEAR(results1[2], results2[2], 1e-2);
    EXPECT_NEAR(results1[3], results2[3], 1e-1);
    EXPECT_NEAR(results1[4], results2[4], 20);
    EXPECT_NEAR(results1[5], results2[5], 400);

    // from sympy import *
    // x = Symbol('x')
    // r = tan(x)
    // valx = 1.32
    // print(lambdify([x], r.diff(x))(valx))
    // print(lambdify([x], r.diff(x).diff(x))(valx))
    // print(lambdify([x], r.diff(x).diff(x).diff(x))(valx))
    // print(lambdify([x], r.diff(x).diff(x).diff(x).diff(x))(valx))
    // print(lambdify([x], r.diff(x).diff(x).diff(x).diff(x).diff(x))(valx))
    // print(lambdify([x],
    // r.diff(x).diff(x).diff(x).diff(x).diff(x).diff(x))(valx))

    EXPECT_NEAR(results1[0], 16.2361239501548, 1e-15);
    EXPECT_NEAR(results1[1], 126.75047699671539, 1e-15);
    EXPECT_NEAR(results1[2], 1516.7258297481221, 1e-12);
    EXPECT_NEAR(results1[3], 24188.235555132123, 1e-11);
    EXPECT_NEAR(results1[4], 482230.28055004874, 1e-10);
    EXPECT_NEAR(results1[5], 11536771.351260751, 1e-8);
}

TEST(UnivariateFunctions, InverseSine) {

    std::array<double, 6> results1;
    double val = 0.32;
    double res = asin_t<double>::v(val);
    asin_t<double>::d<6>(res, val, results1);

    std::function<double(double)> lambdainput = [](double d) {
        return std::asin(d);
    };
    double epsilon = 0.01;
    auto results2 = finite_differences(val, epsilon, lambdainput);
    EXPECT_NEAR(results1[0], results2[0], 1e-13);
    EXPECT_NEAR(results1[1], results2[1], 1e-11);
    EXPECT_NEAR(results1[2], results2[2], 5e-9);
    EXPECT_NEAR(results1[3], results2[3], 1e-6);
    EXPECT_NEAR(results1[4], results2[4], 1e-5);
    EXPECT_NEAR(results1[5], results2[5], 1e-2);
}

TEST(UnivariateFunctions, HyperbolicSine) {

    std::array<double, 6> results1;
    double val = 1.32;
    double res = sinh_t<double>::v(val);
    sinh_t<double>::d<6>(res, val, results1);

    std::function<double(double)> lambdainput = [](double d) {
        return std::sinh(d);
    };
    double epsilon = 0.01;
    auto results2 = finite_differences(val, epsilon, lambdainput);
    EXPECT_NEAR(results1[0], results2[0], 1e-14);
    EXPECT_NEAR(results1[1], results2[1], 1e-11);
    EXPECT_NEAR(results1[2], results2[2], 1e-9);
    EXPECT_NEAR(results1[3], results2[3], 1e-6);
    EXPECT_NEAR(results1[4], results2[4], 1e-5);
    EXPECT_NEAR(results1[5], results2[5], 3e-2);
}

TEST(UnivariateFunctions, HyperbolicCosine) {

    std::array<double, 6> results1;
    double val = 1.32;
    double res = cosh_t<double>::v(val);
    cosh_t<double>::d<6>(res, val, results1);

    std::function<double(double)> lambdainput = [](double d) {
        return std::cosh(d);
    };
    double epsilon = 0.01;
    auto results2 = finite_differences(val, epsilon, lambdainput);
    EXPECT_NEAR(results1[0], results2[0], 1e-14);
    EXPECT_NEAR(results1[1], results2[1], 1e-11);
    EXPECT_NEAR(results1[2], results2[2], 1e-9);
    EXPECT_NEAR(results1[3], results2[3], 1e-6);
    EXPECT_NEAR(results1[4], results2[4], 3e-5);
    EXPECT_NEAR(results1[5], results2[5], 3e-2);
}

TEST(UnivariateFunctions, HyperbolicTangent) {

    std::array<double, 6> results1;
    double val = 1.32;
    double res = tanh_t<double>::v(val);
    tanh_t<double>::d<6>(res, val, results1);

    std::function<double(double)> lambdainput = [](double d) {
        return std::tanh(d);
    };
    double epsilon = 0.01;
    auto results2 = finite_differences(val, epsilon, lambdainput);
    EXPECT_NEAR(results1[0], results2[0], 1e-7);
    EXPECT_NEAR(results1[1], results2[1], 1e-6);
    EXPECT_NEAR(results1[2], results2[2], 1e-2);
    EXPECT_NEAR(results1[3], results2[3], 1e-1);
    EXPECT_NEAR(results1[4], results2[4], 20);
    EXPECT_NEAR(results1[5], results2[5], 400);

    // from sympy import *
    // x = Symbol('x')
    // r = tanh(x)
    // valx = 1.32
    // print(lambdify([x], r.diff(x))(valx))
    // print(lambdify([x], r.diff(x).diff(x))(valx))
    // print(lambdify([x], r.diff(x).diff(x).diff(x))(valx))
    // print(lambdify([x], r.diff(x).diff(x).diff(x).diff(x))(valx))
    // print(lambdify([x], r.diff(x).diff(x).diff(x).diff(x).diff(x))(valx))
    // print(lambdify([x],
    // r.diff(x).diff(x).diff(x).diff(x).diff(x).diff(x))(valx))

    EXPECT_NEAR(results1[0], 0.24868562068767242, 1e-15);
    EXPECT_NEAR(results1[1], -0.4311133986962329, 1e-15);
    EXPECT_NEAR(results1[2], 0.6236752551298124, 1e-15);
    EXPECT_NEAR(results1[3], -0.43791315708639567, 1e-14);
    EXPECT_NEAR(results1[4], -1.5967929450597538, 1e-14);
    EXPECT_NEAR(results1[5], 9.234671156165486, 1e-13);
}

TEST(UnivariateFunctions, ErrorFunction) {

    std::array<double, 6> results1;
    double val = 1.32;
    double res = erfc_t<double>::v(val);
    erfc_t<double>::d<6>(res, val, results1);

    std::function<double(double)> lambdainput = [](double d) {
        return std::erfc(d);
    };
    double epsilon = 0.01;
    auto results2 = finite_differences(val, epsilon, lambdainput);
    EXPECT_NEAR(results1[0], results2[0], 1e-14);
    EXPECT_NEAR(results1[1], results2[1], 1e-12);
    EXPECT_NEAR(results1[2], results2[2], 1e-10);
    EXPECT_NEAR(results1[3], results2[3], 1e-7);
    EXPECT_NEAR(results1[4], results2[4], 1e-5);
    EXPECT_NEAR(results1[5], results2[5], 1e-2);
}

} // namespace adhoc4
