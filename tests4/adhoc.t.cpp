#include <adhoc.hpp>
#include <base.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <numeric>

namespace adhoc4 {

// we use f' = 1 + f^2
void TE_tan(const double &value, std::vector<double> &output,
            std::size_t loc = 0) {
    if (loc == 0) {
        output[loc] = tan(value);
        loc++;
    }

    double fval = output[0];
    double fval_sq = fval * fval;

    std::vector<int> fval_coeffs(output.size() + 1, 0);
    fval_coeffs[1] = 1;

    // highest coefficient (1 at the moment) is an odd power
    bool odd = true;
    while (loc != output.size()) {
        for (std::size_t i = odd; i <= (loc + 1); i += 2) {
            int current_coeff = fval_coeffs[i];
            fval_coeffs[i] = 0;
            if (i) {
                int temp = current_coeff * (int)i;
                fval_coeffs[i - 1] += temp;
                fval_coeffs[i + 1] += temp;
            }
        }
        odd = !odd;

        double fval_der = 0;
        for (std::size_t i = (loc + 1); i > (std::size_t)odd; i -= 2) {
            fval_der += fval_coeffs[i];
            fval_der *= fval_sq;
        }

        fval_der += fval_coeffs[odd];

        if (odd)
            fval_der *= fval;

        output[loc] = fval_der;
        loc++;
    }
}

TEST(Base, Tan) {
    std::cout.precision(std::numeric_limits<double>::max_digits10);
    double in = 0.35;
    double lnin = std::tan(in);

    // std::cout << lnin << std::endl;

    std::vector<double> output(6);
    TE_tan(in, output, 0);

    std::cout << output[0] << std::endl;
    std::cout << output[1] << std::endl;
    std::cout << output[2] << std::endl;
    std::cout << output[3] << std::endl;
    std::cout << output[4] << std::endl;
    std::cout << output[5] << std::endl;
}

double der1_1(double x, double epsilon) {
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
}

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
    EXPECT_NEAR(results1[1], results2[1], 1e-13);
    EXPECT_NEAR(results1[2], results2[2], 1e-10);
    EXPECT_NEAR(results1[3], results2[3], 1e-8);
    EXPECT_NEAR(results1[4], results2[4], 1e-5);
    EXPECT_NEAR(results1[5], results2[5], 1e-3);
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
    EXPECT_NEAR(results1[4], results2[4], 1e-5);
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
    EXPECT_NEAR(results1[3], results2[3], 1e-9);
    EXPECT_NEAR(results1[4], results2[4], 1e-6);
    EXPECT_NEAR(results1[5], results2[5], 1e-2);
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
