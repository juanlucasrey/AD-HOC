#include <adhoc.hpp>
#include <base.hpp>

#include <gtest/gtest.h>

#include <cmath>

namespace adhoc3 {

TEST(Base, InvDer) {

    double val = 100.0;
    double res = 1.0 / val;

    {
        auto derivatives = inv_t<double>::d2<1>(res, val);

        EXPECT_NEAR(derivatives[0], -1.0 / (val * val), 1e-15);
    }

    {
        auto derivatives = inv_t<double>::d2<2>(res, val);

        EXPECT_NEAR(derivatives[0], -1.0 / (val * val), 1e-15);
        EXPECT_NEAR(derivatives[1], 2.0 / (val * val * val), 1e-15);
    }

    {
        auto derivatives = inv_t<double>::d2<3>(res, val);

        EXPECT_NEAR(derivatives[0], -1.0 / (val * val), 1e-15);
        EXPECT_NEAR(derivatives[1], 2.0 / (val * val * val), 1e-15);
        EXPECT_NEAR(derivatives[2], -6.0 / (val * val * val * val), 1e-14);
    }

    {
        auto derivatives = inv_t<double>::d2<4>(res, val);

        EXPECT_NEAR(derivatives[0], -1.0 / (val * val), 1e-15);
        EXPECT_NEAR(derivatives[1], 2.0 / (val * val * val), 1e-15);
        EXPECT_NEAR(derivatives[2], -6.0 / (val * val * val * val), 1e-14);
        EXPECT_NEAR(derivatives[3], 24.0 / (val * val * val * val * val),
                    1e-15);
    }
}

TEST(Base, LogDer) {
    double in = 0.5;
    double lnin = std::log(in);
    auto ders = log_t<double>::d2<4>(lnin, in);

    EXPECT_EQ(ders[0], 2.);
    EXPECT_EQ(ders[1], -4);
    EXPECT_EQ(ders[2], 16);
    EXPECT_EQ(ders[3], -96);
}

TEST(Base, SqrtDer) {
    double in = 0.35;
    double lnin = std::sqrt(in);
    auto ders = sqrt_t<double>::d2<5>(lnin, in);

    // from sympy import *
    // x = Symbol('x')
    // f = sqrt(x)
    // valx = 0.35

    // print(lambdify([x], f.diff(x))(valx))
    // print(lambdify([x], f.diff(x).diff(x))(valx))
    // print(lambdify([x], f.diff(x).diff(x).diff(x))(valx))
    // print(lambdify([x], f.diff(x).diff(x).diff(x).diff(x))(valx))
    // print(lambdify([x], f.diff(x).diff(x).diff(x).diff(x).diff(x))(valx))

    EXPECT_NEAR(ders[0], 0.8451542547285166, 1e-13);
    EXPECT_NEAR(ders[1], -1.2073632210407381, 1e-13);
    EXPECT_NEAR(ders[2], 5.174413804460307, 1e-13);
    EXPECT_NEAR(ders[3], -36.9600986032879, 1e-13);
    EXPECT_NEAR(ders[4], 369.60098603287906, 1e-12);
}

TEST(Base, ErfcDer) {
    double in = 0.35;
    double lnin = std::erfc(in);
    auto ders = erfc_t<double>::d2<5>(lnin, in);

    // from sympy import *
    // x = Symbol('x')
    // f = sqrt(x)
    // valx = 0.35

    // print(lambdify([x], f.diff(x))(valx))
    // print(lambdify([x], f.diff(x).diff(x))(valx))
    // print(lambdify([x], f.diff(x).diff(x).diff(x))(valx))
    // print(lambdify([x], f.diff(x).diff(x).diff(x).diff(x))(valx))
    // print(lambdify([x], f.diff(x).diff(x).diff(x).diff(x).diff(x))(valx))

    // EXPECT_NEAR(ders[0], 0.8451542547285166, 1e-13);
    // EXPECT_NEAR(ders[1], -1.2073632210407381, 1e-13);
    // EXPECT_NEAR(ders[2], 5.174413804460307, 1e-13);
    // EXPECT_NEAR(ders[3], -36.9600986032879, 1e-13);
    // EXPECT_NEAR(ders[4], 369.60098603287906, 1e-12);

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    std::cout << ders[0] << std::endl;
    std::cout << ders[1] << std::endl;
    std::cout << ders[2] << std::endl;
    std::cout << ders[3] << std::endl;
    std::cout << ders[4] << std::endl;
}

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

#ifndef __clang__
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
#endif

} // namespace adhoc3
