#include <evaluate_bwd.hpp>
#include <evaluate_fwd.hpp>
#include <init.hpp>

#include "call_price.hpp"

#include <gtest/gtest.h>

namespace adhoc {

TEST(EvaluateBwd, Exp) {
    auto [v0] = Init<1>();
    auto r = exp(v0);

    auto t = Tape3(r, v0);
    t.val(v0) = 1.0;

    evaluate_fwd(t);
    t.der(r) = 1.0;

    evaluate_bwd(t);
    EXPECT_NEAR(t.der(v0), std::exp(1.0), 1e-13);
}

TEST(EvaluateBwd, Cos) {
    auto [v0] = Init<1>();
    auto r = cos(v0);

    auto t = Tape3(r, v0);
    t.val(v0) = 1.0;

    evaluate_fwd(t);
    t.der(r) = 1.0;

    evaluate_bwd(t);
    EXPECT_NEAR(t.der(v0), -std::sin(1.0), 1e-13);
}

TEST(EvaluateBwd, CosExp) {
    auto [v0] = Init<1>();
    auto r = cos(exp(v0));

    auto t = Tape3(r, v0);
    t.val(v0) = 1.0;

    evaluate_fwd(t);
    t.der(r) = 1.0;

    evaluate_bwd(t);
    EXPECT_NEAR(t.der(v0), -1.1166193174450132, 1e-13);
}

TEST(EvaluateBwd, Univariate) {
    auto [v0] = Init<1>();
    auto r = exp(cos(log(v0)));

    auto t = Tape3(r, v0);
    t.val(v0) = 0.5;
    t.der(r) = 1.0;

    evaluate_fwd(t);
    evaluate_bwd(t);

    EXPECT_NEAR(t.der(v0), 2.757914160416556, 1e-13);
}

TEST(EvaluateBwd, Univariate2) {
    auto [v0] = Init<1>();
    auto r = exp(v0);

    auto t = Tape3(r, v0);
    t.val(v0) = 0.5;
    t.der(r) = 1.0;

    evaluate_fwd(t);
    evaluate_bwd(t);

    EXPECT_NEAR(t.der(v0), 1.6487212707001282, 1e-13);
}

TEST(EvaluateBwd, Add) {
    auto [v0, v1] = Init<2>();
    auto r = v0 + v1;

    auto t = Tape3(r, v0, v1);
    t.val(v0) = 2.0;
    t.val(v1) = 3.0;
    t.der(r) = 1.0;

    evaluate_fwd(t);
    evaluate_bwd(t);

    EXPECT_NEAR(t.der(v0), 1.0, 1e-13);
    EXPECT_NEAR(t.der(v1), 1.0, 1e-13);
}

TEST(EvaluateBwd, Mult) {
    auto [v0, v1] = Init<2>();
    auto r = v0 * v1;

    auto t = Tape3(r, v0, v1);
    t.val(v0) = 100.0;
    t.val(v1) = 10.0;
    t.der(r) = 1.0;

    evaluate_fwd(t);
    evaluate_bwd(t);

    EXPECT_NEAR(t.der(v0), 10.0, 1e-13);
    EXPECT_NEAR(t.der(v1), 100.0, 1e-13);
}

TEST(EvaluateBwd, MultSame) {
    auto [v0] = Init<1>();
    auto r = v0 * v0;

    auto t = Tape3(r, v0);
    t.val(v0) = 1.5;
    t.der(r) = 1.0;

    evaluate_fwd(t);
    evaluate_bwd(t);

    EXPECT_NEAR(t.der(v0), 3.0, 1e-13);
}

TEST(EvaluateBwd, Div) {
    auto [v0, v1] = Init<2>();
    auto r = v0 / v1;

    auto t = Tape3(r, v0, v1);
    t.val(v0) = 1.5;
    t.val(v1) = 2.5;
    t.der(r) = 1.0;

    evaluate_fwd(t);
    evaluate_bwd(t);

    // from sympy import *
    // x = Symbol('x')
    // y = Symbol('y')
    // f = x / y
    // print(lambdify([x, y], f.diff(x))(1.5, 2.5))
    // print(lambdify([x, y], f.diff(y))(1.5, 2.5))

    EXPECT_NEAR(t.der(v0), 0.4, 1e-13);
    EXPECT_NEAR(t.der(v1), -0.24, 1e-13);
}

TEST(EvaluateBwd, RepeatVar) {
    auto [v0, v1] = Init<2>();
    auto m = v0 * v1;
    auto r = m * m;

    auto t = Tape3(r, v0, v1);
    t.val(v0) = 0.5;
    t.val(v1) = 2.5;
    t.der(r) = 1.0;

    evaluate_fwd(t);
    evaluate_bwd(t);

    EXPECT_NEAR(t.der(v0), 6.25, 1e-13);
    EXPECT_NEAR(t.der(v1), 1.25, 1e-13);
}

TEST(EvaluateBwd, BivariateCutLeaf1) {
    auto [v0, v1] = Init<2>();
    auto r = cos(v0) * cos(v1);

    auto t = Tape3(r, v0);
    t.val(v0) = 0.5;
    t.val(v1) = 2.;
    t.der(r) = 1.0;

    evaluate_fwd(t);
    evaluate_bwd(t);

    EXPECT_NEAR(t.der(v0), -std::sin(0.5) * std::cos(2.), 1e-13);
}

TEST(EvaluateBwd, BivariateCutLeaf2) {
    auto [v0, v1] = Init<2>();
    auto r = cos(v0) * cos(v1);

    auto t = Tape3(r, v1);
    t.val(v0) = 0.5;
    t.val(v1) = 2.;
    t.der(r) = 1.0;

    evaluate_fwd(t);
    evaluate_bwd(t);

    EXPECT_NEAR(t.der(v1), std::cos(0.5) * -std::sin(2.), 1e-13);
}

TEST(EvaluateBwd, SimpleEvaluate) {
    auto [v0, v1] = Init<2>();
    auto r = (v0 * v1) + (v0);

    auto t = Tape3(r, v0, v1);
    t.val(v0) = 1.0;
    t.val(v1) = 2.0;
    t.der(r) = 1.0;

    evaluate_fwd(t);
    evaluate_bwd(t);

    EXPECT_NEAR(t.der(v0), 3.0, 1e-15);
    EXPECT_NEAR(t.der(v1), 1.0, 1e-15);
}

TEST(EvaluateBwd, BlackScholes) {
    double result = 0.;

    {
        double S = 100.0;
        double K = 102.0;
        double v = 0.15;
        double T = 0.5;
        result = call_price(S, K, v, T);
    }

    auto [S, K, v, T] = Init<4>();
    auto r = call_price(S, K, v, T);

    auto t = Tape3(r, S, K, v, T);
    t.val(S) = 100.0;
    t.val(K) = 102.0;
    t.val(v) = 0.15;
    t.val(T) = 0.5;
    t.der(r) = 1.0;

    evaluate_fwd(t);
    evaluate_bwd(t);

    EXPECT_NEAR(t.der(S), 0.33961663008862131, 1e-14);
    EXPECT_NEAR(t.der(K), -0.38387614859866631, 1e-14);
    EXPECT_NEAR(t.der(v), -30.580208040388474, 1e-14);
    EXPECT_NEAR(t.der(T), -4.5870312060582705, 1e-14);
}

TEST(EvaluateBwd, ComplexEvaluate) {
    auto [val1, val2, val3, val4, val5] = Init<5>();
    auto exp2 = exp(val2);
    auto tmp1 = val1 + exp2;
    auto valsum = tmp1 + val3;
    auto val12 = val1 * val2;
    auto valprod = val12 + val4;
    auto val22 = val2 * val2;
    auto valprod2 = val22 * val5;
    auto valprod3 = valprod * valprod2;
    auto res = valsum * valprod3;

    auto t = Tape3(res, val1, val2, val3);
    t.val(val1) = 1.0;
    t.val(val2) = 2.0;
    t.val(val3) = 2.0;
    t.val(val4) = 3.0;
    t.val(val5) = 4.0;
    t.der(res) = 1.0;

    evaluate_fwd(t);
    evaluate_bwd(t);

    EXPECT_NEAR(t.der(val1), 412.44979516578081, 1e-14);
    EXPECT_NEAR(t.der(val2), 1588.4738734117946, 1e-14);
    EXPECT_NEAR(t.der(val3), 80, 1e-14);

    auto res2 = valprod3 * valsum;
    auto t2 = Tape3(res2, val1, val2, val3);
    t2.val(val1) = 1.0;
    t2.val(val2) = 2.0;
    t2.val(val3) = 2.0;
    t2.val(val4) = 3.0;
    t2.val(val5) = 4.0;
    t2.der(res2) = 1.0;

    evaluate_fwd(t2);
    evaluate_bwd(t2);

    EXPECT_NEAR(t2.der(val1), 412.44979516578081, 1e-14);
    EXPECT_NEAR(t2.der(val2), 1588.4738734117946, 1e-14);
    EXPECT_NEAR(t2.der(val3), 80, 1e-14);
}

// TEST(EvaluateBwd, Derivative2nd) {
//     adouble x(2);
//     adouble y(3);
//     adouble z(5);
//     auto x2 = x * x;
//     auto y2 = y * y;
//     auto z2 = z * z;
//     auto temp = x2 * y2;
//     auto temp2 = temp * z2;
//     auto res = temp * temp2;

//     // f1(temp, temp2) = temp * temp2
//     // df1dtemp = dtemp/dtemp * temp2 + temp * dtemp2/dtemp
//     // df1dtemp_2 = dtemp/dtemp_2 * temp2 + temp * dtemp2/dtemp_2 +
//     dtemp/dtemp * dtemp2/dtemp double dresdtemp = temp2.v(); double
//     dresdtemp2 = temp.v(); double d2_res_d_temp_2 = 0; double
//     d2_res_d_temp2_2 = 0; double d2_res_dtemp2_dtemp = 1;

//     // f2(temp, z2) = temp * z2
//     // f2(temp, z2) = temp * (temp * z2) = f1(temp, f2(temp, z2))
//     // dfdtemp_2 = dtemp/dtemp_2 * (temp * z2) + temp * d(temp * z2)/dtemp_2
//     + dtemp/dtemp * d(temp * z2)/dtemp dresdtemp += dresdtemp2 * z2.v();
//     double d2resdz2 = dresdtemp2 * temp.v();
//     d2_res_d_temp_2 += 0;
//     double d2_res_d_z2_2 = 0;
//     double d2_res_dz2_dtemp = 0;
// }

// from sympy import *
// x = Symbol('x')
// y = Symbol('y')
// z = Symbol('z')
// temp = ((x * x) * (y * y))
// temp2 = temp * (z * z)
// f = temp * temp2
// print(f)
// dfdx = f.diff(x)
// print(dfdx)
// d2fdx2 = dfdx.diff(x)
// print(d2fdx2)

// dfdx = lambdify([x, y, z], dfdx)
// d2fdx2 = lambdify([x, y, z], d2fdx2)
// # dfdy = lambdify([x, y], dfdy)
// # print(dfdx(1.5, 2.5))
// print(dfdx(2, 3, 5))
// print(d2fdx2(2, 3, 5))

} // namespace adhoc