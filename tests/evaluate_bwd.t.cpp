#include <new/evaluate_bwd.hpp>
#include <new/evaluate_fwd.hpp>
#include <new/init.hpp>

#include "call_price.hpp"

#include <gtest/gtest.h>

namespace adhoc2 {

TEST(EvaluateBwd, Exp) {
    auto [v0] = Init<1>();
    auto r = exp(v0);

    auto leaves_and_roots = TapeRootsAndLeafs(r);
    leaves_and_roots.set(v0, 1.0);

    auto intermediate_tape = evaluate_fwd_return_vals(leaves_and_roots);
    auto tape_d = TapeDerivatives(v0, r);
    tape_d.get(r) = 1.0;

    evaluate_bwd(leaves_and_roots, intermediate_tape, tape_d);
    EXPECT_NEAR(tape_d.get(v0), std::exp(1.0), 1e-13);
}

TEST(EvaluateBwd, Cos) {
    auto [v0] = Init<1>();
    auto r = cos(v0);

    auto leaves_and_roots = TapeRootsAndLeafs(r);
    leaves_and_roots.set(v0, 1.0);

    auto intermediate_tape = evaluate_fwd_return_vals(leaves_and_roots);
    auto tape_d = TapeDerivatives(v0, r);
    tape_d.get(r) = 1.0;

    evaluate_bwd(leaves_and_roots, intermediate_tape, tape_d);
    EXPECT_NEAR(tape_d.get(v0), -std::sin(1.0), 1e-13);
}

TEST(EvaluateBwd, CosExp) {
    auto [v0] = Init<1>();
    auto r = cos(exp(v0));

    auto leaves_and_roots = TapeRootsAndLeafs(r);
    leaves_and_roots.set(v0, 1.0);

    auto intermediate_tape = evaluate_fwd_return_vals(leaves_and_roots);
    auto tape_d = TapeDerivatives(v0, r);
    tape_d.get(r) = 1.0;

    evaluate_bwd(leaves_and_roots, intermediate_tape, tape_d);
    EXPECT_NEAR(tape_d.get(v0), -1.1166193174450132, 1e-13);
}

TEST(EvaluateBwd, Univariate) {
    auto [v0] = Init<1>();
    // auto r = exp(cos(v0) * (v1 * v1));
    auto r = exp(cos(log(v0)));

    auto leaves_and_roots = TapeRootsAndLeafs(r);
    leaves_and_roots.set(v0, 0.5);

    auto intermediate_tape = evaluate_fwd_return_vals(leaves_and_roots);
    auto tape_d = TapeDerivatives(v0, r);
    tape_d.get(r) = 1.0;

    evaluate_bwd(leaves_and_roots, intermediate_tape, tape_d);
    EXPECT_NEAR(tape_d.get(v0), 2.757914160416556, 1e-13);
}

TEST(EvaluateBwd, Univariate2) {
    auto [v0] = Init<1>();
    // auto r = exp(cos(v0) * (v1 * v1));
    auto r = exp(v0);

    auto leaves_and_roots = TapeRootsAndLeafs(r);
    leaves_and_roots.set(v0, 0.5);

    auto intermediate_tape = evaluate_fwd_return_vals(leaves_and_roots);
    auto tape_d = TapeDerivatives(v0, r);
    tape_d.get(r) = 1.0;

    evaluate_bwd(leaves_and_roots, intermediate_tape, tape_d);
    EXPECT_NEAR(tape_d.get(v0), 1.6487212707001282, 1e-13);
}

TEST(EvaluateBwd, Add) {
    auto [val1, val2] = Init<2>();
    auto result = val1 + val2;

    auto leaves_and_roots = TapeRootsAndLeafs(result);
    leaves_and_roots.get(val1) = 2.0;
    leaves_and_roots.get(val2) = 3.0;

    auto intermediate_tape = evaluate_fwd_return_vals(leaves_and_roots);
    auto tape_d = TapeDerivatives(result, val1, val2);
    tape_d.get(result) = 1.0;

    evaluate_bwd(leaves_and_roots, intermediate_tape, tape_d);
    EXPECT_NEAR(tape_d.get(val1), 1.0, 1e-13);
    EXPECT_NEAR(tape_d.get(val2), 1.0, 1e-13);
}

TEST(EvaluateBwd, Mult) {
    auto [val1, val2] = Init<2>();
    auto result = val1 * val2;

    auto leaves_and_roots = TapeRootsAndLeafs(result);
    leaves_and_roots.get(val1) = 100.0;
    leaves_and_roots.get(val2) = 10.0;

    auto intermediate_tape = evaluate_fwd_return_vals(leaves_and_roots);
    auto tape_d = TapeDerivatives(result, val1, val2);
    tape_d.get(result) = 1.0;

    evaluate_bwd(leaves_and_roots, intermediate_tape, tape_d);
    EXPECT_NEAR(tape_d.get(val1), 10.0, 1e-13);
    EXPECT_NEAR(tape_d.get(val2), 100.0, 1e-13);
}

TEST(EvaluateBwd, MultSame) {
    auto [val1] = Init<1>();
    auto result = val1 * val1;

    auto leaves_and_roots = TapeRootsAndLeafs(result);
    leaves_and_roots.get(val1) = 1.5;

    auto intermediate_tape = evaluate_fwd_return_vals(leaves_and_roots);
    auto tape_d = TapeDerivatives(result, val1);
    tape_d.get(result) = 1.0;

    evaluate_bwd(leaves_and_roots, intermediate_tape, tape_d);
    EXPECT_NEAR(tape_d.get(val1), 3.0, 1e-13);
}

TEST(EvaluateBwd, Div) {
    auto [val1, val2] = Init<2>();
    auto result = val1 / val2;

    auto leaves_and_roots = TapeRootsAndLeafs(result);
    leaves_and_roots.get(val1) = 1.5;
    leaves_and_roots.get(val2) = 2.5;

    auto intermediate_tape = evaluate_fwd_return_vals(leaves_and_roots);
    auto tape_d = TapeDerivatives(result, val1, val2);
    tape_d.get(result) = 1.0;

    evaluate_bwd(leaves_and_roots, intermediate_tape, tape_d);

    // from sympy import *
    // x = Symbol('x')
    // y = Symbol('y')
    // f = x / y
    // print(lambdify([x, y], f.diff(x))(1.5, 2.5))
    // print(lambdify([x, y], f.diff(y))(1.5, 2.5))

    EXPECT_NEAR(tape_d.get(val1), 0.4, 1e-13);
    EXPECT_NEAR(tape_d.get(val2), -0.24, 1e-13);
}

TEST(EvaluateBwd, RepeatVar) {
    auto [val1, val2] = Init<2>();
    auto m = val1 * val2;
    auto result = m * m;

    auto leaves_and_roots = TapeRootsAndLeafs(result);
    leaves_and_roots.get(val1) = 0.5;
    leaves_and_roots.get(val2) = 2.5;

    auto intermediate_tape = evaluate_fwd_return_vals(leaves_and_roots);
    auto tape_d = TapeDerivatives(result, val1, val2);
    tape_d.get(result) = 1.0;

    evaluate_bwd(leaves_and_roots, intermediate_tape, tape_d);
    EXPECT_NEAR(tape_d.get(val1), 6.25, 1e-13);
    EXPECT_NEAR(tape_d.get(val2), 1.25, 1e-13);
}

TEST(EvaluateBwd, BivariateCutLeaf1) {
    auto [v0, v1] = Init<2>();
    // auto r = exp(cos(v0) * (v1 * v1));
    auto r = cos(v0) * cos(v1);

    auto leaves_and_roots = TapeRootsAndLeafs(r);
    leaves_and_roots.set(v0, 0.5);
    leaves_and_roots.set(v1, 2.);

    auto intermediate_tape = evaluate_fwd_return_vals(leaves_and_roots);
    auto tape_d = TapeDerivatives(v0, r);
    tape_d.get(r) = 1.0;

    evaluate_bwd(leaves_and_roots, intermediate_tape, tape_d);
    EXPECT_NEAR(tape_d.get(v0), -std::sin(0.5) * std::cos(2.), 1e-13);
}

TEST(EvaluateBwd, BivariateCutLeaf2) {
    auto [v0, v1] = Init<2>();
    // auto r = exp(cos(v0) * (v1 * v1));
    auto r = cos(v0) * cos(v1);

    auto leaves_and_roots = TapeRootsAndLeafs(r);
    leaves_and_roots.set(v0, 0.5);
    leaves_and_roots.set(v1, 2.);

    auto intermediate_tape = evaluate_fwd_return_vals(leaves_and_roots);
    auto tape_d = TapeDerivatives(v1, r);
    tape_d.get(r) = 1.0;

    evaluate_bwd(leaves_and_roots, intermediate_tape, tape_d);
    EXPECT_NEAR(tape_d.get(v1), std::cos(0.5) * -std::sin(2.), 1e-13);
}

TEST(EvaluateBwd, SimpleEvaluate) {
    auto [v0, v1] = Init<2>();
    auto r = (v0 * v1) + (v0);

    auto leaves_and_roots = TapeRootsAndLeafs(r);
    leaves_and_roots.set(v0, 1.0);
    leaves_and_roots.set(v1, 2.0);

    auto intermediate_tape = evaluate_fwd_return_vals(leaves_and_roots);
    auto tape_d = TapeDerivatives(v0, v1, r);
    tape_d.get(r) = 1.0;

    evaluate_bwd(leaves_and_roots, intermediate_tape, tape_d);
    EXPECT_NEAR(tape_d.get(v0), 3.0, 1e-15);
    EXPECT_NEAR(tape_d.get(v1), 1.0, 1e-15);

    // EXPECT_NEAR(tape_d.get(v0), -std::sin(0.5) * std::cos(2.), 1e-13);
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
    auto result_adhoc = call_price(S, K, v, T);

    auto leaves_and_roots = TapeRootsAndLeafs(result_adhoc);
    leaves_and_roots.set(S, 100.0);
    leaves_and_roots.set(K, 102.0);
    leaves_and_roots.set(v, 0.15);
    leaves_and_roots.set(T, 0.5);

    auto intermediate_tape = evaluate_fwd_return_vals(leaves_and_roots);
    auto tape_d = TapeDerivatives(S, K, v, T, result_adhoc);

    tape_d.get(result_adhoc) = 1.0;

    evaluate_bwd(leaves_and_roots, intermediate_tape, tape_d);

    EXPECT_NEAR(tape_d.get(S), 0.33961663008862131, 1e-14);
    EXPECT_NEAR(tape_d.get(K), -0.38387614859866631, 1e-14);
    EXPECT_NEAR(tape_d.get(v), -30.580208040388474, 1e-14);
    EXPECT_NEAR(tape_d.get(T), -4.5870312060582705, 1e-14);
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

    auto leaves_and_roots = TapeRootsAndLeafs(res);
    leaves_and_roots.get(val1) = 1.0;
    leaves_and_roots.get(val2) = 2.0;
    leaves_and_roots.get(val3) = 2.0;
    leaves_and_roots.get(val4) = 3.0;
    leaves_and_roots.get(val5) = 4.0;
    auto intermediate_tape = evaluate_fwd_return_vals(leaves_and_roots);
    auto tape_d = TapeDerivatives(val1, val2, val3, res);

    tape_d.get(res) = 1.0;
    evaluate_bwd(leaves_and_roots, intermediate_tape, tape_d);

    EXPECT_NEAR(tape_d.get(val1), 412.44979516578081, 1e-14);
    EXPECT_NEAR(tape_d.get(val2), 1588.4738734117946, 1e-14);
    EXPECT_NEAR(tape_d.get(val3), 80, 1e-14);

    auto res2 = valprod3 * valsum;
    auto leaves_and_roots2 = TapeRootsAndLeafs(res2);
    leaves_and_roots2.get(val1) = 1.0;
    leaves_and_roots2.get(val2) = 2.0;
    leaves_and_roots2.get(val3) = 2.0;
    leaves_and_roots2.get(val4) = 3.0;
    leaves_and_roots2.get(val5) = 4.0;
    auto intermediate_tape2 = evaluate_fwd_return_vals(leaves_and_roots2);
    auto tape_d2 = TapeDerivatives(val1, val2, val3, res2);

    tape_d2.get(res2) = 1.0;
    evaluate_bwd(leaves_and_roots2, intermediate_tape2, tape_d2);

    EXPECT_NEAR(tape_d2.get(val1), 412.44979516578081, 1e-14);
    EXPECT_NEAR(tape_d2.get(val2), 1588.4738734117946, 1e-14);
    EXPECT_NEAR(tape_d2.get(val3), 80, 1e-14);
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

} // namespace adhoc2