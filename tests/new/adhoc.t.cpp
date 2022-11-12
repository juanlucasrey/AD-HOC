#include <new/adhoc.hpp>
#include <new/constants_type.hpp>
#include <new/evaluate_bwd.hpp>
#include <new/evaluate_fwd.hpp>
#include <new/tape_size_bwd.hpp>
#include <new/tape_size_fwd.hpp>
#include <new/tape_static.hpp>
#include <new/utils.hpp>

#include <gtest/gtest.h>

namespace adhoc {

TEST(AdHoc, UnivariateExp) {

    auto tape = CreateTapeInitial<1>();
    auto [S] = tape.getalias();
    tape.set(S, 100.0);

    auto result_adhoc = exp(S);

    auto intermediate_tape = CreateTapeIntermediate(result_adhoc);
    evaluate_fwd(tape, intermediate_tape);

    double result2 = intermediate_tape.get(result_adhoc);
    const double exp100 = std::exp(100.0);
    EXPECT_EQ(result2, exp100);

    auto derivatives =
        evaluate_bwd(tape, intermediate_tape, result_adhoc, 1.0, S);

    static_assert(derivatives.size() == 1);
    EXPECT_EQ(derivatives[0], exp100);
}

TEST(AdHoc, BivariateMult) {

    auto tape = CreateTapeInitial<2>();
    auto [v1, v2] = tape.getalias();
    tape.set(v1, 100.0);
    tape.set(v2, 10.0);

    auto result_adhoc = v1 * v2;

    auto intermediate_tape = CreateTapeIntermediate(result_adhoc);
    evaluate_fwd(tape, intermediate_tape);

    double result2 = intermediate_tape.get(result_adhoc);
    EXPECT_EQ(result2, 1000.0);

    auto derivatives =
        evaluate_bwd(tape, intermediate_tape, result_adhoc, 1.0, v1, v2);

    static_assert(derivatives.size() == 2);
    EXPECT_EQ(derivatives[0], 10.0);
    EXPECT_EQ(derivatives[1], 100.0);
}

TEST(AdHoc, ComplexEvaluate) {
    auto tape = CreateTapeInitial<5>();
    auto [val1, val2, val3, val4, val5] = tape.getalias();
    tape.set(val1, 1.0);
    tape.set(val2, 2.0);
    tape.set(val3, 2.0);
    tape.set(val4, 3.0);
    tape.set(val5, 4.0);
    auto exp2 = exp(val2);
    auto tmp1 = val1 + exp2;
    auto valsum = tmp1 + val3;
    auto val12 = val1 * val2;
    auto valprod = val12 + val4;
    auto val22 = val2 * val2;
    auto valprod2 = val22 * val5;
    auto valprod3 = valprod * valprod2;
    auto res = valsum * valprod3;

    auto intermediate_tape = CreateTapeIntermediate(res);
    evaluate_fwd(tape, intermediate_tape);

    auto derivatives =
        evaluate_bwd(tape, intermediate_tape, res, 1.0, val1, val2, val3);

    static_assert(derivatives.size() == 3);
    EXPECT_NEAR(derivatives[0], 412.44979516578081, 1e-13);
    EXPECT_NEAR(derivatives[1], 1588.4738734117946, 1e-13);
    EXPECT_EQ(derivatives[2], 80);
}

TEST(AdHoc, SimpleEvaluate) {
    auto tape = CreateTapeInitial<2>();
    auto [val1, val2] = tape.getalias();
    tape.set(val1, 1.0);
    tape.set(val2, 2.0);
    auto valprod3 = (val1 * val2) + (val1);

    auto intermediate_tape = CreateTapeIntermediate(valprod3);
    evaluate_fwd(tape, intermediate_tape);

    auto derivatives =
        evaluate_bwd(tape, intermediate_tape, valprod3, 1.0, val1, val2);

    static_assert(derivatives.size() == 2);
    EXPECT_NEAR(derivatives[0], 3.0, 1e-15);
    EXPECT_NEAR(derivatives[1], 1.0, 1e-15);
}

// TEST(adhoc2, derivativeone) {
//     adouble val1(1.);
//     auto derivatives = evaluate<decltype(val1)>(val1);
//     // static_assert(derivatives.empty());
// }

// TEST(adhoc2, derivativezero) {
//     adouble val1(1.);
//     adouble val2(1.);
//     auto derivatives = evaluate<decltype(val2)>(val1);
//     static_assert(derivatives.empty());
// }

// TEST(adhoc2, derivativeexp) {
//     adouble val1(1.);
//     auto temp = exp(val1);
//     auto derivatives = evaluate<decltype(val1)>(temp);
//     static_assert(derivatives.size() == 1);
// }

// TEST(adhoc2, derivativecos) {
//     adouble val1(1.);
//     auto temp = cos(val1);
//     auto derivatives = evaluate<decltype(val1)>(temp);
//     static_assert(derivatives.size() == 1);
// }

// TEST(adhoc2, derivativeadd) {
//     adouble val1(1.);
//     adouble val2(1.);
//     auto temp = val1 + val2;
//     auto derivatives = evaluate<decltype(val1), decltype(val2)>(temp);
//     static_assert(derivatives.size() == 2);
// }

// TEST(adhoc2, derivativemul) {
//     adouble val1(1.);
//     adouble val2(1.);
//     auto temp = val1 * val2;
//     auto derivatives = evaluate<decltype(val1), decltype(val2)>(temp);
//     static_assert(derivatives.size() == 2);
// }

// TEST(adhoc2, derivative) {
//     adouble val1(1.);
//     auto valexp = exp(val1);
//     auto valcos = cos(valexp);
//     constexpr std::size_t size = tape_size(valcos, val1);
//     static_assert(size == 1);

//     double derivative = valcos.d();
//     derivative *= valexp.d();
//     // www.wolframalpha.com/input?i=d%2Fdx+cos%28exp%28x%29%29+%7C+x%3D1
//     EXPECT_NEAR(derivative, -1.116619317445013, 1e-10);
// }

// TEST(adhoc2, DerivativeUniv) {
//     adouble x(1.);
//     auto result = evaluate<decltype(x)>(cos(exp(x)));

//     static_assert(result.size() == 1);
//     // from sympy import *
//     // x = Symbol('x')
//     // f = cos(exp(x))
//     // dfdx = f.diff(x)
//     // dfdx = lambdify(x, dfdx)
//     // print(dfdx(1))
//     EXPECT_NEAR(result[0], -1.1166193174450132, 1e-10);
// }

// auto makeTemp() {
//     adouble x(1.);
//     auto valexp = exp(x);
//     auto valcos = cos(valexp);
//     return std::tuple<decltype(valcos), decltype(x)>{valcos, x};
// }

// TEST(adhoc2, DerivativeUnivScope) {
//     // this shows why we need deep copies. thisexample would not work on
//     release
//     // if we would use references
//     auto temp = makeTemp();
//     auto result = evaluate<std::tuple_element<1, decltype(temp)>::type>(
//         std::get<0>(temp));

//     static_assert(result.size() == 1);
//     EXPECT_NEAR(result[0], -1.1166193174450132, 1e-10);
// }

// TEST(adhoc2, DerivativeMult) {
//     adouble x(1.5);
//     adouble y(2.5);
//     auto result = evaluate<decltype(x), decltype(y)>(x * y);

//     static_assert(result.size() == 2);
//     // from sympy import *
//     // x = Symbol('x')
//     // y = Symbol('y')
//     // f = x * y
//     // dfdx = f.diff(x)
//     // dfdy = f.diff(y)
//     // dfdx = lambdify([x, y], dfdx)
//     // dfdy = lambdify([x, y], dfdy)
//     // print(dfdx(1.5, 2.5))
//     // print(dfdy(1.5, 2.5))
//     EXPECT_NEAR(result[0], 2.5, 1e-10);
//     EXPECT_NEAR(result[1], 1.5, 1e-10);
//     std::cout << result[0] << std::endl;
//     std::cout << result[1] << std::endl;
// }

// TEST(adhoc2, DerivativeDiv) {
//     adouble x(1.5);
//     adouble y(2.5);
//     auto result = evaluate<decltype(x), decltype(y)>(x / y);

//     static_assert(result.size() == 2);
//     // from sympy import *
//     // x = Symbol('x')
//     // y = Symbol('y')
//     // f = x / y
//     // print(lambdify([x, y], f.diff(x))(1.5, 2.5))
//     // print(lambdify([x, y], f.diff(y))(1.5, 2.5))
//     EXPECT_NEAR(result[0], 0.4, 1e-10);
//     EXPECT_NEAR(result[1], -0.24, 1e-10);
// }

// TEST(adhoc2, DerivativeMultSame) {
//     adouble x(1.5);
//     auto result = evaluate<decltype(x)>(x * x);

//     static_assert(result.size() == 1);
//     // from sympy import *
//     // x = Symbol('x')
//     // y = Symbol('y')
//     // f = x * y
//     // dfdx = f.diff(x)
//     // dfdy = f.diff(y)
//     // dfdx = lambdify([x, y], dfdx)
//     // dfdy = lambdify([x, y], dfdy)
//     // print(dfdx(1.5, 2.5))
//     // print(dfdy(1.5, 2.5))
//     EXPECT_NEAR(result[0], 3.0, 1e-10);

//     auto f = x * x;
//     constexpr std::size_t size = tape_size(f, x);
//     static_assert(size == 1);
// }

// TEST(adhoc2, derivcomplexdexp2) {
//     adouble x(1.);
//     adouble y(2.);
//     adouble z(2.);
//     auto r1 = ((x + exp(y)) + z);
//     auto r2 = ((x * y + adouble(3)) * (y * y * adouble(4)));
//     auto res = r1 * r2;

//     constexpr std::size_t size = tape_size(res, x, y, z);
//     static_assert(size == 5);
//     auto result = evaluate<decltype(x), decltype(y), decltype(z)>(res);
//     static_assert(result.size() == 3);

//     // from sympy import *
//     // x = Symbol('x')
//     // y = Symbol('y')
//     // z = Symbol('z')
//     // f = ((x + exp(y)) + z) * ((x * y + 3) * (y * y * 4))
//     // print(lambdify([x, y, z], f.diff(x))(1, 2, 2))
//     // print(lambdify([x, y, z], f.diff(y))(1, 2, 2))
//     // print(lambdify([x, y, z], f.diff(z))(1, 2, 2))
//     EXPECT_NEAR(result[0], 412.4497951657808, 1e-9);
//     EXPECT_NEAR(result[1], 1588.4738734117946, 1e-9);
//     EXPECT_NEAR(result[2], 80., 1e-9);

//     auto res2 = r2 * r1;
//     constexpr std::size_t size2 = tape_size(res2, x, y, z);
//     static_assert(size2 == 4);
//     auto result2 = evaluate<decltype(x), decltype(y),
//     decltype(z)>(res2, 1.0); static_assert(result2.size() == 3);
//     EXPECT_NEAR(result2[0], 412.4497951657808, 1e-9);
//     EXPECT_NEAR(result2[1], 1588.4738734117946, 1e-9);
//     EXPECT_NEAR(result2[2], 80., 1e-9);
// }

// template <class D> inline auto cdf_n(D const &x) {
//     using std::erfc;
//     using namespace constants;
//     constexpr double minus_one_over_root_two =
//         -1.0 / constexpression::sqrt(2.0);
//     return CD<encode(0.5)>() * erfc(x *
//     CD<encode(minus_one_over_root_two)>());
//     // return CD<0.5>() * erfc(x * CD<minus_one_over_root_two>());
// }

// template <class I1, class I2, class I3, class I4>
// auto call_price(const I1 &S, const I2 &K, const I3 &v, const I4 &T) {
//     using std::log;
//     using std::sqrt;
//     using namespace constants;
//     auto totalvol = v * sqrt(T);
//     auto d1 = log(S / K) / totalvol + totalvol * CD<encode(0.5)>();
//     auto d2 = d1 + totalvol;
//     return S * cdf_n(d1) - K * cdf_n(d2);
// }

// TEST(adhoc2, BlackScholes) {
//     double S = 100.0;
//     double K = 102.0;
//     double v = 0.15;
//     double T = 0.5;
//     double result = call_price(S, K, v, T);
//     double epsilon = 1e-5;
//     double resultp1 = call_price(S + epsilon, K, v, T);
//     double resultm1 = call_price(S - epsilon, K, v, T);
//     std::array<double, 4> fd{};
//     fd[0] = (resultp1 - resultm1) / (2.0 * epsilon);

//     double resultp2 = call_price(S, K + epsilon, v, T);
//     double resultm2 = call_price(S, K - epsilon, v, T);
//     fd[1] = (resultp2 - resultm2) / (2.0 * epsilon);

//     double resultp3 = call_price(S, K, v + epsilon, T);
//     double resultm3 = call_price(S, K, v - epsilon, T);
//     fd[2] = (resultp3 - resultm3) / (2.0 * epsilon);

//     double resultp4 = call_price(S, K, v, T + epsilon);
//     double resultm4 = call_price(S, K, v, T - epsilon);
//     fd[3] = (resultp4 - resultm4) / (2.0 * epsilon);

//     adhoc<ID> aS(S);
//     adhoc<ID> aK(K);
//     adhoc<ID> av(v);
//     adhoc<ID> aT(T);
//     auto result2 = call_price(aS, aK, av, aT);
//     std::cout.precision(std::numeric_limits<double>::max_digits10);
//     auto derivatives =
//         evaluate<decltype(aS), decltype(aK), decltype(av), decltype(aT)>(
//             result2);
//     EXPECT_EQ(result, result2.v());
//     EXPECT_NEAR(fd[0], derivatives[0], 1e-8);
//     EXPECT_NEAR(fd[1], derivatives[1], 1e-8);
//     EXPECT_NEAR(fd[2], derivatives[2], 1e-8);
//     EXPECT_NEAR(fd[3], derivatives[3], 1e-8);
//     std::array<double, 5> tape{};
//     tape[0] = 1;
//     tape[1] = tape[0] * 1.0;
//     tape[0] *= 1.0;
//     tape[2] = tape[0] * 0.44683257413550148;
//     tape[0] *= 100;
//     tape[0] *= 0.5;
//     tape[0] *= (-std::exp(-0.094517515307865022 * 0.094517515307865022) *
//                 1.1283791670955126);
//     tape[0] *= -0.70710678118654757;
//     tape[3] = tape[1] * -102;
//     tape[1] *= 0.48898981923286255;
//     tape[1] *= -1;
//     tape[3] *= 0.5;
//     tape[3] *= (-std::exp(-0.019517515307865021 * 0.019517515307865021) *
//                 1.1283791670955126);
//     tape[3] *= -0.70710678118654757;
//     tape[0] += tape[3] * 1.0;
//     tape[3] *= 1.0;
//     tape[4] = tape[0] * 1.0;
//     tape[0] *= 1.0;
//     tape[3] += tape[0] * (0.1867009606191804 / 0.10606601717798213);
//     tape[0] *= (1.0 / 0.10606601717798213);
//     tape[0] *= (1.0 / 0.98039215686274506);
//     tape[2] += tape[0] * (1.0 / 102);
//     tape[1] += tape[0] * (-0.98039215686274506 / 102);
//     tape[3] += tape[4] * 0.5;
//     tape[4] = tape[3] * 0.70710678118654757;
//     tape[3] *= 0.14999999999999999;
//     tape[3] *= (0.5 * 0.70710678118654757 / 0.5);
//     std::cout << "done" << std::endl;
// }

// TEST(adhoc2, BlackScholes2) {
//     TapeInput<4> tape;
//     auto [S, K, v, T] = tape.getalias();
//     tape.set(S, 100.0);
//     tape.set(K, 102.0);
//     tape.set(v, 0.15);
//     tape.set(T, 0.5);
//     auto result = call_price(S, K, v, T);
//     auto derivatives =
//         evaluate<decltype(S), decltype(K), decltype(v), decltype(T)>(result);
// }

// // TEST(adhoc2, Derivative2nd) {
// //     adouble x(2);
// //     adouble y(3);
// //     adouble z(5);
// //     auto x2 = x * x;
// //     auto y2 = y * y;
// //     auto z2 = z * z;
// //     auto temp = x2 * y2;
// //     auto temp2 = temp * z2;
// //     auto res = temp * temp2;

// //     // f1(temp, temp2) = temp * temp2
// //     // df1dtemp = dtemp/dtemp * temp2 + temp * dtemp2/dtemp
// //     // df1dtemp_2 = dtemp/dtemp_2 * temp2 + temp * dtemp2/dtemp_2 +
// //     dtemp/dtemp * dtemp2/dtemp double dresdtemp = temp2.v(); double
// //     dresdtemp2 = temp.v(); double d2_res_d_temp_2 = 0; double
// //     d2_res_d_temp2_2 = 0; double d2_res_dtemp2_dtemp = 1;

// //     // f2(temp, z2) = temp * z2
// //     // f2(temp, z2) = temp * (temp * z2) = f1(temp, f2(temp, z2))
// //     // dfdtemp_2 = dtemp/dtemp_2 * (temp * z2) + temp * d(temp *
// z2)/dtemp_2
// //     + dtemp/dtemp * d(temp * z2)/dtemp dresdtemp += dresdtemp2 * z2.v();
// //     double d2resdz2 = dresdtemp2 * temp.v();
// //     d2_res_d_temp_2 += 0;
// //     double d2_res_d_z2_2 = 0;
// //     double d2_res_dz2_dtemp = 0;
// // }

// // from sympy import *
// // x = Symbol('x')
// // y = Symbol('y')
// // z = Symbol('z')
// // temp = ((x * x) * (y * y))
// // temp2 = temp * (z * z)
// // f = temp * temp2
// // print(f)
// // dfdx = f.diff(x)
// // print(dfdx)
// // d2fdx2 = dfdx.diff(x)
// // print(d2fdx2)

// // dfdx = lambdify([x, y, z], dfdx)
// // d2fdx2 = lambdify([x, y, z], d2fdx2)
// // # dfdy = lambdify([x, y], dfdy)
// // # print(dfdx(1.5, 2.5))
// // print(dfdx(2, 3, 5))
// // print(d2fdx2(2, 3, 5))

} // namespace adhoc