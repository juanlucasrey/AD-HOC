#include <adhoc.hpp>
#include <constants_type.hpp>
#include <evaluate_bwd.hpp>
#include <evaluate_fwd.hpp>
#include <tape_size_bwd.hpp>
#include <tape_size_fwd.hpp>
#include <tape_static.hpp>
#include <utils.hpp>

#include <gtest/gtest.h>

namespace adhoc {

// TEST(AdHoc, DerivativeOne) {
//     adouble val1(1.);
//     auto derivatives = evaluate<decltype(val1)>(val1);
//     // static_assert(derivatives.empty());
// }

// TEST(AdHoc, DerivativeZero) {
//     adouble val1(1.);
//     adouble val2(1.);
//     auto derivatives = evaluate<decltype(val2)>(val1);
//     static_assert(derivatives.empty());
// }

TEST(AdHoc, DerivativeExp) {
    auto tape = CreateTapeInitial<1>();
    auto [val1] = tape.getalias();
    tape.set(val1, 1.0);
    auto temp = exp(val1);

    auto intermediate_tape = CreateTapeIntermediate(temp);
    evaluate_fwd(tape, intermediate_tape);

    auto derivatives = evaluate_bwd(tape, intermediate_tape, temp, 1.0, val1);

    static_assert(derivatives.size() == 1);
    EXPECT_NEAR(derivatives[0], std::exp(1.0), 1e-15);
}

TEST(AdHoc, DerivativeCos) {
    auto tape = CreateTapeInitial<1>();
    auto [val1] = tape.getalias();
    tape.set(val1, 1.0);
    auto temp = cos(val1);

    auto intermediate_tape = CreateTapeIntermediate(temp);
    evaluate_fwd(tape, intermediate_tape);
    auto derivatives = evaluate_bwd(tape, intermediate_tape, temp, 1.0, val1);
    static_assert(derivatives.size() == 1);
    EXPECT_NEAR(derivatives[0], -std::sin(1.0), 1e-15);
}

TEST(AdHoc, DerivativeAdd) {
    auto tape = CreateTapeInitial<2>();
    auto [val1, val2] = tape.getalias();
    tape.set(val1, 2.0);
    tape.set(val2, 3.0);
    auto temp = val1 + val2;

    auto intermediate_tape = CreateTapeIntermediate(temp);
    evaluate_fwd(tape, intermediate_tape);

    auto derivatives =
        evaluate_bwd(tape, intermediate_tape, temp, 1.0, val1, val2);

    static_assert(derivatives.size() == 2);
    EXPECT_NEAR(derivatives[0], 1.0, 1e-15);
    EXPECT_NEAR(derivatives[1], 1.0, 1e-15);
}

TEST(AdHoc, DerivativeMul) {
    auto tape = CreateTapeInitial<2>();
    auto [val1, val2] = tape.getalias();
    tape.set(val1, 2.0);
    tape.set(val2, 3.0);
    auto temp = val1 * val2;

    auto intermediate_tape = CreateTapeIntermediate(temp);
    evaluate_fwd(tape, intermediate_tape);

    auto derivatives =
        evaluate_bwd(tape, intermediate_tape, temp, 1.0, val1, val2);

    static_assert(derivatives.size() == 2);
    EXPECT_NEAR(derivatives[0], 3.0, 1e-15);
    EXPECT_NEAR(derivatives[1], 2.0, 1e-15);
}

TEST(AdHoc, DerivativeExpCos) {
    auto tape = CreateTapeInitial<1>();
    auto [val1] = tape.getalias();
    tape.set(val1, 1.0);
    auto valexp = exp(val1);
    auto valcos = cos(valexp);

    auto intermediate_tape = CreateTapeIntermediate(valcos);
    evaluate_fwd(tape, intermediate_tape);
    auto derivatives = evaluate_bwd(tape, intermediate_tape, valcos, 1.0, val1);
    static_assert(derivatives.size() == 1);

    // from sympy import *
    // x = Symbol('x')
    // f = cos(exp(x))
    // print(lambdify([x], f.diff(x))(1.0))

    EXPECT_NEAR(derivatives[0], -1.1166193174450132, 1e-15);
}

TEST(AdHoc, DerivativeDiv) {
    auto tape = CreateTapeInitial<2>();
    auto [val1, val2] = tape.getalias();
    tape.set(val1, 1.5);
    tape.set(val2, 2.5);
    auto temp = val1 / val2;

    auto intermediate_tape = CreateTapeIntermediate(temp);
    evaluate_fwd(tape, intermediate_tape);

    auto derivatives =
        evaluate_bwd(tape, intermediate_tape, temp, 1.0, val1, val2);

    static_assert(derivatives.size() == 2);
    // from sympy import *
    // x = Symbol('x')
    // y = Symbol('y')
    // f = x / y
    // print(lambdify([x, y], f.diff(x))(1.5, 2.5))
    // print(lambdify([x, y], f.diff(y))(1.5, 2.5))
    EXPECT_NEAR(derivatives[0], 0.4, 1e-10);
    EXPECT_NEAR(derivatives[1], -0.24, 1e-10);
}

TEST(AdHoc, DerivativeMulSame) {
    auto tape = CreateTapeInitial<1>();
    auto [val1] = tape.getalias();
    tape.set(val1, 1.5);
    auto temp = val1 * val1;

    auto intermediate_tape = CreateTapeIntermediate(temp);
    evaluate_fwd(tape, intermediate_tape);

    auto derivatives = evaluate_bwd(tape, intermediate_tape, temp, 1.0, val1);

    static_assert(derivatives.size() == 1);
    EXPECT_NEAR(derivatives[0], 3.0, 1e-15);
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

    // from sympy import *
    // x = Symbol('x')
    // y = Symbol('y')
    // z = Symbol('z')
    // f = ((x + exp(y)) + z) * ((x * y + 3) * (y * y * 4))
    // print(lambdify([x, y, z], f.diff(x))(1, 2, 2))
    // print(lambdify([x, y, z], f.diff(y))(1, 2, 2))
    // print(lambdify([x, y, z], f.diff(z))(1, 2, 2))

    static_assert(derivatives.size() == 3);
    EXPECT_NEAR(derivatives[0], 412.44979516578081, 1e-13);
    EXPECT_NEAR(derivatives[1], 1588.4738734117946, 1e-13);
    EXPECT_EQ(derivatives[2], 80);

    // TODO: const is needed! should not be
    constexpr std::size_t size =
        detail::tape_size(args<decltype(res) const>{}, args<>{},
                          args<decltype(val1) const, decltype(val2) const,
                               decltype(val3) const>{});
    static_assert(size == 5);

    auto res2 = valprod3 * valsum;
    constexpr std::size_t size2 =
        detail::tape_size(args<decltype(res2) const>{}, args<>{},
                          args<decltype(val1) const, decltype(val2) const,
                               decltype(val3) const>{});
    // same calculation, diferent size. order matters!
    static_assert(size2 == 4);

    auto intermediate_tape2 = CreateTapeIntermediate(res2);
    evaluate_fwd(tape, intermediate_tape2);

    auto derivatives2 =
        evaluate_bwd(tape, intermediate_tape2, res2, 1.0, val1, val2, val3);
    static_assert(derivatives2.size() == 3);
    EXPECT_NEAR(derivatives2[0], 412.44979516578081, 1e-13);
    EXPECT_NEAR(derivatives2[1], 1588.4738734117946, 1e-13);
    EXPECT_EQ(derivatives2[2], 80);
}

// TEST(adhoc2, Derivative2nd) {
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