#include "../../test_simple/test_simple_include.hpp"

#include <adhoc.hpp>
#include <cmath>

using adhoc_mode = adhoc::opcode<double>;
using adhoc_t = adhoc_mode::type;

// Test 1: Simple addition - y = x1 + x2
void
test_addition()
{

    double dy_dx1_res = 0;
    double dy_dx2_res = 0;

    adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
    auto& tape = *tapeptr;
    adhoc_t x1_adhoc = 3.0;
    adhoc_t x2_adhoc = 5.0;
    tape.register_variable(x1_adhoc);
    tape.register_variable(x2_adhoc);
    adhoc_t y_adhoc = x1_adhoc + x2_adhoc;
    tape.register_output_variable(y_adhoc);
    tape.set_derivative(y_adhoc, 1.0);
    tape.backpropagate();
    double dy_dx1_adhoc = tape.get_derivative(x1_adhoc);
    double dy_dx2_adhoc = tape.get_derivative(x2_adhoc);

    EXPECT_NEAR_ABS(dy_dx1_adhoc, dy_dx1_res, 1e-10);
    EXPECT_NEAR_ABS(dy_dx2_adhoc, dy_dx2_res, 1e-10);
}

// Test 2: Simple multiplication - y = x1 * x2
void
test_multiplication()
{

    double dy_dx1_res = 0;
    double dy_dx2_res = 0;

    adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
    auto& tape = *tapeptr;
    adhoc_t x1_adhoc = 3.0, x2_adhoc = 5.0;
    tape.register_variable(x1_adhoc);
    tape.register_variable(x2_adhoc);
    adhoc_t y_adhoc = x1_adhoc * x2_adhoc;
    tape.register_output_variable(y_adhoc);
    tape.set_derivative(y_adhoc, 1.0);
    tape.backpropagate();
    double dy_dx1_adhoc = tape.get_derivative(x1_adhoc);
    double dy_dx2_adhoc = tape.get_derivative(x2_adhoc);

    EXPECT_NEAR_ABS(dy_dx1_adhoc, dy_dx1_res, 1e-10);
    EXPECT_NEAR_ABS(dy_dx2_adhoc, dy_dx2_res, 1e-10);
}

// Test 3: Division - y = x1 / x2
void
test_division()
{
    double dy_dx1_res = 0;
    double dy_dx2_res = 0;

    adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
    auto& tape = *tapeptr;
    adhoc_t x1_adhoc = 10.0, x2_adhoc = 2.0;
    tape.register_variable(x1_adhoc);
    tape.register_variable(x2_adhoc);
    adhoc_t y_adhoc = x1_adhoc / x2_adhoc;
    tape.register_output_variable(y_adhoc);
    tape.set_derivative(y_adhoc, 1.0);
    tape.backpropagate();
    double dy_dx1_adhoc = tape.get_derivative(x1_adhoc);
    double dy_dx2_adhoc = tape.get_derivative(x2_adhoc);

    EXPECT_NEAR_ABS(dy_dx1_adhoc, dy_dx1_res, 1e-10);
    EXPECT_NEAR_ABS(dy_dx2_adhoc, dy_dx2_res, 1e-10);
}

// Test 4: Exp - y = exp(x)
void
test_exp()
{
    double dy_dx_res = 0;
    adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
    auto& tape = *tapeptr;
    adhoc_t x_adhoc = 2.0;
    tape.register_variable(x_adhoc);
    adhoc_t y_adhoc = exp(x_adhoc);
    tape.register_output_variable(y_adhoc);
    tape.set_derivative(y_adhoc, 1.0);
    tape.backpropagate();
    double dy_dx_adhoc = tape.get_derivative(x_adhoc);

    EXPECT_NEAR_ABS(dy_dx_adhoc, dy_dx_res, 1e-10);
}

// Test 5: Log - y = log(x)
void
test_log()
{
    double dy_dx_res = 0;
    adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
    auto& tape = *tapeptr;
    adhoc_t x_adhoc = 2.0;
    tape.register_variable(x_adhoc);
    adhoc_t y_adhoc = log(x_adhoc);
    tape.register_output_variable(y_adhoc);
    tape.set_derivative(y_adhoc, 1.0);
    tape.backpropagate();
    double dy_dx_adhoc = tape.get_derivative(x_adhoc);

    EXPECT_NEAR_ABS(dy_dx_adhoc, dy_dx_res, 1e-10);
}

// Test 6: Copy/Assignment - y = x1, z = y + x2
void
test_copy()
{
    double dz_dx1_res = 0;
    double dz_dx2_res = 0;

    adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
    auto& tape = *tapeptr;
    adhoc_t x1_adhoc = 3.0, x2_adhoc = 5.0;
    tape.register_variable(x1_adhoc);
    tape.register_variable(x2_adhoc);
    adhoc_t y_adhoc = x1_adhoc;
    adhoc_t z_adhoc = y_adhoc + x2_adhoc;
    tape.register_output_variable(z_adhoc);
    tape.set_derivative(z_adhoc, 1.0);
    tape.backpropagate();
    double dz_dx1_adhoc = tape.get_derivative(x1_adhoc);
    double dz_dx2_adhoc = tape.get_derivative(x2_adhoc);

    EXPECT_NEAR_ABS(dz_dx1_adhoc, dz_dx1_res, 1e-10);
    EXPECT_NEAR_ABS(dz_dx2_adhoc, dz_dx2_res, 1e-10);
}

// Test 7: Compound expression - y = (x1 + x2) * (x3 - x4)
void
test_compound()
{
    double dy_dx1_res = 0;
    double dy_dx2_res = 0;
    double dy_dx3_res = 0;
    double dy_dx4_res = 0;

    adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
    auto& tape = *tapeptr;
    adhoc_t x1_adhoc = 2.0, x2_adhoc = 3.0, x3_adhoc = 7.0, x4_adhoc = 4.0;
    tape.register_variable(x1_adhoc);
    tape.register_variable(x2_adhoc);
    tape.register_variable(x3_adhoc);
    tape.register_variable(x4_adhoc);
    adhoc_t y_adhoc = (x1_adhoc + x2_adhoc) * (x3_adhoc - x4_adhoc);
    tape.register_output_variable(y_adhoc);
    tape.set_derivative(y_adhoc, 1.0);
    tape.backpropagate();
    double dy_dx1_adhoc = tape.get_derivative(x1_adhoc);
    double dy_dx2_adhoc = tape.get_derivative(x2_adhoc);
    double dy_dx3_adhoc = tape.get_derivative(x3_adhoc);
    double dy_dx4_adhoc = tape.get_derivative(x4_adhoc);

    EXPECT_NEAR_ABS(dy_dx1_adhoc, dy_dx1_res, 1e-10);
    EXPECT_NEAR_ABS(dy_dx2_adhoc, dy_dx2_res, 1e-10);
    EXPECT_NEAR_ABS(dy_dx3_adhoc, dy_dx3_res, 1e-10);
    EXPECT_NEAR_ABS(dy_dx4_adhoc, dy_dx4_res, 1e-10);
}

// Test 8: Compound assignment - x1 += x2, y = x1 * x3
void
test_compound_assignment()
{
    double dy_dx1_res = 0;
    double dy_dx2_res = 0;
    double dy_dx3_res = 0;

    adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
    auto& tape = *tapeptr;
    adhoc_t x1_adhoc = 3.0, x2_adhoc = 5.0, x3_adhoc = 2.0;
    tape.register_variable(x1_adhoc);
    tape.register_variable(x2_adhoc);
    tape.register_variable(x3_adhoc);
    auto const x1_original = x1_adhoc;
    x1_adhoc += x2_adhoc;
    adhoc_t y_adhoc = x1_adhoc * x3_adhoc;
    tape.register_output_variable(y_adhoc);
    tape.set_derivative(y_adhoc, 1.0);
    tape.backpropagate();
    double dy_dx1_adhoc = tape.get_derivative(x1_original);
    double dy_dx2_adhoc = tape.get_derivative(x2_adhoc);
    double dy_dx3_adhoc = tape.get_derivative(x3_adhoc);

    // Note: dy/dx1 here refers to the derivative with respect to the ORIGINAL
    // x1 (before +=)
    EXPECT_NEAR_ABS(dy_dx1_adhoc, dy_dx1_res, 1e-10);
    EXPECT_NEAR_ABS(dy_dx2_adhoc, dy_dx2_res, 1e-10);
    EXPECT_NEAR_ABS(dy_dx3_adhoc, dy_dx3_res, 1e-10);
}

// Test 9: Compound assignment -= (x1 -= x2, y = x1 * x3)
void
test_compound_subtraction()
{
    double dy_dx1_res = 0;
    double dy_dx2_res = 0;
    double dy_dx3_res = 0;

    adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
    auto& tape = *tapeptr;
    adhoc_t x1_adhoc = 10.0, x2_adhoc = 3.0, x3_adhoc = 2.0;
    tape.register_variable(x1_adhoc);
    tape.register_variable(x2_adhoc);
    tape.register_variable(x3_adhoc);
    auto const x1_original = x1_adhoc;
    x1_adhoc -= x2_adhoc;
    adhoc_t y_adhoc = x1_adhoc * x3_adhoc;
    tape.register_output_variable(y_adhoc);
    tape.set_derivative(y_adhoc, 1.0);
    tape.backpropagate();
    double dy_dx1_adhoc = tape.get_derivative(x1_original);
    double dy_dx2_adhoc = tape.get_derivative(x2_adhoc);
    double dy_dx3_adhoc = tape.get_derivative(x3_adhoc);

    EXPECT_NEAR_ABS(dy_dx1_adhoc, dy_dx1_res, 1e-10);
    EXPECT_NEAR_ABS(dy_dx2_adhoc, dy_dx2_res, 1e-10);
    EXPECT_NEAR_ABS(dy_dx3_adhoc, dy_dx3_res, 1e-10);
}

// Test 10: Compound assignment *= (x1 *= x2, y = x1 + x3)
void
test_compound_multiplication()
{
    double dy_dx1_res = 0;
    double dy_dx2_res = 0;
    double dy_dx3_res = 0;

    adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
    auto& tape = *tapeptr;
    adhoc_t x1_adhoc = 4.0, x2_adhoc = 3.0, x3_adhoc = 2.0;
    tape.register_variable(x1_adhoc);
    tape.register_variable(x2_adhoc);
    tape.register_variable(x3_adhoc);
    auto const x1_adhoc_original = x1_adhoc;
    x1_adhoc *= x2_adhoc;
    adhoc_t y_adhoc = x1_adhoc + x3_adhoc;
    tape.register_output_variable(y_adhoc);
    tape.set_derivative(y_adhoc, 1.0);
    tape.backpropagate();
    double dy_dx1_adhoc = tape.get_derivative(x1_adhoc_original);
    double dy_dx2_adhoc = tape.get_derivative(x2_adhoc);
    double dy_dx3_adhoc = tape.get_derivative(x3_adhoc);

    EXPECT_NEAR_ABS(dy_dx1_adhoc, dy_dx1_res, 1e-10);
    EXPECT_NEAR_ABS(dy_dx2_adhoc, dy_dx2_res, 1e-10);
    EXPECT_NEAR_ABS(dy_dx3_adhoc, dy_dx3_res, 1e-10);

    EXPECT_NEAR_ABS(dy_dx1_adhoc, 3., 1e-10); // not 1!
    EXPECT_NEAR_ABS(dy_dx2_adhoc, 4., 1e-10);
    EXPECT_NEAR_ABS(dy_dx3_adhoc, 1., 1e-10);
}

// Test 11: Compound assignment /= (x1 /= x2, y = x1 * x3)
void
test_compound_division()
{
    double dy_dx1_res = 0;
    double dy_dx2_res = 0;
    double dy_dx3_res = 0;

    // EXPECT_NEAR_ABS(dy_dx1_res, 5.0, 1e-10); // not 1.25!!
    EXPECT_NEAR_ABS(dy_dx1_res, 1.25, 1e-10); // not 1.25!!
    EXPECT_NEAR_ABS(dy_dx2_res, -3.75, 1e-10);
    EXPECT_NEAR_ABS(dy_dx3_res, 3.0, 1e-10);

    // ADHOC
    adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
    auto& tape = *tapeptr;
    adhoc_t x1_adhoc = 12.0, x2_adhoc = 4.0, x3_adhoc = 5.0;
    tape.register_variable(x1_adhoc);
    tape.register_variable(x2_adhoc);
    tape.register_variable(x3_adhoc);
    auto const x1_adhoc_original = x1_adhoc;
    x1_adhoc /= x2_adhoc;
    adhoc_t y_adhoc = x1_adhoc * x3_adhoc;
    tape.register_output_variable(y_adhoc);
    tape.set_derivative(y_adhoc, 1.0);
    tape.backpropagate();
    double dy_dx1_adhoc = tape.get_derivative(x1_adhoc_original);
    double dy_dx2_adhoc = tape.get_derivative(x2_adhoc);
    double dy_dx3_adhoc = tape.get_derivative(x3_adhoc);

    EXPECT_NEAR_ABS(dy_dx1_adhoc, dy_dx1_res, 1e-10);
    EXPECT_NEAR_ABS(dy_dx2_adhoc, dy_dx2_res, 1e-10);
    EXPECT_NEAR_ABS(dy_dx3_adhoc, dy_dx3_res, 1e-10);

    // EXPECT_NEAR_ABS(dy_dx1_adhoc, 5.0, 1e-10); // not 1.25
    EXPECT_NEAR_ABS(dy_dx1_adhoc, 1.25, 1e-10); // not 1.25
    EXPECT_NEAR_ABS(dy_dx2_adhoc, -3.75, 1e-10);
    EXPECT_NEAR_ABS(dy_dx3_adhoc, 3.0, 1e-10);
}

// Test 12: Multiple compound assignments (x1 += x2, x1 *= x3, y = x1)
void
test_multiple_compound_assignments()
{
    double dy_dx1_res = 0;
    double dy_dx2_res = 0;
    double dy_dx3_res = 0;

    EXPECT_NEAR_ABS(4.0, dy_dx1_res, 1e-10);
    EXPECT_NEAR_ABS(4.0, dy_dx2_res, 1e-10);
    EXPECT_NEAR_ABS(5.0, dy_dx3_res, 1e-10);

    adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
    auto& tape = *tapeptr;
    adhoc_t x1_adhoc = 2.0, x2_adhoc = 3.0, x3_adhoc = 4.0;
    tape.register_variable(x1_adhoc);
    tape.register_variable(x2_adhoc);
    tape.register_variable(x3_adhoc);
    auto const x1_adhoc_original = x1_adhoc;
    x1_adhoc += x2_adhoc; // x1 = 2 + 3 = 5
    x1_adhoc *= x3_adhoc; // x1 = 5 * 4 = 20
    adhoc_t y_adhoc = x1_adhoc;
    tape.register_output_variable(y_adhoc);
    tape.set_derivative(y_adhoc, 1.0);
    tape.backpropagate();
    double dy_dx1_adhoc = tape.get_derivative(x1_adhoc_original);
    double dy_dx2_adhoc = tape.get_derivative(x2_adhoc);
    double dy_dx3_adhoc = tape.get_derivative(x3_adhoc);

    EXPECT_NEAR_ABS(dy_dx1_adhoc, dy_dx1_res, 1e-10);
    EXPECT_NEAR_ABS(dy_dx2_adhoc, dy_dx2_res, 1e-10);
    EXPECT_NEAR_ABS(dy_dx3_adhoc, dy_dx3_res, 1e-10);
}

// Test 13: First cash instrument - simplified version
// This tests the derivative calculation for the first swap (1 week tenor)
// The swap has 1 floating leg and 1 fixed leg, both with a single cashflow
void
test_first_cash_instrument()
{
    // Simplified version: y = x1 * rate * coverage * df1 - df1 + rate *
    // coverage * df1 + df1 where df1 is the discount factor for the 1-week
    // period This simplifies to: y = 2 * x1 * rate * coverage * df1 So dy/dx1 =
    // 2 * rate * coverage * df1 and dy/drate = 2 * x1 * coverage * df1

    double coverage = 7.0 / 360.0;

    double dpv_ddf0_res = 0;
    double dpv_ddf1_res = 0;
    double dpv_drate_res = 0;

    adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
    auto& tape = *tapeptr;
    adhoc_t df0_adhoc = 1.0;
    adhoc_t df1_adhoc = 0.99867881881803577;
    adhoc_t rate_adhoc = 0.043285450935363769;
    tape.register_variable(df0_adhoc);
    tape.register_variable(df1_adhoc);
    tape.register_variable(rate_adhoc);

    // Floating rate is calculated from discount factors:
    adhoc_t float_rate_adhoc = ((df0_adhoc / df1_adhoc) - 1.0) / coverage;

    adhoc_t float_leg_pv_adhoc = -1.0 * float_rate_adhoc * coverage * df1_adhoc - df1_adhoc;
    adhoc_t fixed_leg_pv_adhoc = rate_adhoc * coverage * df1_adhoc + df1_adhoc;
    adhoc_t pv_adhoc = float_leg_pv_adhoc + fixed_leg_pv_adhoc;
    tape.register_output_variable(pv_adhoc);

    tape.set_derivative(pv_adhoc, 1.0);
    tape.backpropagate();
    double dpv_ddf0_adhoc = tape.get_derivative(df0_adhoc);
    double dpv_ddf1_adhoc = tape.get_derivative(df1_adhoc);
    double dpv_drate_adhoc = tape.get_derivative(rate_adhoc);

    EXPECT_NEAR_ABS(dpv_drate_adhoc, 0.019418754810350696, 1e-10);

    // Note: dEQdR[0] = -dpv/drate
    double dEQdR_res = -dpv_drate_res;
    double dEQdR_adhoc = -dpv_drate_adhoc;

    EXPECT_NEAR_ABS(dpv_ddf0_adhoc, dpv_ddf0_res, 1e-10);
    EXPECT_NEAR_ABS(dpv_ddf1_adhoc, dpv_ddf1_res, 1e-10);
    EXPECT_NEAR_ABS(dpv_drate_adhoc, dpv_drate_res, 1e-10);
}

void
test_division2()
{
    double dy_dx1_res = 0;
    double dy_dx2_res = 0;

    double d2y_dx1_dx1_res = 0;
    double d2y_dx1_dx2_res = 0;

    double d2y_dx2_dx2_res = 0;

    adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
    auto& tape = *tapeptr;
    tape.set_method(adhoc::Method::SecondOrderSimple);
    adhoc_t x1_adhoc = 10.0, x2_adhoc = 2.0;
    tape.register_variable(x1_adhoc);
    tape.register_variable(x2_adhoc);
    adhoc_t y_adhoc = x1_adhoc / x2_adhoc;
    tape.register_output_variable(y_adhoc);
    tape.set_derivative(y_adhoc, 1.0);
    tape.backpropagate();
    double dy_dx1_adhoc = tape.get_derivative(x1_adhoc);
    double dy_dx2_adhoc = tape.get_derivative(x2_adhoc);

    EXPECT_NEAR_ABS(dy_dx1_adhoc, dy_dx1_res, 1e-10);
    EXPECT_NEAR_ABS(dy_dx2_adhoc, dy_dx2_res, 1e-10);

    double d2y_dx1_dx1_adhoc = tape.get_derivative(x1_adhoc, x1_adhoc);
    double d2y_dx1_dx2_adhoc = tape.get_derivative(x1_adhoc, x2_adhoc);
    double d2y_dx2_dx2_adhoc = tape.get_derivative(x2_adhoc, x2_adhoc);

    EXPECT_NEAR_ABS(d2y_dx1_dx1_adhoc, d2y_dx1_dx1_res, 1e-10);
    EXPECT_NEAR_ABS(d2y_dx1_dx2_adhoc, d2y_dx1_dx2_res, 1e-10);
    EXPECT_NEAR_ABS(d2y_dx2_dx2_adhoc, d2y_dx2_dx2_res, 1e-10);
}

void
test_exp2()
{
    adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
    auto& tape = *tapeptr;
    tape.set_method(adhoc::Method::SecondOrderSimple);
    adhoc_t x_adhoc = 2.0;
    tape.register_variable(x_adhoc);
    adhoc_t y_adhoc{ exp(x_adhoc) };
    tape.register_output_variable(y_adhoc);
    tape.set_derivative(y_adhoc, 1.0);
    tape.backpropagate();
    double dy_dx_adhoc = tape.get_derivative(x_adhoc);
    double d2y_dx2_adhoc = tape.get_derivative(x_adhoc, x_adhoc);
    EXPECT_NEAR_ABS(y_adhoc.get_value(), 7.38905609893065, 1e-10);
    EXPECT_NEAR_ABS(dy_dx_adhoc, 7.38905609893065, 1e-10);
    EXPECT_NEAR_ABS(d2y_dx2_adhoc, 7.38905609893065, 1e-10);

    // from sympy import *
    // X1 = Symbol("X1")
    // result = exp(X1)
    // def diff(f, diffs):
    //     fdiff = f
    //     for i in diffs:
    //         fdiff = fdiff.diff(i)
    //     return fdiff
    // X10 = 2.0
    // print(lambdify([X1], result)(X10))
    // print(lambdify([X1], diff(result, (X1,)))(X10))
    // print(lambdify([X1], diff(result, (X1, X1)))(X10))
}

void
test_log2()
{
    adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
    auto& tape = *tapeptr;
    tape.set_method(adhoc::Method::SecondOrderSimple);
    adhoc_t x_adhoc = 2.0;
    tape.register_variable(x_adhoc);
    adhoc_t y_adhoc{ log(x_adhoc) };
    tape.register_output_variable(y_adhoc);
    tape.set_derivative(y_adhoc, 1.0);
    tape.backpropagate();
    double dy_dx_adhoc = tape.get_derivative(x_adhoc);
    double d2y_dx2_adhoc = tape.get_derivative(x_adhoc, x_adhoc);
    EXPECT_NEAR_ABS(y_adhoc.get_value(), 0.6931471805599453, 1e-10);
    EXPECT_NEAR_ABS(dy_dx_adhoc, 0.5, 1e-10);
    EXPECT_NEAR_ABS(d2y_dx2_adhoc, -0.25, 1e-10);

    // from sympy import *
    // X1 = Symbol("X1")
    // result = log(X1)
    // def diff(f, diffs):
    //     fdiff = f
    //     for i in diffs:
    //         fdiff = fdiff.diff(i)
    //     return fdiff
    // X10 = 2.0
    // print(lambdify([X1], result)(X10))
    // print(lambdify([X1], diff(result, (X1,)))(X10))
    // print(lambdify([X1], diff(result, (X1, X1)))(X10))
}

void
test_register()
{
    // write a test where you first register a variable,
    //  do some calculations with it, register another variable,
    //  do some calculations with both variables, and then backpropagate
    //  to check that the derivatives are correct

    double dy_dx1_res = 0;
    double dy_dx2_res = 0;

    adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
    auto& tape = *tapeptr;
    adhoc_t x1_adhoc = 3.0;
    tape.register_variable(x1_adhoc);

    // First calculation with x1 only
    adhoc_t temp1_adhoc = x1_adhoc * x1_adhoc + 2.0 * x1_adhoc; // temp1 = x1^2 + 2*x1

    // Register second variable
    adhoc_t x2_adhoc = 5.0;
    tape.register_variable(x2_adhoc);

    // Calculations with both variables
    adhoc_t y_adhoc = temp1_adhoc * x2_adhoc + x1_adhoc * x2_adhoc * x2_adhoc; // y = (x1^2 + 2*x1)*x2 + x1*x2^2

    tape.register_output_variable(y_adhoc);
    tape.set_derivative(y_adhoc, 1.0);
    tape.backpropagate();
    double dy_dx1_adhoc = tape.get_derivative(x1_adhoc);
    double dy_dx2_adhoc = tape.get_derivative(x2_adhoc);

    // Analytical derivatives:
    // y = (x1^2 + 2*x1)*x2 + x1*x2^2
    // dy/dx1 = (2*x1 + 2)*x2 + x2^2 = (2*3 + 2)*5 + 25 = 8*5 + 25 = 40 + 25 =
    // 65 dy/dx2 = (x1^2 + 2*x1) + 2*x1*x2 = (9 + 6) + 2*3*5 = 15 + 30 = 45

    EXPECT_NEAR_ABS(dy_dx1_adhoc, dy_dx1_res, 1e-10);
    EXPECT_NEAR_ABS(dy_dx2_adhoc, dy_dx2_res, 1e-10);
}

void
test_abs()
{
    auto check1d = [](double input, double expected_derivative) {
        double dy_dx_res = 0;

        adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
        auto& tape = *tapeptr;
        adhoc_t x_adhoc = input;
        tape.register_variable(x_adhoc);
        adhoc_t y_adhoc = abs(x_adhoc);
        tape.register_output_variable(y_adhoc);
        tape.set_derivative(y_adhoc, 1.0);
        tape.backpropagate();
        double dy_dx_adhoc = tape.get_derivative(x_adhoc);

        EXPECT_NEAR_ABS(dy_dx_adhoc, dy_dx_res, 1e-10);
        EXPECT_NEAR_ABS(dy_dx_adhoc, expected_derivative, 1e-10);
    };
    check1d(-4.0, -1.0);
    check1d(0.0, 1.0); // derivative at 0 is undefined, we expect 1.0 here
    check1d(5.0, 1.0);

    auto check2d = [](double input, double expected_derivative, double expected_second_derivative) {
        adhoc::smart_tape_ptr_t<adhoc_mode> tape2ptr;
        auto& tape2 = *tape2ptr;
        tape2.set_method(adhoc::Method::SecondOrderSimple);
        adhoc_t x_adhoc = input;
        tape2.register_variable(x_adhoc);
        auto y_adhoc = abs(cos(x_adhoc));
        auto res_adhoc = y_adhoc * exp(y_adhoc);

        tape2.register_output_variable(res_adhoc);
        tape2.set_derivative(res_adhoc, 1.0);
        tape2.backpropagate();
        double dy_dx_adhoc = tape2.get_derivative(x_adhoc);
        double d2y_dx2_adhoc = tape2.get_derivative(x_adhoc, x_adhoc);

        EXPECT_NEAR_ABS(dy_dx_adhoc, expected_derivative, 1e-10);
        EXPECT_NEAR_ABS(d2y_dx2_adhoc, expected_second_derivative, 1e-10);
    };
    check2d(-4.0, 2.4060147912038476, 0.84395515018777223);
    check2d(0.0, 0.0, -5.4365636569180902);
    check2d(5.0, 1.6346620580366851, 2.3050919848270324);
}

void
test_sqrt()
{
    auto check1d = [](double input, double expected_derivative) {
        adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
        auto& tape = *tapeptr;
        adhoc_t x_adhoc = input;
        tape.register_variable(x_adhoc);
        adhoc_t y_adhoc = sqrt(x_adhoc);
        tape.register_output_variable(y_adhoc);
        tape.set_derivative(y_adhoc, 1.0);
        tape.backpropagate();
        double dy_dx_adhoc = tape.get_derivative(x_adhoc);

        EXPECT_NEAR_ABS(dy_dx_adhoc, expected_derivative, 1e-10);
    };
    check1d(1.0, 0.5);       // d/dx sqrt(x) = 1/(2*sqrt(x)) = 0.5
    check1d(4.0, 0.25);      // d/dx sqrt(x) = 1/(2*sqrt(4)) = 0.25
    check1d(9.0, 1.0 / 6.0); // d/dx sqrt(x) = 1/(2*sqrt(9)) = 1/6

    auto check2d = [](double input, double expected_derivative, double expected_second_derivative) {
        adhoc::smart_tape_ptr_t<adhoc_mode> tape2ptr;
        auto& tape2 = *tape2ptr;
        tape2.set_method(adhoc::Method::SecondOrderSimple);
        adhoc_t x_adhoc = input;
        tape2.register_variable(x_adhoc);

        auto y_adhoc = sqrt(cos(x_adhoc) + 2.0);
        auto res_adhoc = y_adhoc * exp(y_adhoc);

        tape2.register_output_variable(res_adhoc);
        tape2.set_derivative(res_adhoc, 1.0);
        tape2.backpropagate();
        double dy_dx_adhoc = tape2.get_derivative(x_adhoc);
        double d2y_dx2_adhoc = tape2.get_derivative(x_adhoc, x_adhoc);

        EXPECT_NEAR_ABS(dy_dx_adhoc, expected_derivative, 1e-10);
        EXPECT_NEAR_ABS(d2y_dx2_adhoc, expected_second_derivative, 1e-10);
    };
    check2d(-4.0, -2.2480967722953933, 2.3823277033005299);
    check2d(0.0, 0.0, -4.4577761536301708);
    check2d(0.33, -1.4268301618141512, -4.059825705671301);
    check2d(2.5, -1.7112657152648361, 2.5545714737289607);
    check2d(1.0, -3.3705524251026877, -1.4896779009391881);
}

void
test_pow_c()
{
    auto check1d = [](double input, double exponent, double expected_derivative) {
        adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
        auto& tape = *tapeptr;
        adhoc_t x_adhoc = input;
        tape.register_variable(x_adhoc);
        adhoc_t y_adhoc = pow(x_adhoc, exponent);
        tape.register_output_variable(y_adhoc);
        tape.set_derivative(y_adhoc, 1.0);
        tape.backpropagate();
        double dy_dx_adhoc = tape.get_derivative(x_adhoc);

        EXPECT_NEAR_ABS(dy_dx_adhoc, expected_derivative, 1e-10);
    };
    // d/dx x^b = b * x^(b-1)
    check1d(2.0, 3.0, 12.0);   // d/dx x^3 = 3*x^2 = 3*4 = 12
    check1d(3.0, 2.0, 6.0);    // d/dx x^2 = 2*x = 2*3 = 6
    check1d(4.0, 0.5, 0.25);   // d/dx x^0.5 = 0.5*x^(-0.5) = 0.5/2 = 0.25
    check1d(2.0, -1.0, -0.25); // d/dx x^(-1) = -x^(-2) = -1/4 = -0.25

    auto check2d = [](double input, double exponent, double expected_derivative, double expected_second_derivative) {
        adhoc::smart_tape_ptr_t<adhoc_mode> tape2ptr;
        auto& tape2 = *tape2ptr;
        tape2.set_method(adhoc::Method::SecondOrderSimple);
        adhoc_t x_adhoc = input;
        tape2.register_variable(x_adhoc);

        auto y_adhoc = pow(cos(x_adhoc) + 2.0, exponent);
        auto res_adhoc = y_adhoc * exp(y_adhoc);

        tape2.register_output_variable(res_adhoc);
        tape2.set_derivative(res_adhoc, 1.0);
        tape2.backpropagate();
        double dy_dx_adhoc = tape2.get_derivative(x_adhoc);
        double d2y_dx2_adhoc = tape2.get_derivative(x_adhoc, x_adhoc);

        EXPECT_NEAR_ABS(dy_dx_adhoc, expected_derivative, 1e-10);
        EXPECT_NEAR_ABS(d2y_dx2_adhoc, expected_second_derivative, 1e-10);
    };
    check2d(-4.0, 2.0, -35.117790250535428, 147.07958295968521);
    check2d(0.0, 2.0, 0.0, -486185.03565452306);
    check2d(0.33, 2.0, -108648.38063053414, -76375.121410327905);
    check2d(2.5, 3.0, -39.361779410815053, 230.86327275804263);
    check2d(1.0, 0.5, -3.3705524251026877, -1.4896779009391881);
}

void
test_erf()
{
    auto check1d = [](double input, double expected_derivative) {
        adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
        auto& tape = *tapeptr;
        adhoc_t x_adhoc = input;
        tape.register_variable(x_adhoc);
        adhoc_t y_adhoc = erf(x_adhoc);
        tape.register_output_variable(y_adhoc);
        tape.set_derivative(y_adhoc, 1.0);
        tape.backpropagate();
        double dy_dx_adhoc = tape.get_derivative(x_adhoc);

        EXPECT_NEAR_ABS(dy_dx_adhoc, expected_derivative, 1e-10);
    };
    // d/dx erf(x) = 2/sqrt(pi) * exp(-x^2)
    check1d(0.0, 1.1283791670955126);   // 2/sqrt(pi) * exp(0) = 2/sqrt(pi)
    check1d(1.0, 0.41510749742059439);  // 2/sqrt(pi) * exp(-1)
    check1d(-1.0, 0.41510749742059439); // symmetric function derivative

    check1d(-4.0, 1.2698234671866558E-7);
    check1d(0.0, 1.1283791670955126);
    check1d(0.33, 1.011953111857459);
    check1d(2.5, 0.0021782842303527099);
    check1d(1.0, 0.41510749742059472);

    auto check2d = [](double input, double expected_derivative, double expected_second_derivative) {
        adhoc::smart_tape_ptr_t<adhoc_mode> tape2ptr;
        auto& tape2 = *tape2ptr;
        tape2.set_method(adhoc::Method::SecondOrderSimple);
        adhoc_t x_adhoc = input;
        tape2.register_variable(x_adhoc);

        auto y_adhoc = erf(cos(x_adhoc));
        auto res_adhoc = y_adhoc * exp(y_adhoc);

        tape2.register_output_variable(res_adhoc);
        tape2.set_derivative(res_adhoc, 1.0);
        tape2.backpropagate();
        double dy_dx_adhoc = tape2.get_derivative(x_adhoc);
        double d2y_dx2_adhoc = tape2.get_derivative(x_adhoc, x_adhoc);

        EXPECT_NEAR_ABS(dy_dx_adhoc, expected_derivative, 1e-10);
        EXPECT_NEAR_ABS(d2y_dx2_adhoc, expected_second_derivative, 1e-10);
    };
    check2d(-4.0, 1.2698234671866558E-7, 0.0000010158587737493247);
    check2d(0.0, 1.1283791670955126, 0.);
    check2d(0.33, 1.011953111857459, -0.6678890538259229);
    check2d(2.5, 0.0021782842303527099, -0.01089142115176355);
    check2d(1.0, 0.41510749742059472, -0.83021499484118944);
}

void
test_erfc()
{
    auto check1d = [](double input, double expected_derivative) {
        adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
        auto& tape = *tapeptr;
        adhoc_t x_adhoc = input;
        tape.register_variable(x_adhoc);
        adhoc_t y_adhoc = erfc(x_adhoc);
        tape.register_output_variable(y_adhoc);
        tape.set_derivative(y_adhoc, 1.0);
        tape.backpropagate();
        double dy_dx_adhoc = tape.get_derivative(x_adhoc);

        EXPECT_NEAR_ABS(dy_dx_adhoc, expected_derivative, 1e-10);
    };
    // d/dx erf(x) = 2/sqrt(pi) * exp(-x^2)
    check1d(0.0, 1.1283791670955126);   // 2/sqrt(pi) * exp(0) = 2/sqrt(pi)
    check1d(1.0, 0.41510749742059439);  // 2/sqrt(pi) * exp(-1)
    check1d(-1.0, 0.41510749742059439); // symmetric function derivative

    check1d(-4.0, 1.2698234671866558E-7);
    check1d(0.0, 1.1283791670955126);
    check1d(0.33, 1.011953111857459);
    check1d(2.5, 0.0021782842303527099);
    check1d(1.0, 0.41510749742059472);

    auto check2d = [](double input, double expected_derivative, double expected_second_derivative) {
        adhoc::smart_tape_ptr_t<adhoc_mode> tape2ptr;
        auto& tape2 = *tape2ptr;
        tape2.set_method(adhoc::Method::SecondOrderSimple);
        adhoc_t x_adhoc = input;
        tape2.register_variable(x_adhoc);

        auto y_adhoc = erfc(cos(x_adhoc));
        auto res_adhoc = y_adhoc * exp(y_adhoc);

        tape2.register_output_variable(res_adhoc);
        tape2.set_derivative(res_adhoc, 1.0);
        tape2.backpropagate();
        double dy_dx_adhoc = tape2.get_derivative(x_adhoc);
        double d2y_dx2_adhoc = tape2.get_derivative(x_adhoc, x_adhoc);

        EXPECT_NEAR_ABS(dy_dx_adhoc, expected_derivative, 1e-10);
        EXPECT_NEAR_ABS(d2y_dx2_adhoc, expected_second_derivative, 1e-10);
    };
    check2d(-4.0, 1.2698234671866558E-7, 0.0000010158587737493247);
    check2d(0.0, 1.1283791670955126, 0.);
    check2d(0.33, 1.011953111857459, -0.6678890538259229);
    check2d(2.5, 0.0021782842303527099, -0.01089142115176355);
    check2d(1.0, 0.41510749742059472, -0.83021499484118944);
}

void
test_skip()
{

    auto check2d = [](double input, double expected_derivative) {
        adhoc::smart_tape_ptr_t<adhoc_mode> tape2ptr;
        auto& tape2 = *tape2ptr;
        tape2.set_method(adhoc::Method::SecondOrderSimple);
        adhoc_t x_adhoc = input;
        tape2.register_variable(x_adhoc);

        auto y_adhoc = erf(cos(x_adhoc));
        auto res_adhoc = y_adhoc * exp(y_adhoc);
        auto skip_adhoc = y_adhoc * exp(res_adhoc);

        tape2.register_output_variable(res_adhoc);
        tape2.set_derivative(res_adhoc, 1.0);
        tape2.backpropagate();
        double dy_dx_adhoc = tape2.get_derivative(x_adhoc);

        EXPECT_NEAR_ABS(dy_dx_adhoc, expected_derivative, 1e-10);
    };
    check2d(-4.0, 1.2698234671866558E-7);
    check2d(0.0, 1.1283791670955126);
    check2d(0.33, 1.011953111857459);
    check2d(2.5, 0.0021782842303527099);
    check2d(1.0, 0.41510749742059472);
}

void
test_backprop_stages()
{
    auto check2d = [](double input, double expected_derivative) -> void {
        double dy_dx_adhoc = 0;
        {
            using mode = adhoc::opcode<double>;
            using D = mode::type;
            adhoc::smart_tape_ptr_t<mode> tape;

            D x_adhoc = input;
            tape->register_variable(x_adhoc);

            auto y_adhoc = erf(cos(x_adhoc));
            auto res_adhoc = y_adhoc * exp(y_adhoc);
            auto skip_adhoc = y_adhoc * exp(res_adhoc);

            tape->register_output_variable(res_adhoc);
            tape->set_derivative(res_adhoc, 1.0);
            tape->backpropagate_and_reset_to<true>(3);
            tape->backpropagate();

            dy_dx_adhoc = tape->get_derivative(x_adhoc);
        }

        EXPECT_NEAR_ABS(dy_dx_adhoc, expected_derivative, 1e-10);
    };
    check2d(-4.0, 1.2698234671866558E-7);
    check2d(0.0, 1.1283791670955126);
    check2d(0.33, 1.011953111857459);
    check2d(2.5, 0.0021782842303527099);
    check2d(1.0, 0.41510749742059472);
}

void
test_max()
{
    auto check1d = [](double input1, double input2, double expected_derivative1, double expected_derivative2) -> void {
        adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
        auto& tape = *tapeptr;
        adhoc_t x_adhoc1 = input1;
        adhoc_t x_adhoc2 = input2;
        tape.register_variable(x_adhoc1);
        tape.register_variable(x_adhoc2);
        adhoc_t y_adhoc = max(x_adhoc1, x_adhoc2);
        tape.register_output_variable(y_adhoc);
        tape.set_derivative(y_adhoc, 1.0);
        tape.backpropagate();
        double dy_dx_adhoc1 = tape.get_derivative(x_adhoc1);
        double dy_dx_adhoc2 = tape.get_derivative(x_adhoc2);

        EXPECT_NEAR_ABS(dy_dx_adhoc1, expected_derivative1, 1e-10);
        EXPECT_NEAR_ABS(dy_dx_adhoc2, expected_derivative2, 1e-10);
    };
    check1d(-4.0, -1.0, 1.2698234671866558E-7, 0.0);
    check1d(0.0, 1.0, 1.1283791670955126, 0.0); // derivative at 0 is undefined, we expect 1.0 here
    check1d(5.0, 1.0, 0.0, 1.0);
    check1d(3.0, 3.0, 0.5, 0.5);

    auto check2d = [](double input1,
                      double input2,
                      double expected_derivative1,
                      double expected_derivative2,
                      double expected_derivative3,
                      double expected_derivative4,
                      double expected_derivative5) {
        double dy_dx_res1, dy_dx_res2, d2y_dx2_res11, d2y_dx2_res12, d2y_dx2_res22;

        adhoc::smart_tape_ptr_t<adhoc_mode> tape2ptr;
        auto& tape2 = *tape2ptr;
        tape2.set_method(adhoc::Method::SecondOrderSimple);
        adhoc_t x_adhoc1 = input1;
        adhoc_t x_adhoc2 = input2;
        tape2.register_variable(x_adhoc1);
        tape2.register_variable(x_adhoc2);

        auto y_adhoc = max(cos(x_adhoc1), cos(x_adhoc2));
        auto res_adhoc = y_adhoc * exp(y_adhoc);

        tape2.register_output_variable(res_adhoc);
        tape2.set_derivative(res_adhoc, 1.0);
        tape2.backpropagate();
        double dy_dx_adhoc1 = tape2.get_derivative(x_adhoc1);
        double dy_dx_adhoc2 = tape2.get_derivative(x_adhoc2);
        double d2y_dx2_adhoc11 = tape2.get_derivative(x_adhoc1, x_adhoc1);
        double d2y_dx2_adhoc12 = tape2.get_derivative(x_adhoc1, x_adhoc2);
        double d2y_dx2_adhoc22 = tape2.get_derivative(x_adhoc2, x_adhoc2);

        EXPECT_NEAR_ABS(dy_dx_adhoc1, expected_derivative1, 1e-10);
        EXPECT_NEAR_ABS(dy_dx_adhoc2, expected_derivative2, 1e-10);
        EXPECT_NEAR_ABS(d2y_dx2_adhoc11, expected_derivative3, 1e-10);
        EXPECT_NEAR_ABS(d2y_dx2_adhoc12, expected_derivative4, 1e-10);
        EXPECT_NEAR_ABS(d2y_dx2_adhoc22, expected_derivative5, 1e-10);
    };
    check2d(-4.0, 2., 0, 0, 0, 0, 0);
    check2d(0.0, 0.0, 0, 0, 0, 0, 0);
    check2d(5.0, 1.6, 0, 0, 0, 0, 0);
}

void
test_min()
{
    auto check1d = [](double input1, double input2, double expected_derivative1, double expected_derivative2) -> void {
        adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
        auto& tape = *tapeptr;
        adhoc_t x_adhoc1 = input1;
        adhoc_t x_adhoc2 = input2;
        tape.register_variable(x_adhoc1);
        tape.register_variable(x_adhoc2);
        adhoc_t y_adhoc = min(x_adhoc1, x_adhoc2);
        tape.register_output_variable(y_adhoc);
        tape.set_derivative(y_adhoc, 1.0);
        tape.backpropagate();
        double dy_dx_adhoc1 = tape.get_derivative(x_adhoc1);
        double dy_dx_adhoc2 = tape.get_derivative(x_adhoc2);

        EXPECT_NEAR_ABS(dy_dx_adhoc1, expected_derivative1, 1e-10);
        EXPECT_NEAR_ABS(dy_dx_adhoc2, expected_derivative2, 1e-10);
    };
    check1d(-4.0, -1.0, 0, 0);
    check1d(0.0, 1.0, 0, 0); // derivative at 0 is undefined, we expect 1.0 here
    check1d(5.0, 1.0, 0, 0);
    check1d(3.0, 3.0, 0, 0);

    auto check2d = [](double input1,
                      double input2,
                      double expected_derivative1,
                      double expected_derivative2,
                      double expected_derivative3,
                      double expected_derivative4,
                      double expected_derivative5) {
        double dy_dx_res1, dy_dx_res2, d2y_dx2_res11, d2y_dx2_res12, d2y_dx2_res22;

        adhoc::smart_tape_ptr_t<adhoc_mode> tape2ptr;
        auto& tape2 = *tape2ptr;
        tape2.set_method(adhoc::Method::SecondOrderSimple);
        adhoc_t x_adhoc1 = input1;
        adhoc_t x_adhoc2 = input2;
        tape2.register_variable(x_adhoc1);
        tape2.register_variable(x_adhoc2);

        auto y_adhoc = min(cos(x_adhoc1), cos(x_adhoc2));
        auto res_adhoc = y_adhoc * exp(y_adhoc);

        tape2.register_output_variable(res_adhoc);
        tape2.set_derivative(res_adhoc, 1.0);
        tape2.backpropagate();
        double dy_dx_adhoc1 = tape2.get_derivative(x_adhoc1);
        double dy_dx_adhoc2 = tape2.get_derivative(x_adhoc2);
        double d2y_dx2_adhoc11 = tape2.get_derivative(x_adhoc1, x_adhoc1);
        double d2y_dx2_adhoc12 = tape2.get_derivative(x_adhoc1, x_adhoc2);
        double d2y_dx2_adhoc22 = tape2.get_derivative(x_adhoc2, x_adhoc2);

        EXPECT_NEAR_ABS(dy_dx_adhoc1, expected_derivative1, 1e-10);
        EXPECT_NEAR_ABS(dy_dx_adhoc2, expected_derivative2, 1e-10);
        EXPECT_NEAR_ABS(d2y_dx2_adhoc11, expected_derivative3, 1e-10);
        EXPECT_NEAR_ABS(d2y_dx2_adhoc12, expected_derivative4, 1e-10);
        EXPECT_NEAR_ABS(d2y_dx2_adhoc22, expected_derivative5, 1e-10);
    };
    check2d(-4.0, 2., -0.283943, 0.909297, 0.168414, -0.283943, 0.817556);
    check2d(0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    check2d(5.0, 1.6, 0.993297, 0.999576, 0.986548, 0.993297, 0.999576);
}

void
test_expm1()
{
    auto check1d = [](double input, double expected_derivative) -> void {
        adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
        auto& tape = *tapeptr;
        adhoc_t x_adhoc = input;
        tape.register_variable(x_adhoc);
        adhoc_t y_adhoc = expm1(x_adhoc);
        tape.register_output_variable(y_adhoc);
        tape.set_derivative(y_adhoc, 1.0);
        tape.backpropagate();
        double dy_dx_adhoc = tape.get_derivative(x_adhoc);

        EXPECT_NEAR_ABS(dy_dx_adhoc, expected_derivative, 1e-10);
    };
    check1d(-4.0, exp(-4.0));
    check1d(0.0, 1.0);
    check1d(0.33, exp(0.33));
    check1d(2.5, exp(2.5));
    check1d(1.0, exp(1.0));

    // second order
    auto check2d = [](double input, double expected_derivative1, double expected_derivative2) {
        adhoc::smart_tape_ptr_t<adhoc_mode> tape2ptr;
        auto& tape2 = *tape2ptr;
        tape2.set_method(adhoc::Method::SecondOrderSimple);
        adhoc_t x_adhoc = input;
        tape2.register_variable(x_adhoc);

        auto y_adhoc = expm1(x_adhoc);
        auto res_adhoc = y_adhoc * exp(y_adhoc);

        tape2.register_output_variable(res_adhoc);
        tape2.set_derivative(res_adhoc, 1.0);
        tape2.backpropagate();
        double dy_dx_adhoc = tape2.get_derivative(x_adhoc);
        double d2y_dx2_adhoc = tape2.get_derivative(x_adhoc, x_adhoc);

        EXPECT_NEAR_ABS(dy_dx_adhoc, expected_derivative1, 3e-8);
        EXPECT_NEAR_ABS(d2y_dx2_adhoc, expected_derivative2, 3e-8);
    };
    check2d(-4.0, exp(-4.0), exp(-4.0));
    check2d(0.0, 1.0, 1.0);
    check2d(0.33, exp(0.33), exp(0.33));
    check2d(2.5, exp(2.5), exp(2.5));
    check2d(1.0, exp(1.0), exp(1.0));
}

void
test_induced_path()
{
    adhoc::smart_tape_ptr_t<adhoc_mode> tape2ptr;
    auto& tape2 = *tape2ptr;
    tape2.set_method(adhoc::Method::FirstOrderLossyCompressed);

    adhoc::opcode<double>::type x_adhoc = 0.5;
    tape2.register_variable(x_adhoc);
    auto res = exp(erf(erfc(cos(x_adhoc))));
    tape2.register_output_variable(res);
    tape2.set_derivative(res, 1.0);
    tape2.backpropagate();
    double dy_dx_adhoc = tape2.get_derivative(x_adhoc);
    EXPECT_NEAR_ABS(dy_dx_adhoc, 0.34254904163922512, 1e-13);
}

auto
main() -> int
{
    test_compound_division();
    test_induced_path();
    test_expm1();
    test_min();
    test_max();
    test_erfc();
    test_skip();
    test_erf();
    test_pow_c();
    test_abs();
    test_sqrt();

    // order 2 tests
    test_division2();
    test_exp2();
    test_log2();

    // order 1 tests
    test_addition();
    test_multiplication();
    test_division();
    test_exp();
    test_log();
    test_compound();
    test_compound_assignment();
    test_compound_subtraction();
    test_compound_multiplication();
    test_compound_division();
    test_multiple_compound_assignments();
    test_first_cash_instrument();
    test_copy();
    test_register();

    TEST_END;
}
