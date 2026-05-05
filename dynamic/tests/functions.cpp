#include "../../test_simple/test_simple_include.hpp"

#include <adhoc.hpp>
#include <cmath>

using adhoc_mode = adhoc::opcode<double>;
using adhoc_t = adhoc_mode::type;

// ============================================================================
// BASIC ARITHMETIC OPERATIONS
// ============================================================================

// --- Addition ---
// y = x1 + x2
void
test_addition()
{
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

    EXPECT_NEAR_ABS(dy_dx1_adhoc, 1., 1e-10);
    EXPECT_NEAR_ABS(dy_dx2_adhoc, 1., 1e-10);
}

// --- Multiplication ---
// y = x1 * x2
void
test_multiplication()
{
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

    EXPECT_NEAR_ABS(dy_dx1_adhoc, 5.0, 1e-10);
    EXPECT_NEAR_ABS(dy_dx2_adhoc, 3.0, 1e-10);
}

// --- Division ---
// y = x1 / x2
void
test_division_first_order()
{
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

    EXPECT_NEAR_ABS(dy_dx1_adhoc, 0.5, 1e-10);
    EXPECT_NEAR_ABS(dy_dx2_adhoc, -2.5, 1e-10);
}

void
test_division_second_order()
{
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
    double d2y_dx1_dx1_adhoc = tape.get_derivative(x1_adhoc, x1_adhoc);
    double d2y_dx1_dx2_adhoc = tape.get_derivative(x1_adhoc, x2_adhoc);
    double d2y_dx2_dx2_adhoc = tape.get_derivative(x2_adhoc, x2_adhoc);

    EXPECT_NEAR_ABS(dy_dx1_adhoc, 0.5, 1e-10);
    EXPECT_NEAR_ABS(dy_dx2_adhoc, -2.5, 1e-10);
    EXPECT_NEAR_ABS(d2y_dx1_dx1_adhoc, 0., 1e-10);
    EXPECT_NEAR_ABS(d2y_dx1_dx2_adhoc, -0.25, 1e-10);
    EXPECT_NEAR_ABS(d2y_dx2_dx2_adhoc, 2.5, 1e-10);
}

void
test_addition_itself()
{
    adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
    auto& tape = *tapeptr;
    adhoc_t x1_adhoc = 3.0;
    tape.register_variable(x1_adhoc);
    adhoc_t y_adhoc = x1_adhoc + x1_adhoc;
    tape.register_output_variable(y_adhoc);
    tape.set_derivative(y_adhoc, 1.0);
    tape.backpropagate();
    double dy_dx1_adhoc = tape.get_derivative(x1_adhoc);

    EXPECT_NEAR_ABS(dy_dx1_adhoc, 2., 1e-10);
}

void
test_substraction_itself()
{
    adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
    auto& tape = *tapeptr;
    adhoc_t x1_adhoc = 3.0;
    tape.register_variable(x1_adhoc);
    adhoc_t y_adhoc = x1_adhoc - x1_adhoc;
    tape.register_output_variable(y_adhoc);
    tape.set_derivative(y_adhoc, 1.0);
    tape.backpropagate();
    double dy_dx1_adhoc = tape.get_derivative(x1_adhoc);

    EXPECT_NEAR_ABS(dy_dx1_adhoc, 0., 1e-10);
}

void
test_multiplication_itself()
{
    adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
    auto& tape = *tapeptr;
    adhoc_t x1_adhoc = 3.0;
    tape.register_variable(x1_adhoc);
    adhoc_t y_adhoc = x1_adhoc * x1_adhoc;
    tape.register_output_variable(y_adhoc);
    tape.set_derivative(y_adhoc, 1.0);
    tape.backpropagate();
    double dy_dx1_adhoc = tape.get_derivative(x1_adhoc);

    EXPECT_NEAR_ABS(dy_dx1_adhoc, 6., 1e-10);
}

// ============================================================================
// TRANSCENDENTAL FUNCTIONS
// ============================================================================

// --- Exponential ---
// y = exp(x)
void
test_exp_first_order()
{
    adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
    auto& tape = *tapeptr;
    adhoc_t x_adhoc = 2.0;
    tape.register_variable(x_adhoc);
    adhoc_t y_adhoc = exp(x_adhoc);
    tape.register_output_variable(y_adhoc);
    tape.set_derivative(y_adhoc, 1.0);
    tape.backpropagate();
    double dy_dx_adhoc = tape.get_derivative(x_adhoc);

    EXPECT_NEAR_ABS(dy_dx_adhoc, std::exp(2.0), 1e-10);
}

void
test_exp_second_order()
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
}

// --- Expm1 ---
// y = expm1(x) = exp(x) - 1
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
    check2d(-4.0, 0.00012569096013905282, 0.00025368403051559017);
    check2d(0.0, 1.0, 3.0);
    check2d(0.33, 2.8604190310816122, 9.6995897684489946);
    check2d(2.5, 10665171.183325227, 151258725.89737877);
    check2d(1.0, 41.193555674716123, 194.36280518962906);
}

// --- Logarithm ---
// y = log(x)
void
test_log_first_order()
{
    adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
    auto& tape = *tapeptr;
    adhoc_t x_adhoc = 2.0;
    tape.register_variable(x_adhoc);
    adhoc_t y_adhoc = log(x_adhoc);
    tape.register_output_variable(y_adhoc);
    tape.set_derivative(y_adhoc, 1.0);
    tape.backpropagate();
    double dy_dx_adhoc = tape.get_derivative(x_adhoc);

    EXPECT_NEAR_ABS(dy_dx_adhoc, 0.5, 1e-10);
}

void
test_log_second_order()
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
}

// --- Power ---
// y = pow(x, c) where c is a constant exponent
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
test_pow_c_0()
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

        if (expected_derivative == std::numeric_limits<double>::infinity()) {
            EXPECT_TRUE(dy_dx_adhoc == std::numeric_limits<double>::infinity());
        }
        else {
            EXPECT_NEAR_ABS(dy_dx_adhoc, expected_derivative, 1e-10);
        }
    };

    check1d(0.0, 3.0, 0.0);
    check1d(0.0, 1.0001, 0.0);
    check1d(0.0, 1., 1.0);
    check1d(0.0, 0.99, std::numeric_limits<double>::infinity());

    auto check2d = [](double input, double exponent, double expected_derivative, double expected_second_derivative) {
        adhoc::smart_tape_ptr_t<adhoc_mode> tape2ptr;
        auto& tape2 = *tape2ptr;
        tape2.set_method(adhoc::Method::SecondOrderSimple);
        adhoc_t x_adhoc = input;
        tape2.register_variable(x_adhoc);

        auto res_adhoc = pow(x_adhoc, exponent);

        tape2.register_output_variable(res_adhoc);
        tape2.set_derivative(res_adhoc, 1.0);
        tape2.backpropagate();
        double dy_dx_adhoc = tape2.get_derivative(x_adhoc);
        double d2y_dx2_adhoc = tape2.get_derivative(x_adhoc, x_adhoc);

        if (expected_derivative == std::numeric_limits<double>::infinity()) {
            EXPECT_TRUE(dy_dx_adhoc == std::numeric_limits<double>::infinity());
        }
        else {
            EXPECT_NEAR_ABS(dy_dx_adhoc, expected_derivative, 1e-10);
        }

        if (expected_second_derivative == std::numeric_limits<double>::infinity()) {
            EXPECT_TRUE(d2y_dx2_adhoc == std::numeric_limits<double>::infinity());
        }
        else {
            EXPECT_NEAR_ABS(d2y_dx2_adhoc, expected_second_derivative, 1e-10);
        }
    };
    check2d(0.0, 1.0, 1., 0.);
    check2d(0.0, 1.001, 0., std::numeric_limits<double>::infinity());
    check2d(0.0, 1.999, 0., std::numeric_limits<double>::infinity());
    check2d(0.0, 2.0, 0., 2.);
    check2d(0.0, 2.00001, 0., 0.);
}

// --- Square Root ---
// y = sqrt(x)
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

// ============================================================================
// ERROR FUNCTIONS
// ============================================================================

// --- Error Function ---
// y = erf(x)
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
    check1d(0.33, 1.011953111857459);
    check1d(2.5, 0.0021782842303527099);

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
    check2d(-4.0, -0.10386321203665744, 0.41316406932329819);
    check2d(0.0, 0., -1.776624735613332);
    check2d(0.33, -0.61650244812777133, -2.0351208304886157);
    check2d(2.5, -0.043498384426344829, 0.17550835501850212);
    check2d(1.0, -1.9213850008339559, -0.74227271251854887);
}

// --- Complementary Error Function ---
// y = erfc(x) = 1 - erf(x)
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
    // d/dx erfc(x) = -2/sqrt(pi) * exp(-x^2)
    check1d(0.0, -1.1283791670955126);   // -2/sqrt(pi) * exp(0) = -2/sqrt(pi)
    check1d(1.0, -0.41510749742059472);  // -2/sqrt(pi) * exp(-1)
    check1d(-1.0, -0.41510749742059472); // symmetric function derivative
    check1d(-4.0, -1.2698234671866558E-7);
    check1d(0.33, -1.011953111857459);
    check1d(2.5, -0.0021782842303527099);

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
    check2d(-4.0, 7.6305592867620433, -8.28211505327228);
    check2d(0.0, 0., 0.56223827919315505);
    check2d(0.33, 0.21142902393984789, 0.8052344954627304);
    check2d(2.5, 5.5696316027782018, -10.09525571155049);
    check2d(1.0, 1.5984475702829104, 4.3978013617139542);
}

// ============================================================================
// OTHER MATH FUNCTIONS
// ============================================================================

// --- Absolute Value ---
// y = abs(x)
void
test_abs()
{
    auto check1d = [](double input, double expected_derivative) {
        adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
        auto& tape = *tapeptr;
        adhoc_t x_adhoc = input;
        tape.register_variable(x_adhoc);
        adhoc_t y_adhoc = abs(x_adhoc);
        tape.register_output_variable(y_adhoc);
        tape.set_derivative(y_adhoc, 1.0);
        tape.backpropagate();
        double dy_dx_adhoc = tape.get_derivative(x_adhoc);

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

// --- Maximum ---
// y = max(x1, x2)
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
    check1d(-4.0, -1.0, 0., 1.);
    check1d(0.0, 1.0, 0., 1.);
    check1d(5.0, 1.0, 1., 0.);
    check1d(3.0, 3.0, 1., 0.);

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
    check2d(-4.0, 2., 0., -0.35017031350019029, 0., 0., 1.0240250054925624);
    check2d(0.0, 0.0, 0., 0., -5.4365636569180902, 0., 0.);
    check2d(5.0, 1.6, 1.6346620580366851, 0., 2.3050919848270324, 0., 0.);
}

// --- Minimum ---
// y = min(x1, x2)
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
    check1d(-4.0, -1.0, 1., 0.);
    check1d(0.0, 1.0, 1., 0.);
    check1d(5.0, 1.0, 0., 1.);
    check1d(3.0, 3.0, 1., 0.);

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
    check2d(-4.0, 2., -0.13634271204947493, 0., 0.51885675996032143, 0., 0.);
    check2d(0.0, 0.0, 0., 0., -5.4365636569180902, 0., 0.);
    check2d(5.0, 1.6, 0., -0.94246139260191286, 0., 0., 1.9399852799398956);
}

// ============================================================================
// COMPOUND OPERATIONS & ASSIGNMENTS
// ============================================================================

// --- Compound Expression ---
// y = (x1 + x2) * (x3 - x4)
void
test_compound_expression()
{
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

    EXPECT_NEAR_ABS(dy_dx1_adhoc, 3., 1e-10);
    EXPECT_NEAR_ABS(dy_dx2_adhoc, 3., 1e-10);
    EXPECT_NEAR_ABS(dy_dx3_adhoc, 5., 1e-10);
    EXPECT_NEAR_ABS(dy_dx4_adhoc, -5., 1e-10);
}

// --- Compound Assignment += ---
// x1 += x2, y = x1 * x3
void
test_compound_assignment_add()
{
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
    EXPECT_NEAR_ABS(dy_dx1_adhoc, 2., 1e-10);
    EXPECT_NEAR_ABS(dy_dx2_adhoc, 2., 1e-10);
    EXPECT_NEAR_ABS(dy_dx3_adhoc, 8., 1e-10);
}

// --- Compound Assignment -= ---
// x1 -= x2, y = x1 * x3
void
test_compound_assignment_sub()
{
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

    EXPECT_NEAR_ABS(dy_dx1_adhoc, 2., 1e-10);
    EXPECT_NEAR_ABS(dy_dx2_adhoc, -2., 1e-10);
    EXPECT_NEAR_ABS(dy_dx3_adhoc, 7., 1e-10);
}

// --- Compound Assignment *= ---
// x1 *= x2, y = x1 + x3
void
test_compound_assignment_mul()
{
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

    EXPECT_NEAR_ABS(dy_dx1_adhoc, 3., 1e-10); // not 1!
    EXPECT_NEAR_ABS(dy_dx2_adhoc, 4., 1e-10);
    EXPECT_NEAR_ABS(dy_dx3_adhoc, 1., 1e-10);
}

// --- Compound Assignment /= ---
// x1 /= x2, y = x1 * x3
void
test_compound_assignment_div()
{
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

    EXPECT_NEAR_ABS(dy_dx1_adhoc, 1.25, 1e-10);
    EXPECT_NEAR_ABS(dy_dx2_adhoc, -3.75, 1e-10);
    EXPECT_NEAR_ABS(dy_dx3_adhoc, 3.0, 1e-10);
}

// --- Multiple Compound Assignments ---
// x1 += x2, x1 *= x3, y = x1
void
test_multiple_compound_assignments()
{
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

    EXPECT_NEAR_ABS(dy_dx1_adhoc, 4., 1e-10);
    EXPECT_NEAR_ABS(dy_dx2_adhoc, 4., 1e-10);
    EXPECT_NEAR_ABS(dy_dx3_adhoc, 5., 1e-10);
}

// ============================================================================
// TAPE MECHANICS
// ============================================================================

// --- Copy Assignment ---
// y = x1, z = y + x2
void
test_copy_assignment()
{
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

    EXPECT_NEAR_ABS(dz_dx1_adhoc, 1., 1e-10);
    EXPECT_NEAR_ABS(dz_dx2_adhoc, 1., 1e-10);
}

// --- Delayed Registration ---
// Register variables at different times and verify correct derivatives
void
test_delayed_registration()
{
    // Test where you first register a variable,
    // do some calculations with it, register another variable,
    // do some calculations with both variables, and then backpropagate
    // to check that the derivatives are correct

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
    // dy/dx1 = (2*x1 + 2)*x2 + x2^2 = (2*3 + 2)*5 + 25 = 8*5 + 25 = 40 + 25 = 65
    // dy/dx2 = (x1^2 + 2*x1) + 2*x1*x2 = (9 + 6) + 2*3*5 = 15 + 30 = 45

    EXPECT_NEAR_ABS(dy_dx1_adhoc, 65., 1e-10);
    EXPECT_NEAR_ABS(dy_dx2_adhoc, 45., 1e-10);
}

// --- Skip Unregistered Outputs ---
// Verify that unregistered outputs are correctly skipped during backprop
void
test_skip_unregistered_outputs()
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
    check2d(-4.0, -0.10386321203665744);
    check2d(0.0, 0.);
    check2d(0.33, -0.61650244812777133);
    check2d(2.5, -0.043498384426344829);
    check2d(1.0, -1.9213850008339559);
}

// --- Backprop Stages ---
// Test backpropagation with intermediate reset points
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

            // Save tape position after registering x_adhoc
            auto pos = tape->get_position();

            auto y_adhoc = erf(cos(x_adhoc));
            auto res_adhoc = y_adhoc * exp(y_adhoc);
            auto skip_adhoc = y_adhoc * exp(res_adhoc);

            tape->register_output_variable(res_adhoc);
            tape->set_derivative(res_adhoc, 1.0);
            tape->backpropagate_and_reset_to<true>(pos);
            tape->backpropagate();

            dy_dx_adhoc = tape->get_derivative(x_adhoc);
        }

        EXPECT_NEAR_ABS(dy_dx_adhoc, expected_derivative, 1e-10);
    };
    check2d(-4.0, -0.10386321203665747);
    check2d(0.0, 0.);
    check2d(0.33, -0.61650244812777133);
    check2d(2.5, -0.043498384426344829);
    check2d(1.0, -1.9213850008339559);
}

// --- Induced Path ---
// Test FirstOrderLossyCompressed method
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

// ============================================================================
// INTEGRATION / REAL-WORLD SCENARIOS
// ============================================================================

// --- First Cash Instrument ---
// Simplified version of a swap derivative calculation
// Tests the derivative calculation for the first swap (1 week tenor)
// The swap has 1 floating leg and 1 fixed leg, both with a single cashflow
void
test_first_cash_instrument()
{
    // Simplified version: y = x1 * rate * coverage * df1 - df1 + rate *
    // coverage * df1 + df1 where df1 is the discount factor for the 1-week
    // period This simplifies to: y = 2 * x1 * rate * coverage * df1 So dy/dx1 =
    // 2 * rate * coverage * df1 and dy/drate = 2 * x1 * coverage * df1

    double coverage = 7.0 / 360.0;
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

    EXPECT_NEAR_ABS(dpv_ddf0_adhoc, -1., 1e-10);
    EXPECT_NEAR_ABS(dpv_ddf1_adhoc, 1.0008416615459654, 1e-10);
    EXPECT_NEAR_ABS(dpv_drate_adhoc, 0.019418754810350696, 1e-10);
}

auto
main() -> int
{

    // ---- Basic Arithmetic ----
    test_addition();
    test_multiplication();
    test_division_first_order();
    test_division_second_order();
    test_addition_itself();
    test_substraction_itself();
    test_multiplication_itself();

    // ---- Transcendental Functions ----
    test_exp_first_order();
    test_exp_second_order();
    test_expm1();
    test_log_first_order();
    test_log_second_order();
    test_pow_c();
    test_pow_c_0();
    test_sqrt();

    // ---- Error Functions ----
    test_erf();
    test_erfc();

    // ---- Other Math Functions ----
    test_abs();
    test_max();
    test_min();

    // ---- Compound Operations ----
    test_compound_expression();
    test_compound_assignment_add();
    test_compound_assignment_sub();
    test_compound_assignment_mul();
    test_compound_assignment_div();
    test_multiple_compound_assignments();

    // ---- Tape Mechanics ----
    test_copy_assignment();
    test_delayed_registration();
    test_skip_unregistered_outputs();
    test_backprop_stages();
    // test_induced_path();

    // ---- Integration Tests ----
    test_first_cash_instrument();

    TEST_END;
}
