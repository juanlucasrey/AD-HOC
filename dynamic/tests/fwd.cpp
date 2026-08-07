#include <tape.hpp>
#include <test_simple_include.hpp>

#include "mini_mc.hpp"

void
test_float()
{
    using adhoc_mode = adhoc::opcode<float>;
    using adhoc_t = adhoc_mode::type;

    adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
    auto& tape = *tapeptr;
    tape.configure(adhoc::Method::Bwd, 2, 1);
    adhoc_t x1_adhoc = 3.0;
    adhoc_t x2_adhoc = 5.0;

    unsigned int temp = 3;
    adhoc_t x3_adhoc = x2_adhoc * temp;
    auto const valx1 = x1_adhoc.get_value();
    static_assert(std::is_same_v<decltype(valx1), float const>);

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

using adhoc_mode = adhoc::opcode<double>;
using adhoc_t = adhoc_mode::type;

void
test_2to2()
{
    double dy1_dx1_adhoc_res = 0.0;
    double dy1_dx2_adhoc_res = 0.0;
    double dy2_dx1_adhoc_res = 0.0;
    double dy2_dx2_adhoc_res = 0.0;
    {
        // bwd
        adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
        auto& tape = *tapeptr;
        tape.configure(adhoc::Method::Bwd, 2, 2);
        adhoc_t x1_adhoc = 3.0;
        adhoc_t x2_adhoc = 5.0;
        tape.register_variable(x1_adhoc);
        tape.register_variable(x2_adhoc);
        adhoc_t v1 = (x1_adhoc + x2_adhoc) * exp(x1_adhoc * x2_adhoc);
        adhoc_t v2 = erfc(v1) + log(x2_adhoc);
        adhoc_t y1_adhoc = v1 + v2;
        adhoc_t y2_adhoc = v1 * v2;
        tape.register_output_variable(y1_adhoc);
        tape.register_output_variable(y2_adhoc);
        tape.set_derivative(y1_adhoc, 1.0, 0);
        tape.set_derivative(y2_adhoc, 1.0, 1);
        tape.backpropagate();
        dy1_dx1_adhoc_res = tape.get_derivative(x1_adhoc, 0);
        dy1_dx2_adhoc_res = tape.get_derivative(x2_adhoc, 0);
        dy2_dx1_adhoc_res = tape.get_derivative(x1_adhoc, 1);
        dy2_dx2_adhoc_res = tape.get_derivative(x2_adhoc, 1);
    }
    {
        adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
        auto& tape = *tapeptr;
        tape.configure(adhoc::Method::Fwd, 2, 2);
        adhoc_t x1_adhoc = 3.0;
        adhoc_t x2_adhoc = 5.0;
        tape.register_variable(x1_adhoc);
        tape.register_variable(x2_adhoc);
        adhoc_t v1 = (x1_adhoc + x2_adhoc) * exp(x1_adhoc * x2_adhoc);
        adhoc_t v2 = erfc(v1) + log(x2_adhoc);
        adhoc_t y1_adhoc = v1 + v2;
        adhoc_t y2_adhoc = v1 * v2;
        tape.register_output_variable(y1_adhoc);
        tape.register_output_variable(y2_adhoc);
        tape.set_derivative(x1_adhoc, 1.0, 0);
        tape.set_derivative(x2_adhoc, 1.0, 1);
        tape.backpropagate();
        double dy1_dx1_adhoc = tape.get_derivative(0, 0);
        double dy1_dx2_adhoc = tape.get_derivative(1, 0);
        double dy2_dx1_adhoc = tape.get_derivative(0, 1);
        double dy2_dx2_adhoc = tape.get_derivative(1, 1);
        EXPECT_NEAR_ABS(dy1_dx1_adhoc_res, dy1_dx1_adhoc, 1e-10);
        EXPECT_NEAR_ABS(dy1_dx2_adhoc_res, dy1_dx2_adhoc, 1e-10);
        EXPECT_NEAR_ABS(dy2_dx1_adhoc_res, dy2_dx1_adhoc, 1e-10);
        EXPECT_NEAR_ABS(dy2_dx2_adhoc_res, dy2_dx2_adhoc, 1e-10);
    }
}

void
test_addition()
{
    adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
    auto& tape = *tapeptr;
    tape.configure(adhoc::Method::Fwd, 2, 1);
    adhoc_t x1_adhoc = 3.0;
    adhoc_t x2_adhoc = 5.0;
    tape.register_variable(x1_adhoc);
    tape.register_variable(x2_adhoc);
    adhoc_t y_adhoc = x1_adhoc + x2_adhoc;
    tape.register_output_variable(y_adhoc);

    tape.set_derivative(x1_adhoc, 1.0, 0);
    tape.set_derivative(x2_adhoc, 1.0, 1);
    tape.backpropagate();

    double dy_dx1_adhoc = tape.get_derivative(y_adhoc, 0);
    double dy_dx2_adhoc = tape.get_derivative(y_adhoc, 1);

    EXPECT_NEAR_ABS(dy_dx1_adhoc, 1., 1e-10);
    EXPECT_NEAR_ABS(dy_dx2_adhoc, 1., 1e-10);
}

void
test_multiplication()
{
    adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
    auto& tape = *tapeptr;
    tape.configure(adhoc::Method::Fwd, 2, 1);
    adhoc_t x1_adhoc = 3.0, x2_adhoc = 5.0;
    tape.register_variable(x1_adhoc);
    tape.register_variable(x2_adhoc);
    adhoc_t y_adhoc = x1_adhoc * x2_adhoc;
    tape.register_output_variable(y_adhoc);

    tape.set_derivative(x1_adhoc, 1.0, 0);
    tape.set_derivative(x2_adhoc, 1.0, 1);
    tape.backpropagate();

    double dy_dx1_adhoc = tape.get_derivative(y_adhoc, 0);
    double dy_dx2_adhoc = tape.get_derivative(y_adhoc, 1);

    EXPECT_NEAR_ABS(dy_dx1_adhoc, 5.0, 1e-10);
    EXPECT_NEAR_ABS(dy_dx2_adhoc, 3.0, 1e-10);
}

void
test_division()
{
    adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
    auto& tape = *tapeptr;
    tape.configure(adhoc::Method::Fwd, 2, 1);
    adhoc_t x1_adhoc = 10.0, x2_adhoc = 2.0;
    tape.register_variable(x1_adhoc);
    tape.register_variable(x2_adhoc);
    adhoc_t y_adhoc = x1_adhoc / x2_adhoc;
    tape.register_output_variable(y_adhoc);

    tape.set_derivative(x1_adhoc, 1.0, 0);
    tape.set_derivative(x2_adhoc, 1.0, 1);
    tape.backpropagate();

    double dy_dx1_adhoc = tape.get_derivative(y_adhoc, 0);
    double dy_dx2_adhoc = tape.get_derivative(y_adhoc, 1);

    EXPECT_NEAR_ABS(dy_dx1_adhoc, 0.5, 1e-10);
    EXPECT_NEAR_ABS(dy_dx2_adhoc, -2.5, 1e-10);
}

void
test_first_cash_instrument()
{
    double coverage = 7.0 / 360.0;
    adhoc::smart_tape_ptr_t<adhoc_mode> tapeptr;
    auto& tape = *tapeptr;
    tape.configure(adhoc::Method::Fwd, 3, 1);
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

    tape.set_derivative(df0_adhoc, 1.0, 0);
    tape.set_derivative(df1_adhoc, 1.0, 1);
    tape.set_derivative(rate_adhoc, 1.0, 2);
    tape.backpropagate();

    double dpv_ddf0_adhoc = tape.get_derivative(pv_adhoc, 0);
    double dpv_ddf1_adhoc = tape.get_derivative(pv_adhoc, 1);
    double dpv_drate_adhoc = tape.get_derivative(pv_adhoc, 2);

    EXPECT_NEAR_ABS(dpv_ddf0_adhoc, -1., 1e-10);
    EXPECT_NEAR_ABS(dpv_ddf1_adhoc, 1.0008416615459654, 1e-10);
    EXPECT_NEAR_ABS(dpv_drate_adhoc, 0.019418754810350696, 1e-10);
}

void
test_checkpoint_fd_vs_fwd()
{
    constexpr std::size_t num_paths = 1000;
    double x1_val = 1.5, x2_val = 2.0, x3_val = 0.5;

    double res_fd = 0.;
    double dx1_fd = 0.;
    double dx2_fd = 0.;
    double dx3_fd = 0.;

    {
        constexpr double bump = 1e-6;
        double x1 = x1_val;
        double x2 = x2_val;
        double x3 = x3_val;

        res_fd = compute_result_branch(x1, x2, x3, num_paths);

        // Central finite differences using the same templated function
        dx1_fd =
          (compute_result_branch(x1 + bump, x2, x3, num_paths) - compute_result_branch(x1 - bump, x2, x3, num_paths)) /
          (2.0 * bump);
        dx2_fd =
          (compute_result_branch(x1, x2 + bump, x3, num_paths) - compute_result_branch(x1, x2 - bump, x3, num_paths)) /
          (2.0 * bump);
        dx3_fd =
          (compute_result_branch(x1, x2, x3 + bump, num_paths) - compute_result_branch(x1, x2, x3 - bump, num_paths)) /
          (2.0 * bump);
    }

    double res_adhoc = 0.;
    double dx1_adhoc = 0.;
    double dx2_adhoc = 0.;
    double dx3_adhoc = 0.;

    {
        using adhoc_t = adhoc_t;
        // Create tape
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tapeptr;
        auto& tape = *tapeptr;
        tape.configure(adhoc::Method::Fwd, 3, 1);

        // Initial input variables
        adhoc_t x1, x2, x3;
        x1 = x1_val;
        x2 = x2_val;
        x3 = x3_val;

        // Register inputs
        tape.register_variable(x1);
        tape.register_variable(x2);
        tape.register_variable(x3);

        auto res_adhoc2 = compute_result_branch(x1, x2, x3, num_paths);
        res_adhoc = res_adhoc2.get_value();

        // tape.backpropagate();

        dx1_adhoc = tape.get_derivative(0, 0);
        dx2_adhoc = tape.get_derivative(1, 0);
        dx3_adhoc = tape.get_derivative(2, 0);
    }

    // check that the values are close in between methods
    EXPECT_NEAR_ABS(res_adhoc, res_fd, 1e-8);
    EXPECT_NEAR_ABS(dx1_adhoc, dx1_fd, 1e-8);
    EXPECT_NEAR_ABS(dx2_adhoc, dx2_fd, 1e-8);
    EXPECT_NEAR_ABS(dx3_adhoc, dx3_fd, 1e-8);

    EXPECT_NEAR_ABS(3.421138662549827, res_adhoc, 1e-12);
    EXPECT_NEAR_ABS(1.5027371453017027, dx1_adhoc, 1e-12);
    EXPECT_NEAR_ABS(0.80747795083939378, dx2_adhoc, 1e-12);
    EXPECT_NEAR_ABS(0.11992672558600066, dx3_adhoc, 1e-12);
}

auto
main() -> int
{
    test_2to2();
    test_addition();
    test_multiplication();
    test_division();
    test_first_cash_instrument();
    test_checkpoint_fd_vs_fwd();

    test_float();
    TEST_END;
}
