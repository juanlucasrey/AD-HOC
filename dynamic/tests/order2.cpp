#include <tape.hpp>
#include <test_simple_include.hpp>

using adhoc_mode = adhoc::opcode<double>;
using adhoc_t = adhoc_mode::type;

void
order2_simple()
{
    std::vector<adhoc::Method> const methods = { adhoc::Method::SecondOrderSimple, adhoc::Method::SecondOrderLossy };
    auto const safe_method = methods[0];

    double x1_val = 1.5, x2_val = 2.0, x3_val = 0.5;

    double res = 0.;
    double d1 = 0.;
    double d2 = 0.;
    double d3 = 0.;

    double d11 = 0.;
    double d12 = 0.;
    double d13 = 0.;

    double d22 = 0.;
    double d23 = 0.;

    double d33 = 0.;

    for (auto m : methods) {
        using adhoc_t = adhoc_t;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tapeptr;
        auto& tape = *tapeptr;

        tape.set_method(m);

        // Initial input variables
        adhoc_t x1, x2, x3;
        x1 = x1_val;
        x2 = x2_val;
        x3 = x3_val;

        // Register inputs
        tape.register_variable(x1);
        tape.register_variable(x2);
        tape.register_variable(x3);

        double z1 = 0.8;
        double z2 = 0.2;

        auto y_init = x1 * x2;
        auto y_path = y_init * z1 + x1 * z2 + exp(x2 * z1 * 0.1);
        auto y = exp(y_path + x3 * z1 * z2);

        tape.register_output_variable(y);
        tape.set_derivative(y, 1.0);
        tape.backpropagate();

        if (m == safe_method) {
            res = y.get_value();
            d1 = tape.get_derivative(x1);
            d2 = tape.get_derivative(x2);
            d3 = tape.get_derivative(x3);
            d11 = tape.get_derivative(x1, x1);
            d12 = tape.get_derivative(x1, x2);
            d13 = tape.get_derivative(x1, x3);
            d22 = tape.get_derivative(x2, x2);
            d23 = tape.get_derivative(x2, x3);
            d33 = tape.get_derivative(x3, x3);
        }
        else {
            double _res = y.get_value();
            double _d1 = tape.get_derivative(x1);
            double _d2 = tape.get_derivative(x2);
            double _d3 = tape.get_derivative(x3);
            double _d11 = tape.get_derivative(x1, x1);
            double _d12 = tape.get_derivative(x1, x2);
            double _d13 = tape.get_derivative(x1, x3);
            double _d22 = tape.get_derivative(x2, x2);
            double _d23 = tape.get_derivative(x2, x3);
            double _d33 = tape.get_derivative(x3, x3);
            EXPECT_NEAR_ABS(res, _res, 1e-13);
            EXPECT_NEAR_ABS(d1, _d1, 1e-13);
            EXPECT_NEAR_ABS(d2, _d2, 1e-13);
            EXPECT_NEAR_ABS(d3, _d3, 1e-13);
            EXPECT_NEAR_ABS(d11, _d11, 1e-13);
            EXPECT_NEAR_ABS(d12, _d12, 1e-13);
            EXPECT_NEAR_ABS(d13, _d13, 1e-13);
            EXPECT_NEAR_ABS(d22, _d22, 1e-13);
            EXPECT_NEAR_ABS(d23, _d23, 1e-13);
            EXPECT_NEAR_ABS(d33, _d33, 1e-13);
        }
    }
}

auto
main() -> int
{
    order2_simple();
    TEST_END;
}
