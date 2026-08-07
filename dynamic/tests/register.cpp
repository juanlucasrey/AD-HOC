#include <tape.hpp>
#include <test_simple_include.hpp>

namespace {

void
register_input()
{
    using mode_t = adhoc::opcode<float>;
    using type = mode_t::type;

    {
        adhoc::smart_tape_ptr_t<mode_t> tapeptr;
        auto& tape = *tapeptr;
        tape.configure(adhoc::Method::Bwd, 2, 1);
        type X = 1.0;
        tape.register_variable(X);
        auto X1 = X;
        X *= 2.0;
        tape.register_variable(X);
        auto X2 = X;
        X *= 2.0;
        tape.register_output_variable(X);
        tape.set_derivative(X, 1.0);

        tape.backpropagate();

        double val = tape.get_derivative(X1);
        EXPECT_EQUAL(4., tape.get_derivative(X1)); // should be 0??
        // EXPECT_EQUAL(0., tape.get_derivative(X1)); // 0!!
        EXPECT_EQUAL(2., tape.get_derivative(X2));
        EXPECT_EQUAL(4., X.get_value());
    }

    {
        adhoc::smart_tape_ptr_t<mode_t> tape;
        type X = 1.0;
        tape->register_variable(X);
        auto X1 = X;
        X *= 2.0;
        // tape->register_variable(X);
        auto X2 = X;
        X *= 2.0;
        tape->register_output_variable(X);
        tape->set_derivative(X, 1.0);

        tape->backpropagate();

        EXPECT_EQUAL(4., tape->get_derivative(X1)); // 4!!
        EXPECT_EQUAL(2., tape->get_derivative(X2));
        EXPECT_EQUAL(4., X.get_value());
    }
}

void
register_output()
{
    using mode_t = adhoc::opcode<float>;
    using type = mode_t::type;

    {
        adhoc::smart_tape_ptr_t<mode_t> tape;

        type Constant = std::nan("");
        type X = 1.0;
        // tape->register_variable(X);
        type X2 = X * 2.;
        type Y = X2;
        tape->register_output_variable(Y);
        tape->set_derivative(Y, 1.0);
        tape->backpropagate();

        EXPECT_EQUAL(2., Y.get_value());
        EXPECT_EQUAL(0., tape->get_derivative(X));        // should be 0
        EXPECT_EQUAL(0., tape->get_derivative(X2));       // should be 0
        EXPECT_EQUAL(0., tape->get_derivative(Constant)); // should be 0
    }

    {
        adhoc::smart_tape_ptr_t<mode_t> tape;

        type Constant = std::nan("");
        type X = 1.0;
        tape->register_variable(X);
        type X2 = X * 2.;
        type Y = X2;
        tape->register_output_variable(Y);
        tape->set_derivative(Y, 1.0);
        tape->backpropagate();

        EXPECT_EQUAL(2., Y.get_value());
        EXPECT_EQUAL(2., tape->get_derivative(X));
        EXPECT_EQUAL(1., tape->get_derivative(X2));
        EXPECT_EQUAL(0., tape->get_derivative(Constant));
    }

    {
        adhoc::smart_tape_ptr_t<mode_t> tapeptr;
        auto& tape = *tapeptr;
        tape.configure(adhoc::Method::Bwd, 1, 2);

        type Constant = std::nan("");
        type X = 1.0;
        tape.register_variable(X);
        type X2 = X * 2.;
        type Y = X2;
        tape.register_output_variable(Y);
        type Y2 = Y * 2.;
        tape.register_output_variable(Y2);
        tape.set_derivative(Y2, 1.0);
        tape.backpropagate();

        EXPECT_EQUAL(2., Y.get_value());
        EXPECT_EQUAL(4., tape.get_derivative(X));
        EXPECT_EQUAL(2., tape.get_derivative(X2));
        EXPECT_EQUAL(0., tape.get_derivative(Constant));
    }
}

void
register_output_constant()
{
    const std::array<adhoc::Method, 14> AllMethods = {
        adhoc::Method::Bwd,
        adhoc::Method::BwdBuffer,
        adhoc::Method::BwdBufferCompressed,
        adhoc::Method::BwdBufferPathReuse,
        adhoc::Method::BwdBufferCompressedPathReuse,
        adhoc::Method::BwdBufferCompressedPathReuseV,
        adhoc::Method::SecondOrderSimple,
        adhoc::Method::SecondOrderSimd8_stdmap,
        adhoc::Method::SecondOrderSimd8_stdunorderedmap,
        adhoc::Method::SecondOrderSimd8_ankerl,
        adhoc::Method::SecondOrderSimd8_boost,
        adhoc::Method::SecondOrderLossy,
        adhoc::Method::SecondOrderVLossy /* ,
         adhoc::Method::Fwd */
    };
    using mode_t = adhoc::opcode<float>;
    using type = mode_t::type;

    // loop over all adhoc::Method
    for (const auto& method : AllMethods) {
        adhoc::smart_tape_ptr_t<mode_t> tapeptr;
        auto& tape = *tapeptr;
        tape.configure(method, 1, 2);

        type Constant = std::nan("");
        type X = 1.0;
        tape.register_variable(X);
        type X2 = X * 2.;
        type Y = X2;
        type Y2 = 3;
        tape.register_output_variable(Y);
        tape.register_output_variable(Y2);
        tape.set_derivative(Y, 1.0);
        tape.backpropagate();

        EXPECT_EQUAL(2., Y.get_value());
        EXPECT_EQUAL(2., tape.get_derivative(X));
        // EXPECT_EQUAL(1., tape.get_derivative(X2));
        // EXPECT_EQUAL(0., tape.get_derivative(Constant));
    }
}

} // namespace

auto
main() -> int
{
    register_output_constant();
    register_input();
    register_output();

    TEST_END;
}
