#include <tape.hpp>
#include <test_simple_include.hpp>

void
test_tape_no_default()
{
    using adhoc_mode2 = adhoc::
      opcode<double, "main", double, adhoc::EnumVectorType::Valuecompression, adhoc::IdxVectorType::BitCompression>;
    using adhoc_t2 = adhoc_mode2::type;

    adhoc::smart_tape_ptr_t<adhoc_mode2> tapeptr;
    auto& tape = *tapeptr;
    tape.configure(adhoc::Method::Bwd, 2, 1);
    adhoc_t2 x1_adhoc = 3.0;
    adhoc_t2 x2_adhoc = 5.0;
    tape.register_variable(x1_adhoc);
    tape.register_variable(x2_adhoc);
    adhoc_t2 y_adhoc = x1_adhoc + x2_adhoc;
    tape.register_output_variable(y_adhoc);
    tape.set_derivative(y_adhoc, 1.0);
    tape.backpropagate();
    double dy_dx1_adhoc = tape.get_derivative(x1_adhoc);
    double dy_dx2_adhoc = tape.get_derivative(x2_adhoc);

    EXPECT_NEAR_ABS(dy_dx1_adhoc, 1., 1e-10);
    EXPECT_NEAR_ABS(dy_dx2_adhoc, 1., 1e-10);
}

auto
main() -> int
{
    test_tape_no_default();

    TEST_END;
}
