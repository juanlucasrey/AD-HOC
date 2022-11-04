#include <tape_input.hpp>

#include <gtest/gtest.h>

namespace adhoc2 {

TEST(adhoc2, tapestatic2) {
    adhoc<ID> v1(1);
    adhoc<ID> v2(2);
    adhoc<ID> v3(3);
    adhoc<ID> v4(4);
    auto mul = v1 * v2;
    auto tape = TapeInput(v1, v2);

    auto tape2 = TapeInput(v1);
    // auto tape3 = TapeInput(v1, v1); // compilation error
    // auto tape4 = TapeInput(v1, mul); // compilation error
    // auto tape5 = TapeInput(mul, v1); // compilation error

    tape.set(v1, 2.);
    tape.set(v2, 3.);
    // tape.set(v3, 2.); // compilation error

    // addvar(tape, v3); // compilation error, no discard
    auto tapeb = addvar(tape, v3);
    tapeb.set(v3, 4.);
}

} // namespace adhoc2
