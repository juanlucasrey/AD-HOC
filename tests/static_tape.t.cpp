#include <tape_static.hpp>

#include <gtest/gtest.h>

namespace adhoc {

TEST(StaticTape, Input) {
    auto tape = CreateTapeInitial<4>();
    auto [v0, v1, v2, v3] = tape.getalias();

    tape.set(v0, 1.0);
    tape.set(v1, 2.0);
    tape.set(v2, 3.0);
    tape.set(v3, 4.0);

    detail::adouble_aux<3> temp;
    tape.set(temp, 5.0);
}

} // namespace adhoc
