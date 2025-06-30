#include "../../../test_simple/test_simple_include.hpp"

#include <tools/circular_buffer.hpp>

#include <cstdint>
#include <numeric>

int main() {
    adhoc::circular_buffer<std::uint32_t, 10> cb;
    std::iota(cb.data().begin(), cb.data().end(), 5);
    ++cb;
    EXPECT_EQUAL(cb.at(), 6);
    EXPECT_EQUAL(cb.template at<9>(), 5);

    adhoc::circular_buffer<std::uint32_t, 10> cb2;
    std::iota(cb2.data().begin(), cb2.data().end(), 6);

    cb2.at<9>() = 5;
    EXPECT_EQUAL((cb == cb2), true);

    cb2.at<9>() = 6;
    EXPECT_EQUAL((cb == cb2), false);

    TEST_END;
}
