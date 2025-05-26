#include "../../../test_simple/test_simple_include.hpp"

#include "../../include/tools/circular_buffer.hpp"

#include <cstdint>
#include <numeric>

int main() {
    adhoc::circular_buffer<std::uint32_t, 10> cb;
    std::iota(cb.data().begin(), cb.data().end(), 5);
    ++cb;
    EXPECT_EQUAL(cb.at(), 6);
    EXPECT_EQUAL(cb.at<9>(), 5);

    adhoc::circular_buffer<std::uint32_t, 10> cb2;
    std::iota(cb2.data().begin(), cb2.data().end(), 6);

    cb2.at<9>() = 5;
    EXPECT_EQUAL((cb == cb2), true);

    cb2.at<9>() = 6;
    EXPECT_EQUAL((cb == cb2), false);

    static_assert(adhoc::is_power_of_two(1), "1 is a power of 2");
    static_assert(adhoc::is_power_of_two(8), "8 is a power of 2");
    static_assert(!adhoc::is_power_of_two(0), "0 is not a power of 2");
    static_assert(!adhoc::is_power_of_two(10), "10 is not a power of 2");

    TEST_END;
}
