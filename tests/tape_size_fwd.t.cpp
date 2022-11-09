#include <adhoc2.hpp>
#include <tape_size_fwd.hpp>

#include <gtest/gtest.h>

namespace adhoc2 {

TEST(adhoc2, tapesize_fwd) {
    adouble val1(1.);
    constexpr std::size_t size = tape_size_fwd(val1, val1);
    static_assert(size == 1);
}

TEST(adhoc2, tapesizefwdzero) {
    adouble val1(1.);
    adouble val2(1.);
    constexpr std::size_t size = tape_size_fwd(val1, val2);
    static_assert(size == 0);
}

TEST(adhoc2, tapefwdbi) {
    adouble val1(1.);
    adouble val2(2.);
    auto valprod = val1 * val2;
    constexpr std::size_t size = tape_size_fwd(valprod, val1, val2);
    static_assert(size == 3);
}

TEST(adhoc2, tapefwddiv) {
    adouble val1(2.0);
    adouble val2(3.0);
    auto result = val1 / val2;
    constexpr std::size_t size = tape_size_fwd(result, val2);
    static_assert(size == 2);

    constexpr std::size_t size2 = tape_size_fwd(result, val1, val2);
    static_assert(size2 == 3);
}

} // namespace adhoc2