#include <adhoc2.hpp>
#include <tape_size.hpp>
#include <utils.hpp>

#include <gtest/gtest.h>

namespace adhoc2 {

TEST(adhoc2, tapesize) {
    adouble val1(1.);
    constexpr std::size_t size = tape_size(val1, val1);
    static_assert(size == 1);
}

TEST(adhoc2, tapemult) {
    adouble val1(1.);
    adouble val2(2.);
    auto valprod = val1 * val2;
    constexpr std::size_t size = tape_size(valprod, val1, val2);
    static_assert(size == 2);
}

TEST(adhoc2, tapeadd) {
    adouble val1(1.);
    adouble val2(2.);
    auto valprod = val1 + val2;
    constexpr std::size_t size = tape_size(valprod, val1, val2);
    static_assert(size == 2);
}

TEST(adhoc2, tapecomplexd) {
    adouble val1(1.);
    adouble val2(2.);
    adouble val3(2.);
    auto valsum = val1 + val2 + val3;
    auto valprod = val1 * val2;
    auto valprod2 = val2 * val2;
    auto valprod3 = valprod * valprod2;
    auto res = valsum * valprod3;

    constexpr std::size_t size = tape_size(res, val1, val2, val3);
    static_assert(size == 5);

    auto res2 = valprod3 * valsum;
    constexpr std::size_t size2 = tape_size(res2, val1, val2, val3);
    static_assert(size2 == 4);
}

TEST(adhoc2, tapeskip) {
    adouble val1(1.);
    adouble val2(2.);
    adouble val3(2.);
    auto valsum = val1 + val2;
    auto valprod = valsum * val3;
    auto res = valsum * valprod;

    constexpr std::size_t size = tape_size(res, val1, val2, val3);
    static_assert(size == 3);
}

TEST(adhoc2, tapecomplexd2) {
    adouble val1(1.);
    adouble val2(2.);
    adouble val3(2.);
    auto valsum = val1 + mulscalar(val2, 2) + val3;
    auto valprod = val1 * val2 + 3;
    auto valprod2 = val2 * mulscalar(val2, 4);
    auto valprod3 = valprod * valprod2;
    auto res = valsum * valprod3;

    constexpr std::size_t size = tape_size(res, val1, val2, val3);
    static_assert(size == 5);

    auto res2 = valprod3 * valsum;
    constexpr std::size_t size2 = tape_size(res2, val1, val2, val3);
    static_assert(size2 == 4);
}

TEST(adhoc2, tapetrip) {
    adouble val1(1.);
    auto m1 = mulscalar(val1, 2.0);
    auto m2 = mulscalar(val1, 3.0);
    auto valprod = m1 * m2;

    std::size_t size = tape_size(valprod, val1);
    // static_assert(size == 2);
}

TEST(adhoc2, tapetrip2) {
    adouble val1(1.);
    adouble val2(2.);
    auto m = val1 * val2;
    auto valprod = m * m;

    constexpr std::size_t size = tape_size(valprod, val1, val2);
    static_assert(size == 3);
}

} // namespace adhoc2