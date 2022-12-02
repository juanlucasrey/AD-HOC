#include <init.hpp>
#include <tape_size_bwd.hpp>

#include <gtest/gtest.h>

namespace adhoc2 {

TEST(TapeSizeBwd, OneVar) {
    auto [val1] = Init<1>();
    constexpr std::size_t size = tape_size_bwd<decltype(val1)>();
    static_assert(size == 0);
}

TEST(TapeSizeBwd, TapeIsZero) {
    auto [val1, val2] = Init<2>();
    constexpr std::size_t size =
        tape_size_bwd<decltype(val1), decltype(val2)>();
    static_assert(size == 0);
}

TEST(TapeSizeBwd, TapeMul) {
    auto [val1, val2] = Init<2>();
    auto result = val1 * val2;
    constexpr std::size_t size =
        tape_size_bwd<decltype(result), decltype(val1), decltype(val2)>();
    static_assert(size == 0);
}

TEST(TapeSizeBwd, TapeDiv) {
    auto [val1, val2] = Init<2>();
    auto result = val1 / val2;
    constexpr std::size_t size =
        tape_size_bwd<decltype(result), decltype(val1), decltype(val2)>();
    static_assert(size == 0);
}

TEST(TapeSizeBwd, TapeAdd) {
    auto [val1, val2] = Init<2>();
    auto result = val1 + val2;
    constexpr std::size_t size =
        tape_size_bwd<decltype(result), decltype(val1), decltype(val2)>();
    static_assert(size == 0);
}

TEST(TapeSizeBwd, MulTemp) {
    auto [val1, val2, val3] = Init<3>();
    auto m1 = val1 * val2;
    auto m2 = val1 * val3;
    auto result = m1 * m2;

    constexpr std::size_t size =
        tape_size_bwd<decltype(result), decltype(val1)>();
    static_assert(size == 2);
}

TEST(TapeSizeBwd, MulMul) {
    auto [val1, val2] = Init<2>();
    auto m = val1 * val2;
    auto result = m * m;

    constexpr std::size_t size =
        tape_size_bwd<decltype(result), decltype(val1), decltype(val2)>();
    static_assert(size == 1);
}

TEST(TapeSizeBwd, Complex) {
    auto [val1, val2, val3] = Init<3>();
    auto valsum = val1 + val2 + val3;
    auto valprod = val1 * val2;
    auto valprod2 = val2 * val2;
    auto valprod3 = valprod * valprod2;
    auto res = valsum * valprod3;

    constexpr std::size_t size =
        tape_size_bwd<decltype(res), decltype(val1), decltype(val2),
                      decltype(val3)>();
    static_assert(size == 2);
}

TEST(TapeSizeBwd, Skip) {
    auto [val1, val2, val3] = Init<3>();
    auto valsum = val1 + val2;
    auto valprod = valsum * val3;
    auto res = valsum * valprod;

    constexpr std::size_t size =
        tape_size_bwd<decltype(res), decltype(val1), decltype(val2),
                      decltype(val3)>();
    static_assert(size == 2);
}

TEST(TapeSizeBwd, SameCalcDifferentTapeSize) {
    auto [val1, val2, val3, val4, val5, val6] = Init<6>();
    auto valsum = val1 + val2 * val4 + val3;
    auto valprod = val1 * val2 + val5;
    auto valprod2 = val2 * val2 * val6;
    auto valprod3 = valprod * valprod2;
    auto res = valsum * valprod3;

    constexpr std::size_t size =
        tape_size_bwd<decltype(res), decltype(val1), decltype(val2),
                      decltype(val3)>();
    static_assert(size == 2);

    auto res2 = valprod3 * valsum;
    constexpr std::size_t size2 =
        tape_size_bwd<decltype(res2), decltype(val1), decltype(val2),
                      decltype(val3)>();
    static_assert(size2 == 3);
}

} // namespace adhoc2
