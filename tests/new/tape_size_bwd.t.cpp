#include <new/adhoc.hpp>
#include <new/tape_size_bwd.hpp>
#include <new/tape_static.hpp>
#include <new/utils.hpp>

#include <gtest/gtest.h>

namespace adhoc {

TEST(TapeSizeBwd, OneVar) {
    auto tape = CreateTapeInitial<1>();
    auto [val1] = tape.getalias();
    constexpr std::size_t size = detail::tape_size(
        args<decltype(val1) const>{}, args<>{}, args<decltype(val1) const>{});
    static_assert(size == 1);
}

TEST(TapeSizeBwd, TapeIsZero) {
    auto tape = CreateTapeInitial<2>();
    auto [val1, val2] = tape.getalias();
    constexpr std::size_t size = detail::tape_size(
        args<decltype(val1) const>{}, args<>{}, args<decltype(val2) const>{});
    static_assert(size == 0);
}

TEST(TapeSizeBwd, TapeMul) {
    auto tape = CreateTapeInitial<2>();
    auto [val1, val2] = tape.getalias();
    auto valprod = val1 * val2;
    constexpr std::size_t size =
        detail::tape_size(args<decltype(valprod) const>{}, args<>{},
                          args<decltype(val1) const, decltype(val2) const>{});
    static_assert(size == 2);
}

TEST(TapeSizeBwd, TapeDiv) {
    auto tape = CreateTapeInitial<2>();
    auto [val1, val2] = tape.getalias();
    auto result = val1 / val2;

    constexpr std::size_t size = detail::tape_size(
        args<decltype(result) const>{}, args<>{}, args<decltype(val1) const>{});
    static_assert(size == 1);

    constexpr std::size_t size2 =
        detail::tape_size(args<decltype(result) const>{}, args<>{},
                          args<decltype(val1) const, decltype(val2) const>{});
    static_assert(size2 == 2);
}

TEST(TapeSizeBwd, TapeAdd) {
    auto tape = CreateTapeInitial<2>();
    auto [val1, val2] = tape.getalias();
    auto result = val1 + val2;

    constexpr std::size_t size2 =
        detail::tape_size(args<decltype(result) const>{}, args<>{},
                          args<decltype(val1) const, decltype(val2) const>{});
    static_assert(size2 == 2);
}

TEST(TapeSizeBwd, Complex) {
    detail::adouble_aux<0> val1;
    detail::adouble_aux<1> val2;
    detail::adouble_aux<2> val3;
    // auto tape = CreateTapeInitial<3>();
    // auto [val1, val2, val3] = tape.getalias();

    auto valsum = val1 + val2 + val3;
    auto valprod = val1 * val2;
    auto valprod2 = val2 * val2;
    auto valprod3 = valprod * valprod2;
    auto res = valsum * valprod3;

    // constexpr std::size_t size =
    //     detail::tape_size(args<decltype(res) const>{}, args<>{},
    //                       args<decltype(val1) const, decltype(val2) const,
    //                            decltype(val3) const>{});
    constexpr std::size_t size = tape_size(res, val1, val2, val3);
    static_assert(size == 5);

    auto res2 = valprod3 * valsum;
    constexpr std::size_t size2 =
        detail::tape_size(args<decltype(res2) const>{}, args<>{},
                          args<decltype(val1) const, decltype(val2) const,
                               decltype(val3) const>{});
    static_assert(size2 == 4);
}

TEST(TapeSizeBwd, Skip) {
    detail::adouble_aux<0> val1;
    detail::adouble_aux<1> val2;
    detail::adouble_aux<2> val3;
    auto valsum = val1 + val2;
    auto valprod = valsum * val3;
    auto res = valsum * valprod;

    constexpr std::size_t size = tape_size(res, val1, val2, val3);
    static_assert(size == 3);
}

TEST(TapeSizeBwd, ComplexTempVar) {
    detail::adouble_aux<0> val1;
    detail::adouble_aux<1> val2;
    detail::adouble_aux<2> val3;
    auto valsum = val1 + val2 * detail::adouble_aux<3>{} + val3;
    auto valprod = val1 * val2 + detail::adouble_aux<4>{};
    auto valprod2 = val2 * val2 * detail::adouble_aux<5>{};
    auto valprod3 = valprod * valprod2;
    auto res = valsum * valprod3;

    constexpr std::size_t size = tape_size(res, val1, val2, val3);
    static_assert(size == 5);

    auto res2 = valprod3 * valsum;
    constexpr std::size_t size2 = tape_size(res2, val1, val2, val3);
    static_assert(size2 == 4);
}

TEST(TapeSizeBwd, MulTemp) {
    detail::adouble_aux<0> val1;
    auto m1 = val1 * detail::adouble_aux<1>{};
    auto m2 = val1 * detail::adouble_aux<2>{};
    auto valprod = m1 * m2;

    constexpr std::size_t size = tape_size(valprod, val1);
    static_assert(size == 2);
}

TEST(TapeSizeBwd, MulMul) {
    detail::adouble_aux<0> val1;
    detail::adouble_aux<1> val2;
    auto m = val1 * val2;
    auto valprod = m * m;

    constexpr std::size_t size = tape_size(valprod, val1, val2);
    static_assert(size == 2);
}

} // namespace adhoc