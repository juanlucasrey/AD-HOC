#include <adhoc2.hpp>
#include <tape_size.hpp>
#include <utils.hpp>

#include <gtest/gtest.h>

namespace adhoc2 {

TEST(adhoc2, tapesize4) {
    adouble val1(1.);
    constexpr std::size_t size = detail::tape_size(
        args<decltype(val1) const>{}, args<>{}, args<decltype(val1)>{});
    static_assert(size == 1);
}

TEST(adhoc2, tape4mult) {
    adouble val1(1.);
    adouble val2(2.);
    auto valprod = val1 * val2;
    constexpr std::size_t size =
        detail::tape_size(args<decltype(valprod) const>{}, args<>{},
                          args<decltype(val1), decltype(val2)>{});
    static_assert(size == 2);
}

TEST(adhoc2, tape4add) {
    adouble val1(1.);
    adouble val2(2.);
    auto valprod = val1 + val2;
    constexpr std::size_t size =
        detail::tape_size(args<decltype(valprod) const>{}, args<>{},
                          args<decltype(val1), decltype(val2)>{});
    static_assert(size == 2);
}

TEST(adhoc2, tape4complexd) {
    adouble val1(1.);
    adouble val2(2.);
    adouble val3(2.);
    auto valsum = val1 + val2 + val3;
    auto valprod = val1 * val2;
    auto valprod2 = val2 * val2;
    auto valprod3 = valprod * valprod2;
    auto res = valsum * valprod3;

    constexpr std::size_t size = detail::tape_size(
        args<decltype(res) const>{}, args<>{},
        args<decltype(val1), decltype(val2), decltype(val3)>{});
    static_assert(size == 5);

    auto res2 = valprod3 * valsum;
    constexpr std::size_t size2 = detail::tape_size(
        args<decltype(res2) const>{}, args<>{},
        args<decltype(val1), decltype(val2), decltype(val3)>{});
    static_assert(size2 == 4);
}

TEST(adhoc2, tape4) {
    adouble val1(1.);
    constexpr std::size_t size = detail::tape_size(
        args<decltype(val1) const>{}, args<>{}, args<decltype(val1)>{});
    static_assert(size == 1);
}

} // namespace adhoc2