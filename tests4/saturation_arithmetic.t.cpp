#include <adhoc.hpp>
#include <base.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <numeric>

namespace adhoc4 {

namespace {

template <typename T, typename... Types>
concept is_all_same = (... && std::is_same<T, Types>::value);

constexpr auto add_sat(std::size_t lhs) { return lhs; }

template <is_all_same<std::size_t>... Types>
constexpr auto add_sat(std::size_t lhs, Types... args) {
    auto const rhs = add_sat(args...);

    if (std::numeric_limits<std::size_t>::max() - lhs < rhs) {
        return std::numeric_limits<std::size_t>::max();
    }

    if (std::numeric_limits<std::size_t>::max() - rhs < lhs) {
        return std::numeric_limits<std::size_t>::max();
    }

    return lhs + rhs;
}

} // namespace

TEST(SaturationArithmetic, Inverse) {

    constexpr auto max = std::numeric_limits<std::size_t>::max();
    constexpr auto min = std::numeric_limits<std::size_t>::min();

    constexpr auto a =
        add_sat(static_cast<std::size_t>(3), static_cast<std::size_t>(4));
    static_assert(a == 7);

    constexpr auto b = add_sat(max, static_cast<std::size_t>(4));
    static_assert(b == max);

    constexpr auto c = add_sat(max, max, max);
    static_assert(c == max);

    constexpr auto d =
        add_sat(static_cast<std::size_t>(3), static_cast<std::size_t>(3),
                static_cast<std::size_t>(3));
    static_assert(d == 9);
}

} // namespace adhoc4
