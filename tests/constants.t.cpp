#include <constants.hpp>

#include <adhoc2.hpp>

#include <gtest/gtest.h>

namespace adhoc2::constants2 {

TEST(constants, frac) {
    std::cout.precision(std::numeric_limits<double>::max_digits10);
    constexpr double one_over_root_two = Frac<Const<1>, Sqrt<Const<2>>>::v();
    static_assert(one_over_root_two == 0.70710678118654757);
}

TEST(constants, fracop) {
    std::cout.precision(std::numeric_limits<double>::max_digits10);
    constexpr double one_over_root_two = Frac<Const<1>, Sqrt<Const<2>>>::v();
    static_assert(one_over_root_two == 0.70710678118654757);

    double S = 100;
    adhoc<ID> aS(S);
    auto result = aS * Frac<Const<1>, Sqrt<Const<2>>>();
}

TEST(constants, muldouble) {
    std::cout.precision(std::numeric_limits<double>::max_digits10);
    constexpr double one_over_root_two = Frac<Const<1>, Sqrt<Const<2>>>::v();
    static_assert(one_over_root_two == 0.70710678118654757);

    double S = 100;
    double result = S * Frac<Const<1>, Sqrt<Const<2>>>();
}

} // namespace adhoc2::constants2
