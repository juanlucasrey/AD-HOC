#include <differential_operator/expand_univariate.hpp>

#include <init.hpp>
#include <tape2.hpp>

#include <gtest/gtest.h>

namespace adhoc3 {

TEST(DifferentialOperator, ExpandUnivariate) {
    auto [x1, x2] = Init<2>();
    auto res = x1 * x2;
    CalcTree t(res);
    // we create some calculation order type
    decltype(t)::ValuesTupleInverse co;

    constexpr auto result = expand_univariate<decltype(x1), 5>(co);
    constexpr auto result2 = std::tuple<
        decltype(d<5>(x1)), decltype(d<4>(x1) * d(x1)),
        decltype(d<3>(x1) * d<2>(x1)), decltype(d<3>(x1) * pow<2>(d(x1))),
        decltype(pow<2>(d<2>(x1)) * d(x1)), decltype(d<2>(x1) * pow<3>(d(x1))),
        decltype(pow<5>(d(x1)))>{};
    static_assert(std::is_same_v<decltype(result), decltype(result2)>);
}

} // namespace adhoc3
