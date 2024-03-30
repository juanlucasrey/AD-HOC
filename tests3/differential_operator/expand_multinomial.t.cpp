#include <differential_operator/expand_multinomial.hpp>

#include <init.hpp>
#include <tape2.hpp>

#include <gtest/gtest.h>

namespace adhoc3 {

TEST(DifferentialOperator, ExpandMultinomial) {
    auto [x1, x2] = Init<2>();
    auto res = x1 * x2;
    CalcTree t(res);
    // we create some calculation order type
    decltype(t)::ValuesTupleInverse co;

    auto d1 = d(x1);
    auto d2 = d<2>(res);
    auto d3 = d<3>(x2);

    constexpr auto deriv_tuple =
        std::tuple<decltype(d2), decltype(d1), decltype(d3)>{};
    constexpr auto result = expand_multinomial<3>(co, deriv_tuple);

    constexpr auto result2 =
        std::tuple<decltype(pow<3>(d2)), decltype(pow<2>(d2) * d1),
                   decltype(pow<2>(d2) * d3), decltype(d2 * pow<2>(d1)),
                   decltype(d2 * d1 * d3), decltype(d2 * pow<2>(d3)),
                   decltype(pow<3>(d1)), decltype(pow<2>(d1) * d3),
                   decltype(d1 * pow<2>(d3)), decltype(pow<3>(d3))>{};
    static_assert(std::is_same_v<decltype(result), decltype(result2)>);
}

} // namespace adhoc3
