#include <differential_operator/sort.hpp>
#include <init.hpp>
#include <tape2.hpp>

#include <gtest/gtest.h>
#include <tuple>

#include "../type_name.hpp"

namespace adhoc3 {

TEST(DifferentialOperator, OrderedOperators) {
    auto [x1, x2] = Init<2>();
    auto res = x1 * x2;
    CalcTree t(res);
    // we create some calculation order type
    decltype(t)::ValuesTupleInverse co;

    auto d1 = d<2>(x1) * pow<3>(d(x2));
    constexpr auto d1_ordered = order_differential_operator(d1, co);

    auto d2 = d<2>(x2) * pow<3>(d(x1));
    constexpr auto d2_ordered = order_differential_operator(d2, co);

    auto ops = std::make_tuple(d1, d2);
    auto opres = order_differential_operators(ops, co);
    std::cout << type_name2<decltype(opres)>() << std::endl;

    auto ops2 = std::make_tuple(d2, d1);
    auto opres2 = order_differential_operators(ops2, co);
    std::cout << type_name2<decltype(opres2)>() << std::endl;
}

} // namespace adhoc3
