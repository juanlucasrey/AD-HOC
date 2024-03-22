#include <differential_operator/order_differential_operator.hpp>
#include <differential_operator/ordered_mult.hpp>
#include <init.hpp>
#include <tape2.hpp>

#include <gtest/gtest.h>

namespace adhoc3 {

TEST(DifferentialOperator, OrderedMultipicationMixed) {
    auto [x1, x2] = Init<2>();
    auto res = x1 * x2;
    CalcTree t(res);
    // we create some calculation order type
    decltype(t)::ValuesTupleInverse co;

    auto d1 = d<2>(x1) * pow<3>(d(x2));
    constexpr auto d1_ordered = order_differential_operator(d1, co);

    auto d2 = d<2>(x2) * pow<3>(d(x1));
    constexpr auto d2_ordered = order_differential_operator(d2, co);

    auto d1_times_d2 = multiply_ordered(co, d1_ordered, d2_ordered);

    auto result = d<2>(x1) * pow<3>(d(x1)) * d<2>(x2) * pow<3>(d(x2));
    static_assert(std::is_same_v<decltype(result), decltype(d1_times_d2)>);
}

TEST(DifferentialOperator, OrderedMultipicationSame) {
    auto [x1, x2] = Init<2>();
    auto res = x1 * x2;
    CalcTree t(res);
    // we create some calculation order type
    decltype(t)::ValuesTupleInverse co;

    auto d1 = d<2>(x1) * pow<3>(d(x2));
    constexpr auto d1_ordered = order_differential_operator(d1, co);

    auto d1_times_d1 = multiply_ordered(co, d1_ordered, d1_ordered);

    auto result = pow<2>(d1_ordered);
    static_assert(std::is_same_v<decltype(result), decltype(d1_times_d1)>);
}

TEST(DifferentialOperator, OrderedMultipicationEmptyRight) {
    auto [x1, x2] = Init<2>();
    auto res = x1 * x2;
    CalcTree t(res);
    // we create some calculation order type
    decltype(t)::ValuesTupleInverse co;

    auto d1 = d<2>(x1) * pow<3>(d(x2));
    auto d1_ordered = order_differential_operator(d1, co);

    auto d1_times_d1 = multiply_ordered(co, d1_ordered, std::tuple<>{});

    static_assert(std::is_same_v<decltype(d1_ordered), decltype(d1_times_d1)>);
}

TEST(DifferentialOperator, OrderedMultipicationEmptyLeft) {
    auto [x1, x2] = Init<2>();
    auto res = x1 * x2;
    CalcTree t(res);
    // we create some calculation order type
    decltype(t)::ValuesTupleInverse co;

    auto d1 = d<2>(x1) * pow<3>(d(x2));
    auto d1_ordered = order_differential_operator(d1, co);

    auto d1_times_d1 = multiply_ordered(co, std::tuple<>{}, d1_ordered);

    static_assert(std::is_same_v<decltype(d1_ordered), decltype(d1_times_d1)>);
}

TEST(DifferentialOperator, OrderedMultipicationMultiple) {
    auto [x1, x2] = Init<2>();
    auto res = x1 * x2;
    CalcTree t(res);
    // we create some calculation order type
    decltype(t)::ValuesTupleInverse co;

    auto d1 = d<2>(x1) * pow<3>(d(x2));
    constexpr auto d1_ordered = order_differential_operator(d1, co);

    auto d2 = d<2>(x2) * pow<3>(d(x1));
    constexpr auto d2_ordered = order_differential_operator(d2, co);

    auto d3 = d(x1);
    constexpr auto d3_ordered = order_differential_operator(d3, co);

    auto empty = std::tuple<>{};

    auto d1_times_d1 =
        multiply_ordered(co, d1_ordered, d2_ordered, d3_ordered, empty);

    auto result = d<2>(x1) * pow<4>(d(x1)) * d<2>(x2) * pow<3>(d(x2));

    static_assert(std::is_same_v<decltype(d1_times_d1), decltype(result)>);
}

} // namespace adhoc3
