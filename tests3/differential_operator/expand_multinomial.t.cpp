#include <differential_operator/expand.hpp>
#include <differential_operator/expand_multinomial.hpp>

#include "../type_name.hpp"

#include <adhoc.hpp>
#include <constants_type.hpp>
#include <init.hpp>
#include <tape2.hpp>

#include <gtest/gtest.h>
#include <tuple>

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

TEST(DifferentialOperator, ExpandUnivariate2) {

    auto [x1] = Init<1>();
    auto res = cos(x1);
    CalcTree t(res);
    // we create some calculation order type
    decltype(t)::ValuesTupleInverse co;

    auto d1 = d<3>(res) * d(x1);
    is_ordered(co, d1);
    auto outders = std::tuple<decltype(d1)>{};

    auto result = expand(co, outders);
    std::cout << type_name2<decltype(result)>() << std::endl;
}

TEST(DifferentialOperator, ExpandUnivariate3) {

    auto [x1] = Init<1>();
    auto res = cos(x1);
    auto res2 = exp(x1);
    CalcTree t(res, res2);
    // we create some calculation order type
    decltype(t)::ValuesTupleInverse co;

    auto d1 = d<3>(res) * d(x1);
    auto d2 = d(res2);
    is_ordered(co, d1);
    auto outders = std::tuple<decltype(d1), decltype(d2)>{};
    is_ordered(co, outders);

    auto result = expand(co, outders);
    std::cout << type_name2<decltype(result)>() << std::endl;
}

TEST(DifferentialOperator, ExpandSumRightConst) {

    using constants::CD;
    using constants::encode;
    auto [x1] = Init<1>();
    auto res = x1 + CD<encode(0.5)>();
    CalcTree t(res);
    decltype(t)::ValuesTupleInverse co;
    auto d1 = d(res);
    auto d2 = d<2>(res);
    auto d3 = d<3>(res);
    auto outders = std::tuple<decltype(d3), decltype(d2), decltype(d1)>{};
    auto result = expand(co, outders);

    auto res2 = x1;
    CalcTree t2(res2);
    decltype(t2)::ValuesTupleInverse co2;
    auto d12 = d(res2);
    auto d22 = d<2>(res2);
    auto d32 = d<3>(res2);
    auto outders2 = std::tuple<decltype(d32), decltype(d22), decltype(d12)>{};
    auto result2 = expand(co2, outders2);
    static_assert(std::is_same_v<decltype(result), decltype(result2)>);
}

TEST(DifferentialOperator, ExpandSumLeftConst) {

    using constants::CD;
    using constants::encode;
    auto [x1] = Init<1>();
    auto res = CD<encode(0.5)>() + x1;
    CalcTree t(res);
    decltype(t)::ValuesTupleInverse co;
    auto d1 = d(res);
    auto d2 = d<2>(res);
    auto d3 = d<3>(res);
    auto outders = std::tuple<decltype(d3), decltype(d2), decltype(d1)>{};
    auto result = expand(co, outders);

    auto res2 = x1;
    CalcTree t2(res2);
    decltype(t2)::ValuesTupleInverse co2;
    auto d12 = d(res2);
    auto d22 = d<2>(res2);
    auto d32 = d<3>(res2);
    auto outders2 = std::tuple<decltype(d32), decltype(d22), decltype(d12)>{};
    auto result2 = expand(co2, outders2);
    static_assert(std::is_same_v<decltype(result), decltype(result2)>);
}

TEST(DifferentialOperator, ExpandSubRightConst) {

    using constants::CD;
    using constants::encode;
    auto [x1] = Init<1>();
    auto res = x1 - CD<encode(0.5)>();
    CalcTree t(res);
    decltype(t)::ValuesTupleInverse co;
    auto d1 = d(res);
    auto d2 = d<2>(res);
    auto d3 = d<3>(res);
    auto outders = std::tuple<decltype(d3), decltype(d2), decltype(d1)>{};
    auto result = expand(co, outders);

    auto res2 = x1;
    CalcTree t2(res2);
    decltype(t2)::ValuesTupleInverse co2;
    auto d12 = d(res2);
    auto d22 = d<2>(res2);
    auto d32 = d<3>(res2);
    auto outders2 = std::tuple<decltype(d32), decltype(d22), decltype(d12)>{};
    auto result2 = expand(co2, outders2);
    static_assert(std::is_same_v<decltype(result), decltype(result2)>);
}

TEST(DifferentialOperator, ExpandSubLeftConst) {

    using constants::CD;
    using constants::encode;
    auto [x1] = Init<1>();
    auto res = CD<encode(0.5)>() - x1;
    CalcTree t(res);
    decltype(t)::ValuesTupleInverse co;
    auto d1 = d(res);
    auto d2 = d<2>(res);
    auto d3 = d<3>(res);
    auto outders = std::tuple<decltype(d3), decltype(d2), decltype(d1)>{};
    auto result = expand(co, outders);

    auto res2 = x1;
    CalcTree t2(res2);
    decltype(t2)::ValuesTupleInverse co2;
    auto d12 = d(res2);
    auto d22 = d<2>(res2);
    auto d32 = d<3>(res2);
    auto outders2 = std::tuple<decltype(d32), decltype(d22), decltype(d12)>{};
    auto result2 = expand(co2, outders2);
    static_assert(std::is_same_v<decltype(result), decltype(result2)>);
}

TEST(DifferentialOperator, ExpandMulRightConst) {

    using constants::CD;
    using constants::encode;
    auto [x1] = Init<1>();
    auto res = x1 * CD<encode(0.5)>();
    CalcTree t(res);
    decltype(t)::ValuesTupleInverse co;
    auto d1 = d(res);
    auto d2 = d<2>(res);
    auto d3 = d<3>(res);
    auto outders = std::tuple<decltype(d3), decltype(d2), decltype(d1)>{};
    auto result = expand(co, outders);

    auto res2 = x1;
    CalcTree t2(res2);
    decltype(t2)::ValuesTupleInverse co2;
    auto d12 = d(res2);
    auto d22 = d<2>(res2);
    auto d32 = d<3>(res2);
    auto outders2 = std::tuple<decltype(d32), decltype(d22), decltype(d12)>{};
    auto result2 = expand(co2, outders2);
    static_assert(std::is_same_v<decltype(result), decltype(result2)>);
}

TEST(DifferentialOperator, ExpandMulLeftConst) {

    using constants::CD;
    using constants::encode;
    auto [x1] = Init<1>();
    auto res = CD<encode(0.5)>() * x1;
    CalcTree t(res);
    decltype(t)::ValuesTupleInverse co;
    auto d1 = d(res);
    auto d2 = d<2>(res);
    auto d3 = d<3>(res);
    auto outders = std::tuple<decltype(d3), decltype(d2), decltype(d1)>{};
    auto result = expand(co, outders);

    auto res2 = x1;
    CalcTree t2(res2);
    decltype(t2)::ValuesTupleInverse co2;
    auto d12 = d(res2);
    auto d22 = d<2>(res2);
    auto d32 = d<3>(res2);
    auto outders2 = std::tuple<decltype(d32), decltype(d22), decltype(d12)>{};
    auto result2 = expand(co2, outders2);
    static_assert(std::is_same_v<decltype(result), decltype(result2)>);
}

TEST(DifferentialOperator, ExpandDivRightConst) {

    using constants::CD;
    using constants::encode;
    auto [x1] = Init<1>();
    auto res = x1 / CD<encode(0.5)>();
    CalcTree t(res);
    // we create some calculation order type
    decltype(t)::ValuesTupleInverse co;
    auto d1 = d(res);
    auto d2 = d<2>(res);
    auto d3 = d<3>(res);
    auto outders = std::tuple<decltype(d3), decltype(d2), decltype(d1)>{};
    auto result = expand(co, outders);

    auto res2 = x1;
    CalcTree t2(res2);
    decltype(t2)::ValuesTupleInverse co2;
    auto d12 = d(res2);
    auto d22 = d<2>(res2);
    auto d32 = d<3>(res2);
    auto outders2 = std::tuple<decltype(d32), decltype(d22), decltype(d12)>{};
    auto result2 = expand(co2, outders2);
    static_assert(std::is_same_v<decltype(result), decltype(result2)>);
}

TEST(DifferentialOperator, ExpandDivLeftConst) {

    using constants::CD;
    using constants::encode;
    auto [x1] = Init<1>();
    auto res = CD<encode(0.5)>() / x1;
    CalcTree t(res);
    // we create some calculation order type
    decltype(t)::ValuesTupleInverse co;
    auto d1 = d(res);
    auto d2 = d<2>(res);
    auto d3 = d<3>(res);
    auto outders = std::tuple<decltype(d3), decltype(d2), decltype(d1)>{};
    auto result = expand(co, outders);

    auto res2 = exp(x1);
    CalcTree t2(res2);
    decltype(t2)::ValuesTupleInverse co2;
    auto d12 = d(res2);
    auto d22 = d<2>(res2);
    auto d32 = d<3>(res2);
    auto outders2 = std::tuple<decltype(d32), decltype(d22), decltype(d12)>{};
    auto result2 = expand(co2, outders2);
    static_assert(std::is_same_v<decltype(result), decltype(result2)>);
}

TEST(DifferentialOperator, ExpandSum) {
    auto [x1, x2] = Init<2>();
    auto res = x1 + x2;
    CalcTree t(res);
    decltype(t)::ValuesTupleInverse co;
    auto d1 = d(res);
    auto d2 = d<3>(res);
    auto outders = std::tuple<decltype(d2), decltype(d1)>{};
    auto result = expand(co, outders);

    auto result2 = std::make_tuple(d<3>(x1), d(x1), d<3>(x2), d(x2));
    static_assert(std::is_same_v<decltype(result), decltype(result2)>);
}

TEST(DifferentialOperator, ExpandMult) {
    auto [x1, x2] = Init<2>();
    auto res = x1 * x2;
    CalcTree t(res);
    decltype(t)::ValuesTupleInverse co;
    auto d1 = d(res);
    auto d2 = d<3>(res);
    auto outders = std::tuple<decltype(d2), decltype(d1)>{};
    auto result = expand(co, outders);
    auto result2 = std::make_tuple(d<3>(x1), d<2>(x1) * d(x2), d(x1) * d<2>(x2),
                                   d(x1), d<3>(x2), d(x2));
    static_assert(std::is_same_v<decltype(result), decltype(result2)>);
}

TEST(DifferentialOperator, ExpandDiv) {
    auto [x1, x2] = Init<2>();
    auto res = x1 / x2;
    CalcTree t(res);
    decltype(t)::ValuesTupleInverse co;
    auto d1 = d(res);
    auto d2 = d<3>(res);
    auto outders = std::tuple<decltype(d2), decltype(d1)>{};
    auto result = expand(co, outders);

    auto res2 = x1 * exp(x2);
    CalcTree t2(res2);
    decltype(t2)::ValuesTupleInverse co2;
    auto d12 = d(res2);
    auto d22 = d<3>(res2);
    auto outders2 = std::tuple<decltype(d22), decltype(d12)>{};
    auto result2 = expand(co2, outders2);
    static_assert(std::is_same_v<decltype(result), decltype(result2)>);
}

} // namespace adhoc3
