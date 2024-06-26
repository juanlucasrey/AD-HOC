#include <differential_operator/derivative_non_null.hpp>
#include <differential_operator/differential_operator.hpp>
#include <differential_operator/is_sorted.hpp>
#include <differential_operator/less_than.hpp>
#include <differential_operator/merge_ordered.hpp>
#include <differential_operator/select_root_derivatives.hpp>
#include <differential_operator/sort.hpp>
// #include <differential_operator_util.hpp>
#include <init.hpp>
#include <tape2.hpp>

#include "type_name.hpp"

#include <gtest/gtest.h>

namespace adhoc3 {

TEST(DifferentialOperator2, dID) {

    auto [x, y] = Init<2>();
    auto res = (x * y) * (x * cos(y));

    auto d1 = d<2>(x);
    auto d2 = pow<2>(d1);
    auto d3 = pow<2>(d<2>(x));
    static_assert(std::is_same_v<decltype(d2), decltype(d3)>);
}

TEST(DifferentialOperator2, dID2) {

    auto [x, y] = Init<2>();
    auto res = (x * y) * (x * cos(y));

    auto d1 = d<2>(x) * d<2>(x);
    auto d2 = d<2>(x) * d(x);
    auto d3 = d1 * d2;
    auto d4 = pow<3>(d<2>(x)) * d(x);
    static_assert(std::is_same_v<decltype(d3), decltype(d4)>);
}

TEST(DifferentialOperator2, isRoot) {
    auto [x, y] = Init<2>();
    auto res = (x * y) * (x * cos(y));

    auto dx = d(x);
    static_assert(is_root_class_v<decltype(dx)>);
    auto dy = d(y);
    static_assert(is_root_class_v<decltype(dy)>);
    auto dcross = d(y) * d(x);
    static_assert(is_root_class_v<decltype(dcross)>);
    auto dres = d(res);
    static_assert(!is_root_class_v<decltype(dres)>);
    auto dcrosswithres = d(res) * d(x);
    static_assert(!is_root_class_v<decltype(dcrosswithres)>);
}

TEST(DifferentialOperator2, selectRoot) {
    auto [x, y] = Init<2>();
    auto res = (x * y) * (x * cos(y));

    auto dx = d(x);
    auto dy = d(y);
    auto dcross = d(y) * d(x);
    auto dres = d(res);
    auto dcrosswithres = d(res) * d(x);
    auto tu = std::make_tuple(dx, dy, dcross, dres, dcrosswithres);
    auto tus = select_root_derivatives(tu);
    auto tures = std::make_tuple(dx, dy, dcross);
    static_assert(std::is_same_v<decltype(tus), decltype(tures)>);
    auto tuns = select_non_root_derivatives(tu);
    auto turens = std::make_tuple(dres, dcrosswithres);
    static_assert(std::is_same_v<decltype(tuns), decltype(turens)>);
}

TEST(DifferentialOperator2, derivativeCheck) {
    auto [x, y] = Init<2>();
    auto res = (x * y) * (x * cos(y));

    auto dcross = d(y) * d(x);
    // auto dres = d(res);
    // EXPECT_EQ(derivative_non_null(dres, dcross), false);

    auto dres2 = d<2>(res);
    EXPECT_EQ(derivative_non_null(dres2, dcross), true);

    // auto d2x = d<2>(x);
    // auto dres22 = d(res) * d(res);
    // EXPECT_EQ(derivative_non_null(dres22, d2x), false);
}

TEST(DifferentialOperator2, calcTreeOrderUsage) {
    auto [x, y] = Init<2>();
    auto res1 = x * y;
    auto res2 = x + y;
    CalcTree t(res1, res2);

    auto dx = d(x);
    auto dy = d(y);
    auto d1 = d(res1);
    auto d2 = d(res2);

    Tape2 t2(d1, d2, dx, dy);

    decltype(t2)::leaves a;
    std::cout << "leaves" << std::endl;
    std::cout << type_name2<decltype(a)>() << std::endl;

    decltype(t2)::roots a2;
    std::cout << "roots" << std::endl;
    std::cout << type_name2<decltype(a2)>() << std::endl;

    decltype(t)::ValuesTupleInverse b;
    std::cout << "calc tree" << std::endl;
    std::cout << type_name2<decltype(b)>() << std::endl;
}

TEST(DifferentialOperator2, dIDOrdered) {
    auto [x, y] = Init<2>();
    auto res = (x * y) * (x * cos(y));
    CalcTree t(res);
    decltype(t)::ValuesTupleInverse co;

    auto dx1 = pow<3>(d<2>(x));
    auto dy1 = pow<5>(d<4>(y));

    auto der3 = dx1 * dy1;
    auto der3_inv = dy1 * dx1;

    constexpr auto der3ordered = order_differential_operator(der3, co);
    constexpr auto der3_inv_ordered = order_differential_operator(der3_inv, co);

    static_assert(
        std::is_same_v<decltype(der3ordered), decltype(der3_inv_ordered)>);
}

TEST(DifferentialOperator2, dIDOrderedSameId) {
    auto [x, y] = Init<2>();
    auto res = (x * y) * (x * cos(y));
    CalcTree t(res);
    decltype(t)::ValuesTupleInverse co;

    auto dx1 = pow<3>(d<2>(x));
    auto dx2 = pow<5>(d<4>(x));
    auto dy1 = pow<3>(d<2>(y));
    auto dy2 = pow<5>(d<4>(y));

    auto der3 = dx1 * dy1 * dx2 * dy2 * dx1;
    auto der3_inv = dy2 * dx2 * dy1 * dx1 * dx1;

    std::cout << type_name2<decltype(der3)>() << std::endl;
    std::cout << type_name2<decltype(der3_inv)>() << std::endl;

    constexpr auto der3ordered = order_differential_operator(der3, co);
    constexpr auto der3_inv_ordered = order_differential_operator(der3_inv, co);

    std::cout << type_name2<decltype(der3ordered)>() << std::endl;
    std::cout << type_name2<decltype(der3_inv_ordered)>() << std::endl;

    static_assert(
        std::is_same_v<decltype(der3ordered), decltype(der3_inv_ordered)>);
}

TEST(DifferentialOperator2, lessThan) {
    auto [x, y] = Init<2>();
    auto res = (x * y) * (x * cos(y));
    CalcTree t(res);
    decltype(t)::ValuesTupleInverse co;

    auto const dx1 = pow<3>(d<2>(x));
    auto dx2 = pow<5>(d<4>(x));
    auto dx3 = pow<6>(d<4>(x));
    auto dy1 = pow<3>(d<2>(y));
    auto dy2 = pow<5>(d<4>(y));

    auto dxy1 = dx1 * dy1;
    auto dxy2 = dx1 * dy2;

    auto dxy3 = dy1 * dx1;
    auto dxy4 = dy1 * dx2;

    static_assert(!less_than(co, std::tuple<>{}, std::tuple<>{}));
    static_assert(!less_than(co, dx1, dx1));
    static_assert(less_than(co, std::tuple<>{}, dx1));
    static_assert(!less_than(co, dx1, std::tuple<>{}));
    static_assert(!less_than(co, dx1, dy1));
    static_assert(less_than(co, dx1, dx2));
    static_assert(less_than(co, dx2, dx3));
    static_assert(less_than(co, dy1, dy2));
    static_assert(less_than(co, dxy1, dxy2));
    static_assert(less_than(co, dxy3, dxy4));
}

TEST(DifferentialOperator2, merge) {
    auto [x1, x2, x3, x4] = Init<4>();

    auto res = (x1 * x2) * (x3 * cos(x4));
    CalcTree t(res);
    decltype(t)::ValuesTupleInverse co;

    auto dx1 = d(x1);
    auto dx2 = d(x2);
    auto dx3 = d(x3);
    auto dx4 = d(x4);

    auto r = std::tuple<decltype(dx1), decltype(dx2), decltype(dx3),
                        decltype(dx4)>{};
    auto t1 = std::tuple<decltype(dx1), decltype(dx2)>{};
    auto t2 = std::tuple<decltype(dx3), decltype(dx4)>{};

    auto r1 = merge_ordered(co, t1, t2);

    static_assert(std::is_same_v<decltype(r1), decltype(r)>);

    auto dxx1 = d(x1);
    auto dxx2 = d<2>(x1);
    auto dxx3 = d<3>(x1);
    auto dxx4 = d<4>(x1);

    auto rx = std::tuple<decltype(dxx4), decltype(dxx3), decltype(dxx2),
                         decltype(dxx1)>{};
    auto tx1 = std::tuple<decltype(dxx2), decltype(dxx1)>{};
    auto tx2 = std::tuple<decltype(dxx4), decltype(dxx3)>{};

    auto rx1 = merge_ordered(co, tx1, tx2);

    static_assert(std::is_same_v<decltype(rx1), decltype(rx)>);
}

} // namespace adhoc3
