#include <differential_operator_2.hpp>
#include <differential_operator_util.hpp>
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
    std::cout << type_name2<decltype(a)>() << std::endl;

    decltype(t)::ValuesTupleInverse b;
    std::cout << type_name2<decltype(b)>() << std::endl;
}

} // namespace adhoc3