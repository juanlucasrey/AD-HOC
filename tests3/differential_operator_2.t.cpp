#include <differential_operator_2.hpp>
#include <differential_operator_util.hpp>
#include <init.hpp>

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
    static_assert(is_root_derivative(dx));
    auto dy = d(y);
    static_assert(is_root_derivative(dy));
    auto dcross = d(y) * d(x);
    static_assert(is_root_derivative(dcross));
    auto dres = d(res);
    static_assert(!is_root_derivative(dres));
    auto dcrosswithres = d(res) * d(x);
    static_assert(!is_root_derivative(dcrosswithres));
}

} // namespace adhoc3