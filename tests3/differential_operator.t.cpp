#include <calc_order.hpp>
#include <differential_operator.hpp>
#include <init.hpp>
#include <tuple_utils.hpp>
#include <utils.hpp>

#include "type_name.hpp"

#include <gtest/gtest.h>

namespace adhoc3 {

TEST(DifferentialOperator, dID) {
    auto [x, y] = Init<2>();
    auto res = (x * y) * (x * cos(y));
    constexpr auto co = calc_order_t<true>(res);

    std::tuple<der2::p<der2::d<decltype(x), 2>, 3>,
               der2::p<der2::d<decltype(y), 4>, 5>>
        der3{};
    std::tuple<der2::p<der2::d<decltype(y), 4>, 5>,
               der2::p<der2::d<decltype(x), 2>, 3>>
        der3_inv{};

    constexpr auto der3ordered = order_differential_operator(der3, co);
    constexpr auto der3_inv_ordered = order_differential_operator(der3_inv, co);

    static_assert(
        std::is_same<decltype(der3ordered), decltype(der3_inv_ordered)>::value);
}

TEST(DifferentialOperator, DerivativeNonNull2) {
    auto [x, y] = Init<2>();

    auto prod = x * y;
    auto prodf = cos(x) * exp(y);

    std::tuple<der2::p<der2::d<decltype(prod), 1>, 2>,
               der2::p<der2::d<decltype(prodf), 1>, 1>>
        dnode{};

    std::cout << order<decltype(prod), decltype(x)>::call() << std::endl;
    std::cout << order<decltype(prod), decltype(y)>::call() << std::endl;
    std::cout << order<decltype(prodf), decltype(x)>::call() << std::endl;
    std::cout << order<decltype(prodf), decltype(y)>::call() << std::endl;

    std::tuple<der2::p<der2::d<decltype(x), 1>, 3>> din1{};
    std::cout << der_non_null(dnode, din1) << std::endl;

    std::tuple<der2::p<der2::d<decltype(x), 1>, 2>,
               der2::p<der2::d<decltype(y), 1>, 1>>
        din2{};
    std::cout << der_non_null(dnode, din2) << std::endl;

    std::tuple<der2::p<der2::d<decltype(x), 1>, 1>,
               der2::p<der2::d<decltype(y), 1>, 2>>
        din3{};
    std::cout << der_non_null(dnode, din3) << std::endl;

    std::tuple<der2::p<der2::d<decltype(y), 1>, 3>> din4{};
    std::cout << der_non_null(dnode, din4) << std::endl;
}

TEST(DifferentialOperator, DerivativeNonNull3) {
    auto [x, y] = Init<2>();

    auto prod = x * y;

    std::tuple<der2::p<der2::d<decltype(prod), 1>, 1>,
               der2::p<der2::d<decltype(prod), 2>, 1>>
        dnode{};

    auto din1 = diff_op<3>(x);
    std::cout << der_non_null(dnode, din1) << std::endl;

    std::tuple<der2::p<der2::d<decltype(x), 1>, 2>,
               der2::p<der2::d<decltype(y), 1>, 1>>
        din2{};

    auto din2_2 = diff_op(diff_op<3>(x), diff_op<3>(y));
    std::cout << type_name2<decltype(din2_2)>() << std::endl;

    std::cout << der_non_null(dnode, din2) << std::endl;

    std::tuple<der2::p<der2::d<decltype(x), 1>, 1>,
               der2::p<der2::d<decltype(y), 1>, 2>>
        din3{};
    std::cout << der_non_null(dnode, din3) << std::endl;

    auto din4 = diff_op<3>(y);
    std::cout << der_non_null(dnode, din4) << std::endl;
}

} // namespace adhoc3