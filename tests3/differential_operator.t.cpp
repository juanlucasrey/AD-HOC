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

    std::tuple<der::p<der::d<decltype(x), 2>, 3>,
               der::p<der::d<decltype(y), 4>, 5>>
        der3{};
    std::tuple<der::p<der::d<decltype(y), 4>, 5>,
               der::p<der::d<decltype(x), 2>, 3>>
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

    std::tuple<der::p<der::d<decltype(prod), 1>, 2>,
               der::p<der::d<decltype(prodf), 1>, 1>>
        dnode{};

    static_assert(order<decltype(prod), decltype(x)>::call() == 1);
    static_assert(order<decltype(prod), decltype(y)>::call() == 1);
    static_assert(order<decltype(prodf), decltype(x)>::call() ==
                  std::numeric_limits<std::size_t>::max());
    static_assert(order<decltype(prodf), decltype(y)>::call() ==
                  std::numeric_limits<std::size_t>::max());

    std::tuple<der::p<der::d<decltype(x), 1>, 3>> din1{};
    static_assert(der_non_null(dnode, din1));

    std::tuple<der::p<der::d<decltype(x), 1>, 2>,
               der::p<der::d<decltype(y), 1>, 1>>
        din2{};
    static_assert(der_non_null(dnode, din2));

    std::tuple<der::p<der::d<decltype(x), 1>, 1>,
               der::p<der::d<decltype(y), 1>, 2>>
        din3{};
    static_assert(der_non_null(dnode, din3));

    std::tuple<der::p<der::d<decltype(y), 1>, 3>> din4{};
    static_assert(der_non_null(dnode, din4));
}

TEST(DifferentialOperator, DerivativeNonNull3) {
    auto [x, y, z] = Init<3>();

    auto prod = x * y;

    std::tuple<der::p<der::d<decltype(prod), 1>, 1>,
               der::p<der::d<decltype(prod), 2>, 1>>
        dnode{};

    auto din1 = d<3>(x);
    static_assert(!der_non_null(dnode, din1));

    std::tuple<der::p<der::d<decltype(x), 1>, 2>,
               der::p<der::d<decltype(y), 1>, 1>>
        din2{};

    auto din2_2 = d<3>(x) * d<3>(y) * d(z);
    // std::cout << type_name2<decltype(din2_2)>() << std::endl;

    static_assert(der_non_null(dnode, din2));

    std::tuple<der::p<der::d<decltype(x), 1>, 1>,
               der::p<der::d<decltype(y), 1>, 2>>
        din3{};
    static_assert(der_non_null(dnode, din3));

    auto din4 = d<3>(y);
    static_assert(!der_non_null(dnode, din4));

    auto din1_1 = d(x);
    static_assert(der_non_null(dnode, din1_1));
}

} // namespace adhoc3