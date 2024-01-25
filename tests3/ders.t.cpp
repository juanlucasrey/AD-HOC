#include <dependency.hpp>
#include <ders.hpp>
#include <init.hpp>
#include <strict_order.hpp>

// #include "type_name.hpp"

#include <gtest/gtest.h>

namespace adhoc3 {

TEST(Ders, MultiVarDefinition) {
    auto [S, K, T, R] = Init<4>();

    std::cout << detail::type_name<decltype(S)>() << std::endl;

    auto v1 = var(S);
    std::cout << detail::type_name<decltype(v1)>() << std::endl;

    auto v2 = var<2>(K);
    std::cout << detail::type_name<decltype(v2)>() << std::endl;

    auto v3 = var(v1, v2);
    std::cout << detail::type_name<decltype(v3)>() << std::endl;

    auto v3_bis = var(v2, v1);
    static_assert(std::is_same_v<decltype(v3), decltype(v3_bis)>);
    std::cout << detail::type_name<decltype(v3_bis)>() << std::endl;

    auto v4 = var(var<2>(S), var(K));
    std::cout << detail::type_name<decltype(v4)>() << std::endl;

    auto v5 = var(var<2>(S), var<3>(K), var(T));
    std::cout << detail::type_name<decltype(v5)>() << std::endl;

    // auto v6 = var(var<2>(S), var<3>(K), var(T), var<4>(K));
    auto v6 = var(var<2>(S), var<3>(K), var(T), var<4>(R));
    std::cout << detail::type_name<decltype(v6)>() << std::endl;

    // this is better
    // d<1, 2>(S, K);
    // d(S);
    // d(S, K);
    // Tape t<r1, r2>(d<2>(r1), r2, d(v1), d(v1, v2), d<2, 2>(v1, v2));
}

TEST(Ders, DerivativeNonNull) {

    auto [S, K, T, R] = Init<4>();

    auto v2 = var<2>(S);
    auto r2 = S * S;
    static_assert(der_non_null(r2, v2));

    auto v3 = var<3>(S);
    auto r3 = S * S;
    static_assert(!der_non_null(r3, v3));

    auto v4 = var(var<2>(S), var<2>(K));
    auto r4 = S * S * K * K;
    static_assert(der_non_null(r4, v4));

    auto v5 = var(var<15>(S), var<2>(K));
    auto r5 = exp(S) * K * K;
    static_assert(der_non_null(r5, v5));

    auto v6 = var(var<15>(S), var<3>(K));
    auto r6 = exp(S) * K * K;
    static_assert(!der_non_null(r6, v6));
}

// template <std::size_t Input1, std::size_t... Inputs>
// auto constexpr sum_no_overflow2() -> std::size_t {
//     if constexpr (sizeof...(Inputs) == 0) {
//         return Input1;
//     } else {
//         if (std::numeric_limits<std::size_t>::max() - Input1 <
//             sum_no_overflow2<Inputs...>()) {
//             return std::numeric_limits<std::size_t>::max();
//         }

//         if (std::numeric_limits<std::size_t>::max() -
//                 sum_no_overflow2<Inputs...>() <
//             Input1) {
//             return std::numeric_limits<std::size_t>::max();
//         }

//         return Input1 + sum_no_overflow2<Inputs...>();
//     }
// }

// TEST(Ders, DerivativeNonNull2) {
//     auto [x, y] = Init<2>();

//     auto prod = x * y;
//     auto prodf = cos(x) * exp(y);

//     der::m<der::p<der::d<decltype(prod), 1>, 2>,
//            der::p<der::d<decltype(prodf), 1>, 1>>
//         dnode{};

//     static_assert(order<decltype(prod), decltype(x)>::call() == 1);
//     static_assert(order<decltype(prod), decltype(y)>::call() == 1);
//     static_assert(order<decltype(prodf), decltype(x)>::call() ==
//                   std::numeric_limits<std::size_t>::max());
//     static_assert(order<decltype(prodf), decltype(y)>::call() ==
//                   std::numeric_limits<std::size_t>::max());

//     der::m<der::p<der::d<decltype(x), 1>, 3>> din1{};
//     static_assert(der_non_null(dnode, din1));

//     der::m<der::p<der::d<decltype(x), 1>, 2>, der::p<der::d<decltype(y), 1>,
//     1>>
//         din2{};
//     static_assert(der_non_null(dnode, din2));

//     der::m<der::p<der::d<decltype(x), 1>, 1>, der::p<der::d<decltype(y), 1>,
//     2>>
//         din3{};
//     static_assert(der_non_null(dnode, din3));

//     der::m<der::p<der::d<decltype(y), 1>, 3>> din4{};
//     static_assert(der_non_null(dnode, din4));
// }

// TEST(Ders, DerivativeNonNull3) {
//     auto [x, y] = Init<2>();

//     auto prod = x * y;

//     der::m<der::p<der::d<decltype(prod), 1>, 1>,
//            der::p<der::d<decltype(prod), 2>, 1>>
//         dnode{};

//     der::m<der::p<der::d<decltype(x), 1>, 3>> din1{};
//     static_assert(!der_non_null(dnode, din1));

//     der::m<der::p<der::d<decltype(x), 1>, 2>, der::p<der::d<decltype(y), 1>,
//     1>>
//         din2{};
//     static_assert(der_non_null(dnode, din2));

//     der::m<der::p<der::d<decltype(x), 1>, 1>, der::p<der::d<decltype(y), 1>,
//     2>>
//         din3{};
//     static_assert(der_non_null(dnode, din3));

//     der::m<der::p<der::d<decltype(y), 1>, 3>> din4{};
//     static_assert(!der_non_null(dnode, din4));
// }

} // namespace adhoc3