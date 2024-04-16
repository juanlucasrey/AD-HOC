#include <differential_operator/expand_tree.hpp>
#include <init.hpp>
#include <tape2.hpp>
#include <tape_buffer_size.hpp>

#include <differential_operator/less_than.hpp>
#include <differential_operator/sort.hpp>

#include "call_price.hpp"

#include "type_name.hpp"

#include <cmath>

#include <gtest/gtest.h>

namespace adhoc3 {

TEST(TapeBufferSize, First) {
    auto [x1, x2] = Init<2>();
    auto res = exp(x1) * cos(x2) + x1;
    CalcTree t(res);
    decltype(t)::ValuesTupleInverse co;
    auto d1 = d(res);
    auto d2 = d<3>(res);
    auto dersin = std::tuple<decltype(d2), decltype(d1)>{};

    auto do1 = d(x1);
    auto do2 = d<3>(x2);
    auto dersout = std::tuple<decltype(do1), decltype(do2)>{};

    auto result = expand_tree(co, dersout, dersin);
    constexpr auto size = tape_buffer_size(co, result, dersin);
    static_assert(size == 2);
    // std::cout << size << std::endl;
}

// TEST(TapeBufferSize, BS) {
//     auto [S, K, v, T] = Init<4>();

//     auto res = call_price(S, K, v, T);

//     CalcTree t(res);
//     decltype(t)::ValuesTupleInverse co;

//     auto der1 = d(res);
//     auto der2 = d<2>(res);
//     auto dersin = std::tuple<decltype(der2), decltype(der1)>{};
//     auto dersout = std::make_tuple(d(S) * d(v), d(S), d(K), d(v), d(T));

//     auto dd1 = d(S) * d(K);
//     auto dd2 = pow<2>(d(S));

//     auto result = expand_tree(co, dersout, dersin);
//     constexpr auto size = tape_buffer_size(co, result, dersin);
//     static_assert(size == 9);
// }

// TEST(TapeBufferSize, Bug1) {
//     auto [x] = Init<1>();
//     auto res = x * x;
//     CalcTree t(res);
//     decltype(t)::ValuesTupleInverse co;
//     auto der1 = d(res);
//     // auto der2 = d<2>(res);
//     auto dersin = std::tuple<decltype(der1)>{};
//     auto dersout = std::make_tuple(d(x));
//     auto result = expand_tree(co, dersout, dersin);
//     std::cout << type_name2<decltype(result)>() << std::endl;
// }

// template <std::size_t Order, class Id, std::size_t... I>
// constexpr auto generate_operators_mul3(std::index_sequence<I...> /* i */) {
//     if constexpr ((Order / 2) * 2 == Order) {
//         return std::make_tuple(
//             std::tuple<der2::p<1, der2::d<Order, Id>>>{},
//             std::tuple<der2::p<1, der2::d<Order - I - 1, Id>>,
//                        der2::p<1, der2::d<I + 1, Id>>>{}...,
//             std::tuple<der2::p<2, der2::d<Order / 2, Id>>>{});
//     } else {
//         return std::make_tuple(
//             std::tuple<der2::p<1, der2::d<Order, Id>>>{},
//             std::tuple<der2::p<1, der2::d<Order - I - 1, Id>>,
//                        der2::p<1, der2::d<I + 1, Id>>>{}...);
//     }
// }

// TEST(TapeBufferSize, Bug2) {
//     auto [x1] = Init<1>();
//     constexpr std::size_t Order = 2;
//     constexpr std::size_t OrderDiv2 = (Order - 1) / 2;
//     constexpr auto seq = std::make_index_sequence<OrderDiv2>{};
//     constexpr auto restype = generate_operators_mul3<Order,
//     decltype(x1)>(seq); std::cout << type_name2<decltype(restype)>() <<
//     std::endl;
// }

} // namespace adhoc3
