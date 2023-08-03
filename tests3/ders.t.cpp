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

} // namespace adhoc3