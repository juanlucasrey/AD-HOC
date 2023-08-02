#include <ders.hpp>
#include <init.hpp>
#include <strict_order.hpp>

// #include "type_name.hpp"

#include <gtest/gtest.h>

namespace adhoc3 {

TEST(Ders, First) {
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

} // namespace adhoc3