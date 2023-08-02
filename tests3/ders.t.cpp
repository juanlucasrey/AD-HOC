#include <ders.hpp>
#include <init.hpp>

#include "type_name.hpp"

#include <gtest/gtest.h>

namespace adhoc3 {

TEST(Ders, First) {
    auto [S, K, T, R] = Init<4>();

    auto v1 = var(S);
    std::cout << type_name<decltype(v1)>() << std::endl;

    auto v2 = var<2>(K);
    std::cout << type_name<decltype(v2)>() << std::endl;

    auto v3 = var(v1, v2);
    std::cout << type_name<decltype(v3)>() << std::endl;

    auto v4 = var(var<2>(S), var(K));
    std::cout << type_name<decltype(v4)>() << std::endl;

    auto v5 = var(var<2>(S), var<3>(K), var(T));
    std::cout << type_name<decltype(v5)>() << std::endl;

    // auto v8 = var2(var2<2>(S), var2<3>(K), var2(T), var2<4>(K));
    auto v6 = var(var<2>(S), var<3>(K), var(T), var<4>(R));
    std::cout << type_name<decltype(v6)>() << std::endl;
}

} // namespace adhoc3