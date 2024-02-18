#include <differential_operator_2.hpp>
#include <init.hpp>

#include "type_name.hpp"

#include <gtest/gtest.h>

namespace adhoc3 {

TEST(DifferentialOperator2, dID) {

    auto [x, y] = Init<2>();
    auto res = (x * y) * (x * cos(y));

    auto d1 = d<2>(x);

    std::cout << type_name2<decltype(d1)>() << std::endl;

    auto d2 = pow<2>(d1);

    std::cout << type_name2<decltype(d2)>() << std::endl;
}

TEST(DifferentialOperator2, dID2) {

    auto [x, y] = Init<2>();
    auto res = (x * y) * (x * cos(y));

    auto d1 = d<2>(x) * d<2>(x);

    std::cout << type_name2<decltype(d1)>() << std::endl;

    auto d2 = d<2>(x) * d(x);

    std::cout << type_name2<decltype(d2)>() << std::endl;

    auto d3 = d1 * d2;

    std::cout << type_name2<decltype(d3)>() << std::endl;
}

} // namespace adhoc3