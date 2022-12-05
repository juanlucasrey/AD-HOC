#include <adhoc.hpp>
// #include <base.hpp>
// #include <calc_order.hpp>
// #include <init.hpp>
#include <tape.hpp>

// #include "type_name.hpp"

#include <gtest/gtest.h>

namespace adhoc {

TEST(Tape, TapeInitial) {
    double_t<0> val1;
    double_t<1> val2;
    auto res = val1 * val2;

    auto temp =
        detail::Tape2<std::tuple<decltype(val1), decltype(val2), decltype(res)>,
                      std::tuple<decltype(val1), decltype(val2), decltype(res)>,
                      2>{};

    temp.val(val1) = 1;
    temp.der(val1) = 1;
    temp.val(val2) = 1;
}

TEST(Tape, TapeInitial2) {
    double_t<0> val1;
    double_t<1> val2;
    auto res = val1 * val2;

    auto tape = Tape3(val1, res);
    tape.der_check(val1) = 1.0;
    // auto temp =
    //     Tape2<std::tuple<decltype(val1), decltype(val2), decltype(res)>,
    //           std::tuple<decltype(val1), decltype(val2), decltype(res)>,
    //           2>{};

    // temp.val(val1) = 1;
    // temp.der(val1) = 1;
    // temp.val(val2) = 1;
}

// TEST(StaticTape, TapeNodes2) {

//     auto [v0, v1, v2, v3] = Init<4>();
//     auto res1 = v0 * v1 + exp(v2);
//     auto res2 = cos(res1);
//     auto res3 = log(res1);

//     auto const v1_2 = v1;
//     auto calc_c = v1_2 * v2;
//     std::cout << type_name<decltype(calc_c)>() << std::endl;

//     auto calcs = calc_order_t(res3, res2);
//     std::cout << type_name<decltype(v0)>() << std::endl;
//     std::cout << type_name<decltype(v0 * v1)>() << std::endl;
//     auto const mul = v0 * v1;
//     std::cout << type_name<decltype(mul)>() << std::endl;

//     std::cout << type_name<decltype(calcs)>() << std::endl;
// }

} // namespace adhoc
