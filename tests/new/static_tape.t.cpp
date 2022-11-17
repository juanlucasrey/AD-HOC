#include <new/adhoc.hpp>
#include <new/base.hpp>
#include <new/calc_order.hpp>
#include <new/init.hpp>
#include <new/tape_nodes.hpp>

#include "../type_name.hpp"

#include <gtest/gtest.h>

namespace adhoc2 {

TEST(StaticTape, TapeNodes) {

    auto [v0, v1, v2, v3] = Init<4>();
    auto res1 = v0 * v1 + v2 * v3;
    auto res2 = v0 + v1 * v2 + v3;

    auto leaf_and_root_tape = TapeRootsAndLeafs(res1, res2);
    leaf_and_root_tape.set(v0, 1.0);
    leaf_and_root_tape.set(v1, 2.0);
    leaf_and_root_tape.set(v2, 3.0);
    leaf_and_root_tape.set(v3, 4.0);

    auto derivatives = Tape(v0, v1, v2, res1, res2);

    auto calcs = calc_order_t(res1, res2);
    std::cout << type_name<decltype(calcs)>() << std::endl;
}

TEST(StaticTape, TapeNodes2) {

    auto [v0, v1, v2, v3] = Init<4>();
    auto res1 = v0 * v1 + exp(v2);
    auto res2 = cos(res1);
    auto res3 = log(res1);

    auto const v1_2 = v1;
    auto calc_c = v1_2 * v2;
    std::cout << type_name<decltype(calc_c)>() << std::endl;

    auto calcs = calc_order_t(res3, res2);
    std::cout << type_name<decltype(v0)>() << std::endl;
    std::cout << type_name<decltype(v0 * v1)>() << std::endl;
    auto const mul = v0 * v1;
    std::cout << type_name<decltype(mul)>() << std::endl;

    std::cout << type_name<decltype(calcs)>() << std::endl;
}

} // namespace adhoc2
