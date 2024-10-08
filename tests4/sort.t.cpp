#include <calc_tree.hpp>
#include <sort.hpp>
#include <utils/bivariate.hpp>

#include "type_name.hpp"
#include <gtest/gtest.h>

#include <random>
#include <tuple>

namespace adhoc4 {

TEST(BackPropagator, Sort) {
    ADHOC(x);
    ADHOC(y);
    ADHOC(z);
    auto res = x + (y * z);
    auto temp = (z * y);
    CalcTree ct(res);
    using PrimalNodes = decltype(ct)::ValuesTupleInverse;
    constexpr auto primal_nodes_inverted = PrimalNodes{};

    constexpr auto new_bi = detail::sort_bivariate(temp, primal_nodes_inverted);
    std::cout << type_name2<decltype(new_bi)>() << std::endl;

    auto d1 = d(res);
    auto d2 = d(res) * d(x);
    auto d3 = d(res) * d(y);
    auto d4 = d(res) * d(z);
    auto d5 = d<2>(res);
    auto d6 = d<2>(res) * d<2>(x);

    constexpr auto newtuple = sort_differential_operators(
        std::make_tuple(d4, d3, d2, d1, d5, d6), primal_nodes_inverted);
    std::cout << type_name2<decltype(newtuple)>() << std::endl;
}

} // namespace adhoc4
