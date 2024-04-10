#include <differential_operator/expand.hpp>
#include <differential_operator/expand_tree.hpp>

#include <init.hpp>
#include <tape2.hpp>

#include "../type_name.hpp"

#include <gtest/gtest.h>

namespace adhoc3 {

template <bool... Ints> class bool_sequence {};

// template <bool... Ints> class bool_sequence::std::tuple<Ints...>;

TEST(ExpandTree, First) {
    auto [x1, x2] = Init<2>();
    auto res = exp(x1) * cos(x2) + x1;
    CalcTree t(res);
    decltype(t)::ValuesTupleInverse co;
    auto d1 = d(res);
    auto d2 = d<3>(res);
    auto outders = std::tuple<decltype(d2), decltype(d1)>{};

    auto do1 = d(x1);
    auto do2 = d<3>(x2);
    auto dersout = std::tuple<decltype(do1), decltype(do2)>{};

    auto result = expand_tree(co, dersout, outders);

    auto result2 = std::tuple<decltype(d2), decltype(d1),
                              decltype(d<3>(exp(x1) * cos(x2))),
                              decltype(d(exp(x1) * cos(x2))),
                              decltype(d(exp(x1))), decltype(d<3>(cos(x2)))>{};
    static_assert(std::is_same_v<decltype(result), decltype(result2)>);
}

} // namespace adhoc3
