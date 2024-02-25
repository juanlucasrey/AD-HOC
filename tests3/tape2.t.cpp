#include "type_name.hpp"
#include <differential_operator_util.hpp>
#include <init.hpp>
#include <tape2.hpp>

#include <calc_tree.hpp>

#include <gtest/gtest.h>

namespace adhoc3 {

TEST(Tape2, First) {
    auto [x, y] = Init<2>();
    auto res = (x * y) * (x * cos(y));

    auto dx = d(x);
    auto dy = d(y);
    auto dcross = d(y) * d(x);
    auto dres = d(res);
    auto dcrosswithres = d(res) * d(x);

    auto t = Tape2(dx, dy, dcross, dres, dcrosswithres);

    t.set(dx) = 1.;
    t.set(dy) = 2.;
    t.set(dcross) = 3.;
    t.set(dres) = 4.;
    t.set(dcrosswithres) = 5.;

    EXPECT_EQ(t.get(dx), 1.);
    EXPECT_EQ(t.get(dy), 2.);
    EXPECT_EQ(t.get(dcross), 3.);
    EXPECT_EQ(t.get(dres), 4.);
    EXPECT_EQ(t.get(dcrosswithres), 5.);
}

TEST(Tape2, TapeAndTree) {
    auto [x, y] = Init<2>();
    auto res = (x * y) * (x * cos(y));
}

} // namespace adhoc3
