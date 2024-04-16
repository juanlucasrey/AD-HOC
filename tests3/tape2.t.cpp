#include "type_name.hpp"
#include <differential_operator/differential_operator.hpp>
#include <differential_operator/select_root_derivatives.hpp>
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

// TEST(Tape2, TapeAndTree) {
//     auto [x, y] = Init<2>();
//     auto res = (x * y) * (x * cos(y));
//     auto res2 = cos(x) * y;
//     CalcTree ct(res, res2);

//     auto dx = d(x);
//     auto dy = d(y);
//     auto dcross = d(y) * d(x);
//     auto dres = d(res);
//     auto dcrosswithres = d(res) * d(x);
//     auto t = Tape2(dx, dy, dcross, dres, dcrosswithres);
//     t.set(dx) = 1.;
//     t.set(dy) = 2.;
//     t.set(dcross) = 3.;
//     t.set(dres) = 4.;
//     t.set(dcrosswithres) = 5.;

//     t.backpropagate(ct);
// }

TEST(Tape2, TapeAndTreeUnivariate) {
    auto [x] = Init<1>();
    auto res = exp(x);
    CalcTree ct(res);

    auto dx = d(x);
    auto dx2 = pow<2>(d(x));
    auto dx3 = pow<3>(d(x));
    auto dx4 = pow<4>(d(x));
    auto dres = d(res);
    auto dres2 = d<2>(res);
    auto dres3 = d<3>(res);
    auto dres4 = d<4>(res);
    auto t = Tape2(dx, dx2, dx3, dx4, dres, dres2, dres3, dres4);
    t.set(dres) = 1.;
    t.set(dres2) = 1.;
    t.set(dres3) = 1.;
    t.set(dres4) = 1.;
    t.backpropagate(ct);
}

// TEST(Tape2, TapeAndTree2) {
//     auto [x, y] = Init<2>();
//     auto res = (x * y) * (x * cos(y));
//     auto res2 = cos(x) * y;

//     auto t = Tape2(res, res2, d(x), d(y), d(y) * d(x), d(res), d(res2));
//     t.set(x) = 1.0;
//     t.set(y) = 1.0;
//     t.calcprimal();

//     t.set(d(res)) = 1.0;
//     t.set(d(res2)) = 1.0;
//     t.backpropagate();

//     auto dres1 = t.get(d(x));
//     auto dres2 = t.get(d(y));
//     auto dres3 = t.get(d(y) * d(x));
// }

} // namespace adhoc3
