#include "type_name.hpp"
#include <differential_operator/differential_operator.hpp>
#include <differential_operator/select_root_derivatives.hpp>
#include <init.hpp>
#include <tape2.hpp>

#include <calc_tree.hpp>

#include <gtest/gtest.h>

namespace adhoc3 {

TEST(Tape2, TapeAndTreeUnivariateLog) {
    double const val = 3.2;

    auto [x] = Init<1>();
    auto res = log(x);
    CalcTree ct(res);
    ct.set(x) = val;
    ct.evaluate();

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

    // from sympy import *
    // x = Symbol('x')
    // f = ln(x)
    // print(lambdify([x], f.diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x).diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x).diff(x).diff(x))(3.2))

    EXPECT_NEAR(t.get(dx), 0.3125, 1e-15);
    EXPECT_NEAR(t.get(dx2), -0.09765624999999999, 1e-15);
    EXPECT_NEAR(t.get(dx3), 0.061035156249999986, 1e-15);
    EXPECT_NEAR(t.get(dx4), -0.057220458984374986, 1e-15);
}

TEST(Tape2, TapeAndTreeUnivariateLogLog) {
    double const val = 3.2;

    auto [x] = Init<1>();
    auto res = log(log(x));
    CalcTree ct(res);
    ct.set(x) = val;
    ct.evaluate();

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

    // from sympy import *
    // x = Symbol('x')
    // f = ln(ln(x))
    // print(lambdify([x], f.diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x).diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x).diff(x).diff(x))(3.2))

    EXPECT_NEAR(t.get(dx), 0.26866679485199957, 1e-15);
    EXPECT_NEAR(t.get(dx2), -0.15614022004729622, 1e-15);
    EXPECT_NEAR(t.get(dx3), 0.15893019538473163, 1e-15);
    EXPECT_NEAR(t.get(dx4), -0.2307180116199667, 1e-15);
}

TEST(Tape2, TapeAndTreeUnivariateLogLogLog) {
    double const val = 3.2;

    auto [x] = Init<1>();
    auto res = log(log(log(x)));
    CalcTree ct(res);
    ct.set(x) = val;
    ct.evaluate();

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

    // from sympy import *
    // x = Symbol('x')
    // f = ln(ln(ln(x)))
    // print(lambdify([x], f.diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x).diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x).diff(x).diff(x))(3.2))

    EXPECT_NEAR(t.get(dx), 1.7776899535583102, 1e-15);
    EXPECT_NEAR(t.get(dx2), -4.193315943532229, 1e-15);
    EXPECT_NEAR(t.get(dx3), 17.797018671582318, 5e-14);
    EXPECT_NEAR(t.get(dx4), -111.30552385962602, 5e-14);
}

TEST(Tape2, TapeAndTreeUnivariateSumLog) {
    using constants::CD;
    using constants::encode;

    double const val = 3.2;

    auto [x] = Init<1>();
    auto res = log(x) + CD<encode(-0.5)>();
    CalcTree ct(res);
    ct.set(x) = val;
    ct.evaluate();

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

    // from sympy import *
    // x = Symbol('x')
    // f = ln(x) - 0.5
    // print(lambdify([x], f.diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x).diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x).diff(x).diff(x))(3.2))

    EXPECT_NEAR(t.get(dx), 0.3125, 1e-15);
    EXPECT_NEAR(t.get(dx2), -0.09765624999999999, 1e-15);
    EXPECT_NEAR(t.get(dx3), 0.061035156249999986, 1e-15);
    EXPECT_NEAR(t.get(dx4), -0.057220458984374986, 1e-15);
}

TEST(Tape2, TapeAndTreeUnivariateLogSum) {
    using constants::CD;
    using constants::encode;

    double const val = 3.2;

    auto [x] = Init<1>();
    auto res = log(x + CD<encode(0.3)>());
    CalcTree ct(res);
    ct.set(x) = val;
    ct.evaluate();

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

    // from sympy import *
    // x = Symbol('x')
    // f = ln(x + 0.3)
    // print(lambdify([x], f.diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x).diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x).diff(x).diff(x))(3.2))

    EXPECT_NEAR(t.get(dx), 0.2857142857142857, 1e-15);
    EXPECT_NEAR(t.get(dx2), -0.08163265306122448, 1e-15);
    EXPECT_NEAR(t.get(dx3), 0.04664723032069971, 1e-15);
    EXPECT_NEAR(t.get(dx4), -0.03998334027488547, 1e-15);
}

TEST(Tape2, TapeAndTreeUnivariateSumLogSum) {
    using constants::CD;
    using constants::encode;

    double const val = 3.2;

    auto [x] = Init<1>();
    auto res = log(x + CD<encode(0.3)>()) + CD<encode(0.5)>();
    CalcTree ct(res);
    ct.set(x) = val;
    ct.evaluate();

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

    // from sympy import *
    // x = Symbol('x')
    // f = ln(x + 0.3) + 0.5
    // print(lambdify([x], f.diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x).diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x).diff(x).diff(x))(3.2))

    EXPECT_NEAR(t.get(dx), 0.2857142857142857, 1e-15);
    EXPECT_NEAR(t.get(dx2), -0.08163265306122448, 1e-15);
    EXPECT_NEAR(t.get(dx3), 0.04664723032069971, 1e-15);
    EXPECT_NEAR(t.get(dx4), -0.03998334027488547, 1e-15);
}

TEST(Tape2, TapeAndTreeUnivariateSum) {
    double const valx = 3.2;
    double const valy = 0.44;

    auto [x, y] = Init<2>();
    auto res = log(x) + log(y);
    CalcTree ct(res);
    ct.set(x) = valx;
    ct.set(y) = valy;
    ct.evaluate();

    auto dx = d(x);
    auto dy = d(y);
    auto dx2 = pow<2>(d(x));
    auto dy2 = pow<2>(d(y));
    // auto dxy = d(x) * d(y);
    auto dres = d(res);
    auto dres2 = d<2>(res);
    auto t = Tape2(dx2, dx, dy2, dy, dres, dres2);
    t.set(dres) = 1.;
    t.set(dres2) = 1.;
    t.backpropagate(ct);

    // from sympy import *
    // x = Symbol('x')
    // y = Symbol('y')
    // f = log(x) + log(y)
    // print(lambdify([x, y], f.diff(x))(3.2, 0.44))
    // print(lambdify([x, y], f.diff(x).diff(x))(3.2, 0.44))
    // print(lambdify([x, y], f.diff(y))(3.2, 0.44))
    // print(lambdify([x, y], f.diff(y).diff(y))(3.2, 0.44))

    EXPECT_NEAR(t.get(dx), 0.3125, 1e-15);
    EXPECT_NEAR(t.get(dx2), -0.09765624999999999, 1e-15);
    EXPECT_NEAR(t.get(dy), 2.272727272727273, 1e-15);
    EXPECT_NEAR(t.get(dy2), -5.1652892561983474, 1e-15);
}

TEST(Tape2, TapeAndTreeUnivariateSameSum) {
    using constants::CD;
    using constants::encode;

    double const val = 3.2;

    auto [x] = Init<1>();
    auto res = log(x) + log(x);
    CalcTree ct(res);
    ct.set(x) = val;
    ct.evaluate();

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

    // from sympy import *
    // x = Symbol('x')
    // f = ln(x + 0.3)
    // print(lambdify([x], f.diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x).diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x).diff(x).diff(x))(3.2))

    EXPECT_NEAR(t.get(dx), 0.625, 1e-15);
    EXPECT_NEAR(t.get(dx2), -0.19531249999999997, 1e-15);
    EXPECT_NEAR(t.get(dx3), 0.12207031249999997, 1e-15);
    EXPECT_NEAR(t.get(dx4), -0.11444091796874997, 1e-15);
    // std::cout << t.get(dx) << std::endl;
    // std::cout << t.get(dx2) << std::endl;
    // std::cout << t.get(dx3) << std::endl;
    // std::cout << t.get(dx4) << std::endl;
}

TEST(Tape2, TapeAndTreeUnivariateLargeSum) {
    double const valx = 3.2;
    double const valy = 0.44;

    auto [x, y] = Init<2>();
    auto res = log(log(log(x) + log(y)));
    CalcTree ct(res);
    ct.set(x) = valx;
    ct.set(y) = valy;
    ct.evaluate();

    auto dx = d(x);
    auto dy = d(y);
    auto dx2 = pow<2>(d(x));
    auto dy2 = pow<2>(d(y));
    auto dxy = d(x) * d(y);

    auto dx3 = pow<3>(d(x));
    auto dx2y1 = pow<2>(d(x)) * d(y);
    auto dx1y2 = d(x) * pow<2>(d(y));
    auto dy3 = pow<3>(d(y));

    auto dres = d(res);
    auto dres2 = d<2>(res);
    auto dres3 = d<3>(res);

    auto t = Tape2(dx3, dx2y1, dx1y2, dy3, dx2, dx, dy2, dxy, dy, dres, dres2,
                   dres3);
    // auto t = Tape2(dx, dy, dres);
    t.set(dres) = 1.;
    t.set(dres2) = 1.;
    t.set(dres3) = 1.;
    t.backpropagate(ct);

    // from sympy import *
    // x = Symbol('x')
    // y = Symbol('y')
    // f = log(log(log(x) + log(y)))
    // print(lambdify([x, y], f.diff(x))(3.2, 0.44))
    // print(lambdify([x, y], f.diff(y))(3.2, 0.44))
    // print(lambdify([x, y], f.diff(x).diff(x))(3.2, 0.44))
    // print(lambdify([x, y], f.diff(y).diff(y))(3.2, 0.44))
    // print(lambdify([x, y], f.diff(y).diff(x))(3.2, 0.44))
    // print(lambdify([x, y], f.diff(x).diff(x).diff(x))(3.2, 0.44))
    // print(lambdify([x, y], f.diff(x).diff(x).diff(y))(3.2, 0.44))
    // print(lambdify([x, y], f.diff(x).diff(y).diff(y))(3.2, 0.44))
    // print(lambdify([x, y], f.diff(y).diff(y).diff(y))(3.2, 0.44))

    EXPECT_NEAR(t.get_d(dx), -0.8515928225237093, 1e-15);
    EXPECT_NEAR(t.get_d(dy), -6.193402345626976, 1e-15);
    EXPECT_NEAR(t.get_d(dx2), 0.3186619513860184, 1e-15);
    EXPECT_NEAR(t.get_d(dy2), 16.85484701545883, 1e-14);
    EXPECT_NEAR(t.get_d(dxy), 0.3821032316171591, 1e-15);
    EXPECT_NEAR(t.get_d(dx3), -0.8843909678854258, 1e-15);
    EXPECT_NEAR(t.get_d(dx2y1), -4.983470896503015, 5e-14);
    EXPECT_NEAR(t.get_d(dx1y2), -36.24342470184011, 5e-14);
    EXPECT_NEAR(t.get_d(dy3), -340.2014842654683, 5e-13);

    // std::cout << t.get_d(dx) << std::endl;
    // std::cout << t.get_d(dy) << std::endl;
    // std::cout << t.get_d(dx2) << std::endl;
    // std::cout << t.get_d(dy2) << std::endl;
    // std::cout << t.get_d(dxy) << std::endl;
    // std::cout << t.get_d(dx3) << std::endl;
    // std::cout << t.get_d(dx2y1) << std::endl;
    // std::cout << t.get_d(dx1y2) << std::endl;
    // std::cout << t.get_d(dy3) << std::endl;
}

TEST(Tape2, TapeAndTreeUnivariateLargeSumSame) {
    double const valx = 3.2;

    auto [x] = Init<1>();
    auto res = log(log(log(x) + log(x)));
    CalcTree ct(res);
    ct.set(x) = valx;
    ct.evaluate();

    auto dx = d(x);
    auto dx2 = pow<2>(d(x));
    auto dx3 = pow<3>(d(x));
    auto dx4 = pow<4>(d(x));

    auto dres = d(res);
    auto dres2 = d<2>(res);
    auto dres3 = d<3>(res);
    auto dres4 = d<4>(res);

    auto t = Tape2(dx4, dx3, dx2, dx, dres, dres2, dres3, dres4);
    t.set(dres) = 1.;
    t.set(dres2) = 1.;
    t.set(dres3) = 1.;
    t.set(dres4) = 1.;
    t.backpropagate(ct);

    // from sympy import *
    // x = Symbol('x')
    // f = log(log(log(x) + log(x)))
    // print(lambdify([x], f.diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x).diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x).diff(x).diff(x))(3.2))

    EXPECT_NEAR(t.get_d(dx), 0.318220121691175, 1e-15);
    EXPECT_NEAR(t.get_d(dx2), -0.28620301402981957, 1e-15);
    EXPECT_NEAR(t.get_d(dx3), 0.4292459478298339, 1e-15);
    EXPECT_NEAR(t.get_d(dx4), -0.9017492480110629, 1e-15);
}

TEST(Tape2, TapeAndTreeUnivariateSameSub) {
    using constants::CD;
    using constants::encode;

    double const val = 3.2;

    auto [x] = Init<1>();
    auto res = log(x) - log(x);
    CalcTree ct(res);
    ct.set(x) = val;
    ct.evaluate();

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

    // from sympy import *
    // x = Symbol('x')
    // f = log(x) - log(x)
    // print(lambdify([x], f.diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x).diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x).diff(x).diff(x))(3.2))

    EXPECT_NEAR(t.get(dx), 0., 1e-15);
    EXPECT_NEAR(t.get(dx2), 0., 1e-15);
    EXPECT_NEAR(t.get(dx3), 0., 1e-15);
    EXPECT_NEAR(t.get(dx4), 0., 1e-15);
}

TEST(Tape2, TapeAndTreeUnivariateSubLogRightConst) {
    using constants::CD;
    using constants::encode;

    double const val = 3.2;

    auto [x] = Init<1>();
    auto res = log(x) - CD<encode(0.5)>();
    CalcTree ct(res);
    ct.set(x) = val;
    ct.evaluate();

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

    // from sympy import *
    // x = Symbol('x')
    // f = ln(x) - 0.5
    // print(lambdify([x], f.diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x).diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x).diff(x).diff(x))(3.2))

    EXPECT_NEAR(t.get(dx), 0.3125, 1e-15);
    EXPECT_NEAR(t.get(dx2), -0.09765624999999999, 1e-15);
    EXPECT_NEAR(t.get(dx3), 0.061035156249999986, 1e-15);
    EXPECT_NEAR(t.get(dx4), -0.057220458984374986, 1e-15);
    // std::cout << t.get(dx) << std::endl;
    // std::cout << t.get(dx2) << std::endl;
    // std::cout << t.get(dx3) << std::endl;
    // std::cout << t.get(dx4) << std::endl;
}

TEST(Tape2, TapeAndTreeUnivariateSubLogLeftConst) {
    using constants::CD;
    using constants::encode;

    double const val = 3.2;

    auto [x] = Init<1>();
    auto res = CD<encode(0.5)>() - log(x);
    CalcTree ct(res);
    ct.set(x) = val;
    ct.evaluate();

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

    // from sympy import *
    // x = Symbol('x')
    // f = ln(x) - 0.5
    // print(lambdify([x], f.diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x).diff(x))(3.2))
    // print(lambdify([x], f.diff(x).diff(x).diff(x).diff(x))(3.2))

    EXPECT_NEAR(t.get(dx), -0.3125, 1e-15);
    EXPECT_NEAR(t.get(dx2), 0.09765624999999999, 1e-15);
    EXPECT_NEAR(t.get(dx3), -0.061035156249999986, 1e-15);
    EXPECT_NEAR(t.get(dx4), 0.057220458984374986, 1e-15);
    // std::cout << t.get(dx) << std::endl;
    // std::cout << t.get(dx2) << std::endl;
    // std::cout << t.get(dx3) << std::endl;
    // std::cout << t.get(dx4) << std::endl;
}

TEST(Tape2, TapeAndTreeUnivariateLargeSub) {
    double const valx = 3.2;
    double const valy = 0.44;

    auto [x, y] = Init<2>();
    auto res = log(log(log(x) - log(y)));
    CalcTree ct(res);
    ct.set(x) = valx;
    ct.set(y) = valy;
    ct.evaluate();

    auto dx = d(x);
    auto dy = d(y);
    auto dx2 = pow<2>(d(x));
    auto dy2 = pow<2>(d(y));
    auto dxy = d(x) * d(y);

    auto dx3 = pow<3>(d(x));
    auto dx2y1 = pow<2>(d(x)) * d(y);
    auto dx1y2 = d(x) * pow<2>(d(y));
    auto dy3 = pow<3>(d(y));

    auto dres = d(res);
    auto dres2 = d<2>(res);
    auto dres3 = d<3>(res);

    auto t = Tape2(dx3, dx2y1, dx1y2, dy3, dx2, dx, dy2, dxy, dy, dres, dres2,
                   dres3);
    t.set(dres) = 1.;
    t.set(dres2) = 1.;
    t.set(dres3) = 1.;
    t.backpropagate(ct);

    // from sympy import *
    // x = Symbol('x')
    // y = Symbol('y')
    // f = log(log(log(x) - log(y)))
    // print(lambdify([x, y], f.diff(x))(3.2, 0.44))
    // print(lambdify([x, y], f.diff(y))(3.2, 0.44))
    // print(lambdify([x, y], f.diff(x).diff(x))(3.2, 0.44))
    // print(lambdify([x, y], f.diff(y).diff(y))(3.2, 0.44))
    // print(lambdify([x, y], f.diff(y).diff(x))(3.2, 0.44))
    // print(lambdify([x, y], f.diff(x).diff(x).diff(x))(3.2, 0.44))
    // print(lambdify([x, y], f.diff(x).diff(x).diff(y))(3.2, 0.44))
    // print(lambdify([x, y], f.diff(x).diff(y).diff(y))(3.2, 0.44))
    // print(lambdify([x, y], f.diff(y).diff(y).diff(y))(3.2, 0.44))

    EXPECT_NEAR(t.get_d(dx), 0.22986568884243355, 1e-15);
    EXPECT_NEAR(t.get_d(dy), -1.6717504643086079, 1e-15);
    EXPECT_NEAR(t.get_d(dx2), -0.1608750287834944, 1e-15);
    EXPECT_NEAR(t.get_d(dy2), -0.9102266846663642, 1e-14);
    EXPECT_NEAR(t.get_d(dxy), 0.647578189238065, 1e-15);
    EXPECT_NEAR(t.get_d(dx3), 0.18903408490935755, 1e-15);
    EXPECT_NEAR(t.get_d(dx2y1), -0.6435432139612622, 5e-14);
    EXPECT_NEAR(t.get_d(dx1y2), 1.7367770595452485, 5e-14);
    EXPECT_NEAR(t.get_d(dy3), -8.493711866391067, 5e-13);

    // std::cout << t.get_d(dx) << std::endl;
    // std::cout << t.get_d(dy) << std::endl;
    // std::cout << t.get_d(dx2) << std::endl;
    // std::cout << t.get_d(dy2) << std::endl;
    // std::cout << t.get_d(dxy) << std::endl;
    // std::cout << t.get_d(dx3) << std::endl;
    // std::cout << t.get_d(dx2y1) << std::endl;
    // std::cout << t.get_d(dx1y2) << std::endl;
    // std::cout << t.get_d(dy3) << std::endl;
}

TEST(Tape2, TapeAndTreeUnivariateLargeMult) {
    double const valx = 3.2;
    double const valy = 0.44;

    auto [x, y] = Init<2>();
    auto res = log(log(x * y));
    CalcTree ct(res);
    ct.set(x) = valx;
    ct.set(y) = valy;
    ct.evaluate();

    auto dx = d(x);
    auto dy = d(y);
    auto dx2 = pow<2>(d(x));
    auto dy2 = pow<2>(d(y));
    auto dxy = d(x) * d(y);

    auto dx3 = pow<3>(d(x));
    auto dx2y1 = pow<2>(d(x)) * d(y);
    auto dx1y2 = d(x) * pow<2>(d(y));
    auto dy3 = pow<3>(d(y));

    auto dres = d(res);
    auto dres2 = d<2>(res);
    auto dres3 = d<3>(res);

    auto t = Tape2(dx3, dx2y1, dx1y2, dy3, dx2, dx, dy2, dxy, dy, dres, dres2,
                   dres3);
    t.set(dres) = 1.;
    t.set(dres2) = 1.;
    t.set(dres3) = 1.;
    t.backpropagate(ct);

    // from sympy import *
    // x = Symbol('x')
    // y = Symbol('y')
    // f = log(log(x * y))
    // print(lambdify([x, y], f.diff(x))(3.2, 0.44))
    // print(lambdify([x, y], f.diff(y))(3.2, 0.44))
    // print(lambdify([x, y], f.diff(x).diff(x))(3.2, 0.44))
    // print(lambdify([x, y], f.diff(y).diff(y))(3.2, 0.44))
    // print(lambdify([x, y], f.diff(y).diff(x))(3.2, 0.44))
    // print(lambdify([x, y], f.diff(x).diff(x).diff(x))(3.2, 0.44))
    // print(lambdify([x, y], f.diff(x).diff(x).diff(y))(3.2, 0.44))
    // print(lambdify([x, y], f.diff(x).diff(y).diff(y))(3.2, 0.44))
    // print(lambdify([x, y], f.diff(y).diff(y).diff(y))(3.2, 0.44))

    EXPECT_NEAR(t.get_d(dx), 0.9132880282109275, 1e-15);
    EXPECT_NEAR(t.get_d(dy), 6.642094750624928, 1e-15);
    EXPECT_NEAR(t.get_d(dx2), -1.1194975312893185, 1e-15);
    EXPECT_NEAR(t.get_d(dy2), -59.213092564063146, 1e-14);
    EXPECT_NEAR(t.get_d(dxy), -6.066145617988392, 1e-15);
    EXPECT_NEAR(t.get_d(dx3), 2.483878648409331, 1e-15);
    EXPECT_NEAR(t.get_d(dx2y1), 12.975946846207323, 5e-14);
    EXPECT_NEAR(t.get_d(dx1y2), 94.37052251787145, 5e-14);
    EXPECT_NEAR(t.get_d(dy3), 955.4814936029886, 5e-13);

    // std::cout.precision(std::numeric_limits<double>::max_digits10);
    // std::cout << t.get_d(dx) << std::endl;
    // std::cout << t.get_d(dy) << std::endl;
    // std::cout << t.get_d(dx2) << std::endl;
    // std::cout << t.get_d(dy2) << std::endl;
    // std::cout << t.get_d(dxy) << std::endl;
    // std::cout << t.get_d(dx3) << std::endl;
    // std::cout << t.get_d(dx2y1) << std::endl;
    // std::cout << t.get_d(dx1y2) << std::endl;
    // std::cout << t.get_d(dy3) << std::endl;
}

} // namespace adhoc3
