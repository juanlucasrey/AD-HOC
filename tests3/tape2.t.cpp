#include <calc_tree.hpp>
#include <differential_operator/differential_operator.hpp>
#include <differential_operator/select_root_derivatives.hpp>
#include <init.hpp>
#include <tape2.hpp>

#include "call_price.hpp"
#include "type_name.hpp"

#include <gtest/gtest.h>

#include <chrono>
#include <random>

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

TEST(Tape2, TapeAndTreeUnivariateLargeMult2) {
    double const valx = 0.5;
    double const valy = 0.44;

    auto [x, y] = Init<2>();
    auto res = log(log(cos(x) * sin(y)));
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
    // f = log(log(cos(x) * sin(y)))
    // valx = 0.5
    // valy = 0.44

    // print(lambdify([x, y], f.diff(x))(valx, valy))
    // print(lambdify([x, y], f.diff(y))(valx, valy))
    // print(lambdify([x, y], f.diff(x).diff(x))(valx, valy))
    // print(lambdify([x, y], f.diff(y).diff(y))(valx, valy))
    // print(lambdify([x, y], f.diff(y).diff(x))(valx, valy))
    // print(lambdify([x, y], f.diff(x).diff(x).diff(x))(valx, valy))
    // print(lambdify([x, y], f.diff(x).diff(x).diff(y))(valx, valy))
    // print(lambdify([x, y], f.diff(x).diff(y).diff(y))(valx, valy))
    // print(lambdify([x, y], f.diff(y).diff(y).diff(y))(valx, valy))

    EXPECT_NEAR(t.get_d(dx), 0.5551616009854359, 1e-15);
    EXPECT_NEAR(t.get_d(dy), -2.158578010960933, 1e-15);
    EXPECT_NEAR(t.get_d(dx2), 1.0112982558453643, 1e-15);
    EXPECT_NEAR(t.get_d(dy2), 0.9418621704413335, 1e-14);
    EXPECT_NEAR(t.get_d(dxy), 1.1983596244170294, 1e-15);
    EXPECT_NEAR(t.get_d(dx3), -0.4137099502892953, 1e-15);
    EXPECT_NEAR(t.get_d(dx2y1), 1.517682929943283, 5e-14);
    EXPECT_NEAR(t.get_d(dx1y2), 2.0638670240401735, 5e-14);
    EXPECT_NEAR(t.get_d(dy3), -7.638836895792224, 5e-13);

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

TEST(Tape2, TapeAndTreeUnivariateLargeMultConstRight) {
    using constants::CD;
    using constants::encode;
    double const valx = 0.44;

    auto [x] = Init<1>();
    auto res = log(cos(x) * CD<encode(0.5)>());
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
    // f = log(cos(x) * 0.5)
    // valx = 0.44
    // print(lambdify([x], f.diff(x))(valx))
    // print(lambdify([x], f.diff(x).diff(x))(valx))
    // print(lambdify([x], f.diff(x).diff(x).diff(x))(valx))
    // print(lambdify([x], f.diff(x).diff(x).diff(x).diff(x))(valx))

    EXPECT_NEAR(t.get_d(dx), -0.47078052727762176, 1e-15);
    EXPECT_NEAR(t.get_d(dx2), -1.2216343048637954, 1e-15);
    EXPECT_NEAR(t.get_d(dx3), -1.1502432843684172, 1e-15);
    EXPECT_NEAR(t.get_d(dx4), -4.067805029465111, 1e-15);

    // std::cout.precision(std::numeric_limits<double>::max_digits10);
    // std::cout << t.get_d(dx) << std::endl;
    // std::cout << t.get_d(dx2) << std::endl;
    // std::cout << t.get_d(dx3) << std::endl;
    // std::cout << t.get_d(dx4) << std::endl;
}

TEST(Tape2, TapeAndTreeUnivariateLargeMultConstLeft) {
    using constants::CD;
    using constants::encode;
    double const valx = 0.44;

    auto [x] = Init<1>();
    auto res = log(CD<encode(0.5)>() * cos(x));
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
    // f = log(cos(x) * 0.5)
    // valx = 0.44
    // print(lambdify([x], f.diff(x))(valx))
    // print(lambdify([x], f.diff(x).diff(x))(valx))
    // print(lambdify([x], f.diff(x).diff(x).diff(x))(valx))
    // print(lambdify([x], f.diff(x).diff(x).diff(x).diff(x))(valx))

    EXPECT_NEAR(t.get_d(dx), -0.47078052727762176, 1e-15);
    EXPECT_NEAR(t.get_d(dx2), -1.2216343048637954, 1e-15);
    EXPECT_NEAR(t.get_d(dx3), -1.1502432843684172, 1e-15);
    EXPECT_NEAR(t.get_d(dx4), -4.067805029465111, 1e-15);

    // std::cout.precision(std::numeric_limits<double>::max_digits10);
    // std::cout << t.get_d(dx) << std::endl;
    // std::cout << t.get_d(dx2) << std::endl;
    // std::cout << t.get_d(dx3) << std::endl;
    // std::cout << t.get_d(dx4) << std::endl;
}

TEST(Tape2, TapeAndTreeUnivariateLargeMultSame) {
    using constants::CD;
    using constants::encode;
    double const valx = 0.44;

    auto [x] = Init<1>();
    auto res = log(cos(x) * cos(x));
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
    // f = log(cos(x) * cos(x))
    // valx = 0.44

    // print(lambdify([x], f.diff(x))(valx))
    // print(lambdify([x], f.diff(x).diff(x))(valx))
    // print(lambdify([x], f.diff(x).diff(x).diff(x))(valx))
    // print(lambdify([x], f.diff(x).diff(x).diff(x).diff(x))(valx))

    EXPECT_NEAR(t.get_d(dx), -0.9415610545552435, 1e-15);
    EXPECT_NEAR(t.get_d(dx2), -2.443268609727591, 1e-15);
    EXPECT_NEAR(t.get_d(dx3), -2.3004865687368343, 1e-15);
    EXPECT_NEAR(t.get_d(dx4), -8.135610058930222, 1e-15);

    // std::cout.precision(std::numeric_limits<double>::max_digits10);
    // std::cout << t.get_d(dx) << std::endl;
    // std::cout << t.get_d(dx2) << std::endl;
    // std::cout << t.get_d(dx3) << std::endl;
    // std::cout << t.get_d(dx4) << std::endl;
}

// TEST(Tape2, TestJoan) {
//     auto [x, y, z] = Init<3>();
//     auto res = sin(sin(x) * log(x) * sin(y) + cos(z));
//     // auto res = sin(x) * log(x) * sin(y) + cos(z);

//     CalcTree ct(res);
//     ct.set(x) = 0.3;
//     ct.set(y) = 0.52;
//     ct.set(z) = 1.5;
//     ct.evaluate();

//     auto dx = d(x);
//     auto dy = d(y);
//     auto dz = d(z);

//     // auto dz2 = d(z) * d(z);
//     // auto dx2 = d(x) * d(x);
//     // auto dy2 = d(y) * d(y);
//     // auto dxy = d(x) * d(y);
//     // auto dzy = d(z) * d(y);
//     // auto dxz = d(x) * d(z);

//     auto dz3 = d(z) * d(z) * d(z);
//     auto dx3 = d(x) * d(x) * d(x);
//     // auto darb = pow<2>(d(x)) * d(z);

//     auto dres = d(res);
//     // auto dres2 = d<2>(res);
//     auto dres3 = d<3>(res);
//     // auto t = Tape2(dx, dy, dz, dres, dz2, dx2, dy2, dxy, dzy, dxz, dres2,
//     // dz3,
//     //                dx3, dres3);

//     auto t = Tape2(dx, dy, dz, dres, dz3, dx3, dres3);

//     t.set(dres) = 1.;
//     // t.set(dres2) = 1.;
//     t.set(dres3) = 1.;
//     t.backpropagate(ct);

//     std::cout.precision(std::numeric_limits<double>::max_digits10);
//     std::cout << t.get_d(dx) << std::endl;
//     std::cout << t.get_d(dy) << std::endl;
//     std::cout << t.get_d(dz) << std::endl;

//     // std::cout << t.get_d(dz2) << std::endl;
//     // std::cout << t.get_d(dx2) << std::endl;
//     // std::cout << t.get_d(dy2) << std::endl;
//     // std::cout << t.get_d(dxy) << std::endl;
//     // std::cout << t.get_d(dzy) << std::endl;
//     // std::cout << t.get_d(dxz) << std::endl;

//     std::cout << t.get_d(dz3) << std::endl;
//     std::cout << t.get_d(dx3) << std::endl;

//     // std::cout << t.get_d(darb) << std::endl;
//     // std::cout << t.get_d(dx2) << std::endl;
//     // std::cout << t.get_d(dy2) << std::endl;
//     // std::cout << t.get_d(dxy) << std::endl;
//     // std::cout << t.get_d(dx3) << std::endl;
//     // std::cout << t.get_d(dx2y1) << std::endl;
//     // std::cout << t.get_d(dx1y2) << std::endl;
//     // std::cout << t.get_d(dy3) << std::endl;
// }

// TEST(Tape2, ClangBreaking) {
//     auto [x, y, z] = Init<3>();
//     // auto res = sin(sin(x) * log(x) * sin(y) + cos(z));
//     auto res = sin(sin(x) * log(x) * sin(y) + cos(z));

//     CalcTree ct(res);
//     // // ct.set(x) = 0.3;
//     // // ct.set(y) = 0.52;
//     // ct.set(z) = 1.5;
//     // ct.evaluate();

//     // auto dx = d(x);
//     // auto dy = d(y);
//     // auto dz = d(z);
//     // auto dx2 = d(x) * d(x);
//     // auto dy2 = d(y) * d(y);
//     // auto dz2 = d(z) * d(z);
//     // auto dxy = d(x) * d(y);
//     // auto dxz = d(x) * d(z);
//     // auto dyz = d(y) * d(z);

//     auto dz3 = d(z) * d(z) * d(z);
//     // auto dx3 = d(x) * d(x) * d(x);
//     // auto darb = pow<2>(d(x)) * d(z);

//     // auto dres = d(res);
//     // auto dres2 = d<2>(res);
//     auto dres3 = d<3>(res);

//     // auto t = Tape2(dx, dy, dz, dres, dz3, dx3, darb, dres3);
//     // auto t = Tape2(dx, dy, dz, dx2, dy2, dz2, dxy, dxz, dyz, dz3, dres,
//     // dres2,
//     //                dres3);
//     auto t = Tape2(dz3, dres3);

//     // t.set(dres) = 1.;
//     // t.set(dres3) = 1.;
//     t.backpropagate(ct);

//     // from sympy import *
//     // x = Symbol('x')
//     // y = Symbol('y')
//     // z = Symbol('z')
//     // f = sin(sin(x) * log(x) * sin(y) + cos(z))
//     // valx = 0.3
//     // valy = 0.52
//     // valz = 1.5

//     // print(lambdify([x, y, z], f.diff(x))(valx, valy, valz))
//     // print(lambdify([x, y, z], f.diff(y))(valx, valy, valz))
//     // print(lambdify([x, y, z], f.diff(z))(valx, valy, valz))
//     // print(lambdify([x, y, z], f.diff(z).diff(z).diff(z))(valx, valy,
//     // valz))
//     // print(lambdify([x, y, z], f.diff(x).diff(x).diff(x))(valx, valy,
//     // valz))
//     // print(lambdify([x, y, z], f.diff(x).diff(x).diff(z))(valx, valy,
//     // valz))

//     // EXPECT_NEAR(t.get_d(dx), -0.0815897296146379, 1e-15);
//     // EXPECT_NEAR(t.get_d(dy), -0.3070338475892694, 1e-15);
//     // EXPECT_NEAR(t.get_d(dz), -0.9918908259808428, 1e-15);
//     // EXPECT_NEAR(t.get_d(dz3), 2.0012254874538344, 1e-15);
//     // EXPECT_NEAR(t.get_d(dx3), -5.854071971576513, 1e-15);
//     // EXPECT_NEAR(t.get_d(darb), -0.1738608795967254, 1e-15);

//     // std::cout.precision(std::numeric_limits<double>::max_digits10);
//     // std::cout << t.get_d(dx) << std::endl;
//     // std::cout << t.get_d(dy) << std::endl;
//     // std::cout << t.get_d(dz) << std::endl;
//     // std::cout << t.get_d(dz3) << std::endl;
//     // std::cout << t.get_d(dx3) << std::endl;
//     // std::cout << t.get_d(darb) << std::endl;
// }

TEST(Tape2, BSDerivatives) {
    auto [S, K, v, T] = Init<4>();

    auto res = call_price(S, K, v, T);

    CalcTree ct(res);
    ct.set(S) = 100.0;
    ct.set(K) = 102.0;
    ct.set(v) = 0.15;
    ct.set(T) = 0.5;
    ct.evaluate();

    auto dS = d(S);
    auto dK = d(K);
    auto dv = d(v);
    auto dT = d(T);
    auto d2S = pow<2>(d(S));
    auto d2K = pow<2>(d(K));
    auto d2v = pow<2>(d(v));
    auto d2T = pow<2>(d(T));
    auto dSK = d(S) * d(K);
    auto dSv = d(S) * d(v);
    auto dST = d(S) * d(T);
    auto dKv = d(K) * d(v);
    auto dKT = d(K) * d(T);
    auto dvT = d(v) * d(T);

    auto dr1 = d(res);
    auto dr2 = d<2>(res);
    auto dr22 = pow<2>(d(res));

    auto t = Tape2(dS, dK, dv, dT, d2S, d2K, d2v, d2T, dSK, dSv, dST, dKv, dKT,
                   dvT, dr1, dr2);

    t.set(dr1) = 1.;
    t.set(dr2) = 1.;
    t.backpropagate(ct);

    /*from sympy import *
    S = Symbol('S')
    K = Symbol('K')
    v = Symbol('v')
    T = Symbol('T')

    totalvol = v * sqrt(T)
    d1 = log(S / K) / totalvol + totalvol * 0.5
    d2 = d1 + totalvol

    minus_one_over_root_two =-1.0 / sqrt(2.0)
    cdf_n_d1 = 0.5 * erfc(d1 * minus_one_over_root_two)
    cdf_n_d2 = 0.5 * erfc(d2 * minus_one_over_root_two)

    f = S * cdf_n_d1 - K * cdf_n_d2

    valS = 100.0
    valK = 102.0
    valv = 0.15
    valT = 0.5

    print(lambdify([S, K, v, T], f.diff(S))(valS, valK, valv, valT))
    print(lambdify([S, K, v, T], f.diff(K))(valS, valK, valv, valT))
    print(lambdify([S, K, v, T], f.diff(v))(valS, valK, valv, valT))
    print(lambdify([S, K, v, T], f.diff(T))(valS, valK, valv, valT))
    print(lambdify([S, K, v, T], f.diff(S).diff(S))(valS, valK, valv,valT))
    print(lambdify([S, K, v, T], f.diff(K).diff(K))(valS, valK, valv,valT))
    print(lambdify([S, K, v, T], f.diff(v).diff(v))(valS, valK, valv,valT))
    print(lambdify([S, K, v, T], f.diff(T).diff(T))(valS, valK, valv,valT))
    print(lambdify([S, K, v, T], f.diff(S).diff(K))(valS, valK, valv,valT))
    print(lambdify([S, K, v, T], f.diff(S).diff(v))(valS, valK, valv,valT))
    print(lambdify([S, K, v, T], f.diff(S).diff(T))(valS, valK, valv,valT))
    print(lambdify([S, K, v, T], f.diff(K).diff(v))(valS, valK, valv,valT))
    print(lambdify([S, K, v, T], f.diff(K).diff(T))(valS, valK, valv,valT))
    print(lambdify([S, K, v, T], f.diff(v).diff(T))(valS, valK, valv,valT))*/

    EXPECT_NEAR(t.get_d(dS), 0.33961663008862164, 1e-15);
    EXPECT_NEAR(t.get_d(dK), -0.3838761485986675, 1e-14);
    EXPECT_NEAR(t.get_d(dv), -30.580208040388406, 1e-13);
    EXPECT_NEAR(t.get_d(dT), -4.587031206058272, 1e-14);
    EXPECT_NEAR(t.get_d(d2S), 0.11262750926223955, 1e-15);
    EXPECT_NEAR(t.get_d(d2K), 0.10825404581145669, 1e-15);
    EXPECT_NEAR(t.get_d(d2v), 26.405582264938744, 5e-13);
    EXPECT_NEAR(t.get_d(d2T), 5.181156807019395, 5e-14);
    EXPECT_NEAR(t.get_d(dSK), -0.11041912672768583, 5e-13);
    EXPECT_NEAR(t.get_d(dSv), 1.8990518775466825, 5e-13);
    EXPECT_NEAR(t.get_d(dST), 0.28485778163200237, 5e-13);
    EXPECT_NEAR(t.get_d(dKv), -2.161621527402506, 5e-13);
    EXPECT_NEAR(t.get_d(dKT), -0.3242432291103765, 5e-13);
    EXPECT_NEAR(t.get_d(dvT), -26.61937070064767, 5e-13);
}

TEST(Tape2, BSPrimalTime) {
    std::mt19937 generator(123);
    std::uniform_real_distribution<double> stock_distr(90, 110.0);
    std::uniform_real_distribution<double> vol_distr(0.05, 0.3);
    std::uniform_real_distribution<double> time_distr(0.5, 1.5);

    double acc = 0;
    std::size_t iters = 1;
    if (auto env_p = std::getenv("ITERATIONS")) {
        iters = std::stoul(env_p);
    }

    auto [S, K, v, T] = Init<4>();

    auto res = call_price(S, K, v, T);

    using namespace std::chrono;
    time_point<high_resolution_clock> start = high_resolution_clock::now();

    for (std::size_t j = 0; j < iters; ++j) {
        double S = stock_distr(generator);
        double K = stock_distr(generator);
        double v = vol_distr(generator);
        double T = time_distr(generator);
        auto res = call_price(S, K, v, T);
        acc += res;
    }

    std::cout << duration_cast<microseconds>(high_resolution_clock::now() -
                                             start)
                     .count()
              << std::endl;

    acc /= static_cast<double>(iters);
    std::cout << acc << std::endl;
}

TEST(Tape2, BSDerivatives1Time) {
    std::mt19937 generator(123);
    std::uniform_real_distribution<double> stock_distr(90, 110.0);
    std::uniform_real_distribution<double> vol_distr(0.05, 0.3);
    std::uniform_real_distribution<double> time_distr(0.5, 1.5);

    double acc = 0;
    std::size_t iters = 1;
    if (auto env_p = std::getenv("ITERATIONS")) {
        iters = std::stoul(env_p);
    }

    auto [S, K, v, T] = Init<4>();

    auto res = call_price(S, K, v, T);

    CalcTree ct(res);

    auto dS = d(S);
    auto dK = d(K);
    auto dv = d(v);
    auto dT = d(T);

    auto dr = d(res);

    auto t = Tape2(dr, dS, dK, dv, dT);

    using namespace std::chrono;
    time_point<high_resolution_clock> start = high_resolution_clock::now();

    for (std::size_t j = 0; j < iters; ++j) {
        ct.set(S) = stock_distr(generator);
        ct.set(K) = stock_distr(generator);
        ct.set(v) = vol_distr(generator);
        ct.set(T) = time_distr(generator);
        ct.evaluate();

        t.reset_der();
        t.set(dr) = 1.;
        t.backpropagate(ct);

        acc += ct.val(res);
    }

    std::cout << duration_cast<microseconds>(high_resolution_clock::now() -
                                             start)
                     .count()
              << std::endl;

    acc /= static_cast<double>(iters);
    std::cout << acc << std::endl;
}

TEST(Tape2, BSDerivatives2Time) {
    std::mt19937 generator(123);
    std::uniform_real_distribution<double> stock_distr(90, 110.0);
    std::uniform_real_distribution<double> vol_distr(0.05, 0.3);
    std::uniform_real_distribution<double> time_distr(0.5, 1.5);

    double acc = 0;
    std::size_t iters = 1;
    if (auto env_p = std::getenv("ITERATIONS")) {
        iters = std::stoul(env_p);
    }

    auto [S, K, v, T] = Init<4>();

    auto res = call_price(S, K, v, T);

    CalcTree ct(res);

    auto dS = d<2>(S);
    auto dK = d<2>(K);
    auto dv = d<2>(v);
    auto dT = d<2>(T);
    auto d2S = pow<2>(d(S));
    auto d2K = pow<2>(d(K));
    auto d2v = pow<2>(d(v));
    auto d2T = pow<2>(d(T));
    auto dSK = d(S) * d(K);
    auto dSv = d(S) * d(v);
    auto dST = d(S) * d(T);
    auto dKv = d(K) * d(v);
    auto dKT = d(K) * d(T);
    auto dvT = d(v) * d(T);

    auto dr2 = d<2>(res);

    auto t = Tape2(d2S, dS, dK, dv, dT, d2K, d2v, d2T, dSK, dSv, dST, dKv, dKT,
                   dvT, dr2);

    using namespace std::chrono;
    time_point<high_resolution_clock> start = high_resolution_clock::now();

    for (std::size_t j = 0; j < iters; ++j) {
        ct.set(S) = stock_distr(generator);
        ct.set(K) = stock_distr(generator);
        ct.set(v) = vol_distr(generator);
        ct.set(T) = time_distr(generator);
        ct.evaluate();

        t.reset_der();
        t.set(dr2) = 1.;
        t.backpropagate(ct);

        acc += ct.val(res);
    }

    std::cout << duration_cast<microseconds>(high_resolution_clock::now() -
                                             start)
                     .count()
              << std::endl;

    acc /= static_cast<double>(iters);
    std::cout << acc << std::endl;
    /*from sympy import *
    S = Symbol('S')
    K = Symbol('K')
    v = Symbol('v')
    T = Symbol('T')

    totalvol = v * sqrt(T)
    d1 = log(S / K) / totalvol + totalvol * 0.5
    d2 = d1 + totalvol

    minus_one_over_root_two =-1.0 / sqrt(2.0)
    cdf_n_d1 = 0.5 * erfc(d1 * minus_one_over_root_two)
    cdf_n_d2 = 0.5 * erfc(d2 * minus_one_over_root_two)

    f = S * cdf_n_d1 - K * cdf_n_d2

    valS = 100.0
    valK = 102.0
    valv = 0.15
    valT = 0.5

    print(lambdify([S, K, v, T], f.diff(S))(valS, valK, valv, valT))
    print(lambdify([S, K, v, T], f.diff(K))(valS, valK, valv, valT))
    print(lambdify([S, K, v, T], f.diff(v))(valS, valK, valv, valT))
    print(lambdify([S, K, v, T], f.diff(T))(valS, valK, valv, valT))
    print(lambdify([S, K, v, T], f.diff(S).diff(S))(valS, valK, valv,valT))
    print(lambdify([S, K, v, T], f.diff(K).diff(K))(valS, valK, valv,valT))
    print(lambdify([S, K, v, T], f.diff(v).diff(v))(valS, valK, valv,valT))
    print(lambdify([S, K, v, T], f.diff(T).diff(T))(valS, valK, valv,valT))
    print(lambdify([S, K, v, T], f.diff(S).diff(K))(valS, valK, valv,valT))
    print(lambdify([S, K, v, T], f.diff(S).diff(v))(valS, valK, valv,valT))
    print(lambdify([S, K, v, T], f.diff(S).diff(T))(valS, valK, valv,valT))
    print(lambdify([S, K, v, T], f.diff(K).diff(v))(valS, valK, valv,valT))
    print(lambdify([S, K, v, T], f.diff(K).diff(T))(valS, valK, valv,valT))
    print(lambdify([S, K, v, T], f.diff(v).diff(T))(valS, valK, valv,valT))*/

    // EXPECT_NEAR(t.get_d(dS), 0.33961663008862164, 1e-15);
    // EXPECT_NEAR(t.get_d(dK), -0.3838761485986675, 1e-14);
    // EXPECT_NEAR(t.get_d(dv), -30.580208040388406, 1e-13);
    // EXPECT_NEAR(t.get_d(dT), -4.587031206058272, 1e-14);
    // EXPECT_NEAR(t.get_d(d2S), 0.11262750926223955, 1e-15);
    // EXPECT_NEAR(t.get_d(d2K), 0.10825404581145669, 1e-15);
    // EXPECT_NEAR(t.get_d(d2v), 26.405582264938744, 5e-13);
    // EXPECT_NEAR(t.get_d(d2T), 5.181156807019395, 5e-14);
    // EXPECT_NEAR(t.get_d(dSK), -0.11041912672768583, 5e-13);
    // EXPECT_NEAR(t.get_d(dSv), 1.8990518775466825, 5e-13);
    // EXPECT_NEAR(t.get_d(dST), 0.28485778163200237, 5e-13);
    // EXPECT_NEAR(t.get_d(dKv), -2.161621527402506, 5e-13);
    // EXPECT_NEAR(t.get_d(dKT), -0.3242432291103765, 5e-13);
    // EXPECT_NEAR(t.get_d(dvT), -26.61937070064767, 5e-13);
}

// TEST(Tape2, TapeAndTreeUnivariateMultipleOutput) {
//     double const valx = 0.5;
//     double const valy = 0.44;

//     auto [x, y] = Init<2>();
//     auto res1 = log(log(cos(x) * sin(y)));
//     auto res2 = erfc(log(cos(x) * sin(y)));
//     CalcTree ct(res1, res2);
//     ct.set(x) = valx;
//     ct.set(y) = valy;
//     ct.evaluate();

//     auto dx = d(x);
//     auto dy = d(y);
//     auto t = Tape2(dx, dy, dres1, dres2, dx2, dy2, dxy, dres1_2, dres2_2,
//     dx3,
//                    dx2y1, dx1y2, dy3, dres1_3, dres2_3);

//     auto dx2 = pow<2>(d(x));
//     auto dy2 = pow<2>(d(y));
//     auto dxy = d(x) * d(y);

//     auto dx3 = pow<3>(d(x));
//     auto dx2y1 = pow<2>(d(x)) * d(y);
//     auto dx1y2 = d(x) * pow<2>(d(y));
//     auto dy3 = pow<3>(d(y));

//     auto dres1 = d(res1);
//     auto dres2 = d(res2);

//     auto dres1_2 = d<2>(res1);
//     auto dres2_2 = d<2>(res2);

//     auto dres1_3 = d<3>(res1);
//     // auto dres1_2_dres2_1 = d<2>(res1) * d(res2);
//     // auto dres1_1_dres2_2 = d(res1) * d<2>(res2);
//     auto dres2_3 = d<3>(res2);

//     // auto t = Tape2(dx, dy, dres1, dres2);
//     auto t = Tape2(dx, dy, dres1, dres2, dx2, dy2, dxy, dres1_2, dres2_2,
//     dx3,
//                    dx2y1, dx1y2, dy3, dres1_3, dres2_3);
//     std::cout.precision(std::numeric_limits<double>::max_digits10);

//     t.set(dres1) = 1.;
//     t.set(dres2) = 0.;

//     t.set(dres1_2) = 1.;
//     t.set(dres2_2) = 0.;

//     t.set(dres1_3) = 1.;
//     t.set(dres2_3) = 0.;

//     t.backpropagate(ct);
//     std::cout << t.get_d(dx) << std::endl;
//     std::cout << t.get_d(dy) << std::endl;

//     std::cout << t.get_d(dx2) << std::endl;
//     std::cout << t.get_d(dy2) << std::endl;
//     std::cout << t.get_d(dxy) << std::endl;

//     std::cout << t.get_d(dx3) << std::endl;
//     std::cout << t.get_d(dx2y1) << std::endl;
//     std::cout << t.get_d(dx1y2) << std::endl;
//     std::cout << t.get_d(dy3) << std::endl;

//     t.reset_der();

//     t.set(dres1) = 0.;
//     t.set(dres2) = 1.;

//     t.set(dres1_2) = 0.;
//     t.set(dres2_2) = 1.;

//     t.set(dres1_3) = 0.;
//     t.set(dres2_3) = 1.;
//     t.backpropagate(ct);
//     std::cout << t.get_d(dx) << std::endl;
//     std::cout << t.get_d(dy) << std::endl;

//     std::cout << t.get_d(dx2) << std::endl;
//     std::cout << t.get_d(dy2) << std::endl;
//     std::cout << t.get_d(dxy) << std::endl;

//     std::cout << t.get_d(dx3) << std::endl;
//     std::cout << t.get_d(dx2y1) << std::endl;
//     std::cout << t.get_d(dx1y2) << std::endl;
//     std::cout << t.get_d(dy3) << std::endl;
// }

// TEST(Tape2, TapeAndTreeUnivariateMultipleOutputBreaking) {
//     double const valx = 0.5;
//     double const valy = 0.44;

//     auto [x, y] = Init<2>();
//     auto res = log(log(cos(x) * sin(y)));
//     auto res2 = erfc(log(cos(x) * sin(y)));
//     CalcTree ct(res, res2);
//     ct.set(x) = valx;
//     ct.set(y) = valy;
//     ct.evaluate();

//     auto dx = d(x);
//     auto dy = d(y);

//     auto dx2 = pow<2>(d(x));
//     auto dy2 = pow<2>(d(y));
//     auto dxy = d(x) * d(y);

//     auto dx3 = pow<3>(d(x));
//     auto dx2y1 = pow<2>(d(x)) * d(y);
//     auto dx1y2 = d(x) * pow<2>(d(y));
//     auto dy3 = pow<3>(d(y));

//     auto dres = d(res);
//     auto dres1_2 = d<2>(res);
//     auto dres3 = d<3>(res);

//     auto dres2_2 = d<2>(res2);

//     auto t = Tape2(dx3, dx2y1, dx1y2, dy3, dx2, dx, dy2, dxy, dy, dres,
//     dres1_2,
//                    dres3, dres1_2, dres2_2);
//     t.set(dres) = 1.;
//     t.set(dres1_2) = 0.;
//     t.set(dres3) = 1.;
//     t.set(dres2_2) = 1.;
//     t.backpropagate(ct);
// }

} // namespace adhoc3
