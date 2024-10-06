#include <adhoc.hpp>
#include <backpropagator.hpp>
#include <calc_tree.hpp>
#include <differential_operator.hpp>
#include <utils/bivariate.hpp>

#include <utils/tuple.hpp>

#include <sort.hpp>

#include "type_name.hpp"

#include <gtest/gtest.h>

namespace adhoc4 {

TEST(BackPropagator2, UnivariateSingle) {
    ADHOC(x);
    auto res = erfc(x);

    CalcTree ct(res);
    ct.set(x) = 1.2;
    ct.evaluate();

    BackPropagator bp(d(res), d(x), d<2>(x), d<3>(x));
    bp.set(d(res)) = 1.;
    bp.backpropagate2(ct);

    // from sympy import *
    // x = Symbol('x')
    // f = erfc(x)

    // def taylor_coeff(f,x,x0, n):
    //     fdiff = f
    //     for i in range(0,n):
    //         fdiff = fdiff.diff(x)

    //     val = lambdify([x], fdiff)(x0)
    //     val = val / factorial(n)
    //     return val

    // print(taylor_coeff(f, x, 1.2, 1))
    // print(taylor_coeff(f, x, 1.2, 2))
    // print(taylor_coeff(f, x, 1.2, 3))

    EXPECT_NEAR(bp.get(d(x)), -0.267344347003539, 1e-14);
    EXPECT_NEAR(bp.get(d<2>(x)), 0.320813216404247, 1e-14);
    EXPECT_NEAR(bp.get(d<3>(x)), -0.167535790788885, 1e-14);
}

TEST(BackPropagator2, UnivariateDouble) {
    ADHOC(x);
    auto res = log(erfc(x));

    CalcTree ct(res);
    ct.set(x) = 1.2;
    ct.evaluate();

    BackPropagator bp(d(res), d(x), d<2>(x), d<3>(x));
    bp.set(d(res)) = 1.;
    bp.backpropagate2(ct);

    // from sympy import *
    // x = Symbol('x')
    // f = log(erfc(x))

    // def taylor_coeff(f,x,x0, n):
    //     fdiff = f
    //     for i in range(0,n):
    //         fdiff = fdiff.diff(x)

    //     val = lambdify([x], fdiff)(x0)
    //     val = val / factorial(n)
    //     return val

    // print(taylor_coeff(f, x, 1.2, 1))
    // print(taylor_coeff(f, x, 1.2, 2))
    // print(taylor_coeff(f, x, 1.2, 3))

    EXPECT_NEAR(bp.get(d(x)), -2.98089202449025, 1e-14);
    EXPECT_NEAR(bp.get(d<2>(x)), -0.865788201446497, 1e-14);
    EXPECT_NEAR(bp.get(d<3>(x)), -0.0342861937391300, 1e-14);
}

TEST(BackPropagator2, SumSingle) {
    ADHOC(x);
    ADHOC(y);
    // auto res = x + y;
    auto res = log(x + y);

    CalcTree ct(res);
    ct.set(x) = 1.2;
    ct.set(y) = 0.4;
    ct.evaluate();

    BackPropagator bp(d(res), d(x), d(y), d(x) * d(y), d<5>(x));
    bp.set(d(res)) = 1.;
    bp.backpropagate2(ct);

    EXPECT_NEAR(bp.get(d(x)), 0.625, 1e-13);
    EXPECT_NEAR(bp.get(d(y)), 0.625, 1e-13);
    EXPECT_NEAR(bp.get(d(x) * d(y)), -0.390625, 1e-13);
}

TEST(BackPropagator2, MulSingle) {
    ADHOC(x);
    ADHOC(y);
    auto res = log(x * y);

    CalcTree ct(res);
    ct.set(x) = 1.2;
    ct.set(y) = 0.4;
    ct.evaluate();

    BackPropagator bp(d(res), d(x), d(y), d(x) * d(y), d<5>(x));
    bp.set(d(res)) = 1.;
    bp.backpropagate2(ct);

    EXPECT_NEAR(bp.get(d(x)), 0.833333333333333, 1e-14);
    EXPECT_NEAR(bp.get(d(y)), 2.5, 1e-14);
    EXPECT_NEAR(bp.get(d(x) * d(y)), 0., 1e-14);
    EXPECT_NEAR(bp.get(d<5>(x)), 0.0803755144032922, 1e-14);
}

TEST(BackPropagator2, MulSingle2) {
    ADHOC(x);
    ADHOC(y);
    auto res = erfc(x * y);

    CalcTree ct(res);
    ct.set(x) = 1.2;
    ct.set(y) = 0.4;
    ct.evaluate();

    BackPropagator bp(d(res), d(x), d(y), d(x) * d(y), d<5>(x));
    bp.set(d(res)) = 1.;
    bp.backpropagate2(ct);

    // std::cout.precision(std::numeric_limits<double>::max_digits10);
    // std::cout << bp.get(d(x)) << std::endl;
    // std::cout << bp.get(d(y)) << std::endl;
    // std::cout << bp.get(d(x) * d(y)) << std::endl;
    // std::cout << bp.get(d<5>(x)) << std::endl;

    EXPECT_NEAR(bp.get(d(x)), -0.358470648907805, 1e-14);
    EXPECT_NEAR(bp.get(d(y)), -1.07541194672341, 1e-14);
    EXPECT_NEAR(bp.get(d(x) * d(y)), -0.483218434727721, 1e-14);
    EXPECT_NEAR(bp.get(d<5>(x)), -0.000136899199787365, 1e-14);
}

TEST(BackPropagator2, DivSingle) {
    ADHOC(x);
    ADHOC(y);
    auto res = erfc(x / y);

    CalcTree ct(res);
    ct.set(x) = 1.2;
    ct.set(y) = 0.4;
    ct.evaluate();

    BackPropagator bp(d(res), d(x), d(y), d(x) * d(y), d<5>(x));
    bp.set(d(res)) = 1.;
    bp.backpropagate2(ct);

    EXPECT_NEAR(bp.get(d(x)), -0.000348132629866870, 1e-14);
    EXPECT_NEAR(bp.get(d(y)), 0.00104439788960061, 1e-14);
    EXPECT_NEAR(bp.get(d(x) * d(y)), -0.0147956367693420, 1e-14);
    EXPECT_NEAR(bp.get(d<5>(x)), -0.0992721952354747, 1e-14);
}

} // namespace adhoc4
