#include <adhoc.hpp>
#include <backpropagator.hpp>
#include <calc_tree.hpp>
#include <differential_operator.hpp>

#include "call_price.hpp"
#include <gtest/gtest.h>

namespace adhoc4 {

TEST(BlackScholes, Order2) {
    ADHOC(S);
    ADHOC(K);
    ADHOC(v);
    ADHOC(T);

    auto res = call_price2(S, K, v, T);

    CalcTree ct(res);
    ct.set(S) = 1.01;
    ct.set(K) = 1.02;
    ct.set(v) = 0.15;
    ct.set(T) = 0.5;
    ct.evaluate();

    // order 1
    auto dS = d(S);
    auto dK = d(K);
    auto dv = d(v);
    auto dT = d(T);

    // order 2
    auto dSS = d<2>(S);
    auto dKK = d<2>(K);
    auto dvv = d<2>(v);
    auto dTT = d<2>(T);
    auto dSK = d(S) * d(K);
    auto dSv = d(S) * d(v);
    auto dST = d(S) * d(T);
    auto dKv = d(K) * d(v);
    auto dKT = d(K) * d(T);
    auto dvT = d(v) * d(T);

    auto dres = d(res);

    BackPropagator t(dS, dK, dv, dT, dSS, dKK, dvv, dTT, dSK, dSv, dST, dKv,
                     dKT, dvT, dres);

    t.set(dres) = 1.;
    t.backpropagate(ct);

    EXPECT_NEAR(ct.get(res), 1.02586812693369, 1e-14);
    EXPECT_NEAR(t.get(dS), 8.03256922212305, 1e-14);
    EXPECT_NEAR(t.get(dK), -6.94806547785353, 1e-14);
    EXPECT_NEAR(t.get(dv), 1.07376087149468, 1e-14);
    EXPECT_NEAR(t.get(dT), 0.161064130724201, 1e-14);
    EXPECT_NEAR(t.get(dSS), -0.975928505799621, 1e-14);
    EXPECT_NEAR(t.get(dKK), -0.956886455945976, 1e-14);
    EXPECT_NEAR(t.get(dvv), -6.73743074286942, 1e-14);
    EXPECT_NEAR(t.get(dTT), -0.312656322438763, 1e-14);
    EXPECT_NEAR(t.get(dSK), 0.966360579272174, 1e-14);
    EXPECT_NEAR(t.get(dSv), -50.2111752240705, 2e-14);
    EXPECT_NEAR(t.get(dST), -7.53167628361058, 1e-14);
    EXPECT_NEAR(t.get(dKv), 50.7716155370646, 1e-14);
    EXPECT_NEAR(t.get(dKT), 7.61574233055969, 1e-14);
    EXPECT_NEAR(t.get(dvT), 0.0631462600642643, 1e-14);

    // from sympy import *
    // S = Symbol('S')
    // K = Symbol('K')
    // v = Symbol('v')
    // T = Symbol('T')

    // def cdf(x):
    //     return 0.5 * erfc(x * -0.70710678118654746)

    // totalvol = v * sqrt(T)
    // d1 = log(S / K) / totalvol + totalvol * 0.5
    // d2 = d1 + totalvol
    // result = S * cdf(d1) + K * cdf(d2)

    // def diff(f, diffs):
    //     fdiff = f
    //     for i in diffs:
    //         fdiff = fdiff.diff(i)

    //     return fdiff

    // S0 = 1.01
    // K0 = 1.02
    // v0 = 0.15
    // T0 = 0.5

    // print(lambdify([S, K, v, T], result)(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (S,)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (K,)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (v,)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (T,)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (S, S)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (K, K)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (v, v)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (T, T)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (S, K)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (S, v)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (S, T)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (K, v)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (K, T)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (v, T)))(S0, K0, v0, T0))
}

} // namespace adhoc4
