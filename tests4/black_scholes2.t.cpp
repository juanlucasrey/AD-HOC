#include <adhoc.hpp>
#include <backpropagator.hpp>
#include <calc_tree.hpp>
#include <differential_operator.hpp>

#include <gtest/gtest.h>

namespace adhoc4 {

namespace {

template <class I1, class I2, class I3, class I4>
auto call_price2(const I1 &S, const I2 &K, const I3 &v, const I4 &T) {
    using constants::CD;
    using std::erfc;
    using std::log;
    using std::sqrt;
    auto totalvol = v * sqrt(T);
    auto d1 = log(S * sqrt(K)) * log(totalvol) + totalvol * CD<0.5>();
    auto d2 = d1 + totalvol;
    return S * erfc(d1) + K * erfc(d2);
}

} // namespace

TEST(BlackScholes, BS2) {
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

    // from sympy import *
    // S = Symbol('S')
    // K = Symbol('K')
    // v = Symbol('v')
    // T = Symbol('T')

    // totalvol = v * sqrt(T)
    // d1 = log(S * sqrt(K)) * log(totalvol) + totalvol * 0.5
    // d2 = d1 + totalvol
    // result = S * erfc(d1) + K * erfc(d2)

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

    EXPECT_NEAR(ct.get(res), 1.8890465506626162, 1e-14);
    EXPECT_NEAR(t.get(dS), 6.0454412443913581, 1e-14);
    EXPECT_NEAR(t.get(dK), 3.3740304502248355, 1e-14);
    EXPECT_NEAR(t.get(dv), -1.9089022751421003, 1e-14);
    EXPECT_NEAR(t.get(dT), -0.28633534127131505, 1e-14);
    EXPECT_NEAR(t.get(dSS), 1.38790636613215, 1e-14);
    EXPECT_NEAR(t.get(dKK), 0.334319305617748, 1e-14);
    EXPECT_NEAR(t.get(dvv), 2.38267251329635, 1e-14);
    EXPECT_NEAR(t.get(dTT), 0.339945472820483, 1e-14);
    EXPECT_NEAR(t.get(dSK), 4.3981330899714148, 1e-14);
    EXPECT_NEAR(t.get(dSv), -16.278830190533597, 1e-14);
    EXPECT_NEAR(t.get(dST), -2.4418245285800397, 1e-14);
    EXPECT_NEAR(t.get(dKv), -9.1167961040278751, 1e-14);
    EXPECT_NEAR(t.get(dKT), -1.3675194156041814, 1e-14);
    EXPECT_NEAR(t.get(dvT), -1.55150139814765, 1e-14);
}

} // namespace adhoc4
