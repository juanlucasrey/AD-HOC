#include <adhoc.hpp>
#include <backpropagator.hpp>
#include <calc_tree.hpp>
#include <differential_operator.hpp>

#include "call_price.hpp"
#include "test_tools.hpp"

int main() {
    TEST_START;
    using namespace adhoc4;
    ADHOC(S);
    ADHOC(K);
    ADHOC(v);
    ADHOC(T);

    auto res = call_price(S, K, v, T);

    CalcTree ct(res);
    ct.set(S) = 104.;
    ct.set(K) = 98.;
    ct.set(v) = 0.22;
    ct.set(T) = 1.25;
    ct.evaluate();

    // order 1
    auto dS = d(S);
    auto dK = d(K);
    auto dv = d(v);
    auto dT = d(T);

    auto dres = d(res);

    BackPropagator bp(dS, dK, dv, dT, dres);

    bp.set(dres) = 1.;
    bp.backpropagate(ct);

    EXPECT_NEAR(ct.get(res), 13.171437637423516, 1e-13);
    EXPECT_NEAR(bp.get(dS), 0.6422853983194291, 1e-14);
    EXPECT_NEAR(bp.get(dK), -0.54720656926323585, 1e-14);
    EXPECT_NEAR(bp.get(dv), 43.404665363383245, 1e-14);
    EXPECT_NEAR(bp.get(dT), 3.8196105519777257, 1e-14);

    // from sympy import *
    // S = Symbol('S')
    // K = Symbol('K')
    // v = Symbol('v')
    // T = Symbol('T')

    // def cdf(x):
    //     return 0.5 * erfc(x * -0.70710678118654746)

    // totalvol = v * sqrt(T)
    // d1 = log(S / K) / totalvol + totalvol * 0.5
    // d2 = d1 - totalvol
    // result = S * cdf(d1) - K * cdf(d2)

    // def diff(f, diffs):
    //     fdiff = f
    //     for i in diffs:
    //         fdiff = fdiff.diff(i)

    //     return fdiff

    // S0 = 104.
    // K0 = 98.
    // v0 = 0.22
    // T0 = 1.25

    // print(lambdify([S, K, v, T], result)(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (S,)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (K,)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (v,)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (T,)))(S0, K0, v0, T0))

    TEST_END;
}
