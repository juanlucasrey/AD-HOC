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

    // order 2
    auto dSS = d<2>(S);
    auto dSK = d(S) * d(K);
    auto dSv = d(S) * d(v);
    auto dST = d(S) * d(T);
    auto dKK = d<2>(K);
    auto dKv = d(K) * d(v);
    auto dKT = d(K) * d(T);
    auto dvv = d<2>(v);
    auto dvT = d(v) * d(T);
    auto dTT = d<2>(T);

    // order 3
    auto dSSS = d<3>(S);
    auto dSSK = d<2>(S) * d(K);
    auto dSSv = d<2>(S) * d(v);
    auto dSST = d<2>(S) * d(T);
    auto dSKK = d(S) * d<2>(K);
    auto dSKv = d(S) * d(K) * d(v);
    auto dSKT = d(S) * d(K) * d(T);
    auto dSvv = d(S) * d<2>(v);
    auto dSvT = d(S) * d(v) * d(T);
    auto dSTT = d(S) * d<2>(T);
    auto dKKK = d<3>(K);
    auto dKKv = d<2>(K) * d(v);
    auto dKKT = d<2>(K) * d(T);
    auto dKvv = d(K) * d<2>(v);
    auto dKvT = d(K) * d(v) * d(T);
    auto dKTT = d(K) * d<2>(T);
    auto dvvv = d<3>(v);
    auto dvvT = d<2>(v) * d(T);
    auto dvTT = d(v) * d<2>(T);
    auto dTTT = d<3>(T);

    auto dres = d(res);

    BackPropagator bp(dS, dK, dv, dT, dSS, dSK, dSv, dST, dKK, dKv, dKT, dvv,
                      dvT, dTT, dSSS, dSSK, dSSv, dSST, dSKK, dSKv, dSKT, dSvv,
                      dSvT, dSTT, dKKK, dKKv, dKKT, dKvv, dKvT, dKTT, dvvv,
                      dvvT, dvTT, dTTT, dres);

    bp.set(dres) = 1.;
    bp.backpropagate(ct);

    EXPECT_NEAR(ct.get(res), 13.171437637423516, 1e-13);
    EXPECT_NEAR(bp.get(dS), 0.6422853983194291, 1e-14);
    EXPECT_NEAR(bp.get(dK), -0.54720656926323585, 1e-14);
    EXPECT_NEAR(bp.get(dv), 43.404665363383245, 1e-14);
    EXPECT_NEAR(bp.get(dT), 3.8196105519777257, 1e-14);
    EXPECT_NEAR(bp.get(dSS), 0.014592746558426314, 1e-14);
    EXPECT_NEAR(bp.get(dSK), -0.015486180021187103, 1e-14);
    EXPECT_NEAR(bp.get(dSv), -0.20124961109986916, 2e-14);
    EXPECT_NEAR(bp.get(dST), -0.017709965776788628, 1e-14);
    EXPECT_NEAR(bp.get(dKK), 0.01643431349187204, 1e-14);
    EXPECT_NEAR(bp.get(dKv), 0.6564757644670376, 1e-14);
    EXPECT_NEAR(bp.get(dKT), 0.057769867273099296, 1e-14);
    EXPECT_NEAR(bp.get(dvv), 8.53118665253345, 1e-13);
    EXPECT_NEAR(bp.get(dvT), 18.112610570776255, 1e-13);
    EXPECT_NEAR(bp.get(dTT), -1.461778711353871, 1e-14);
    EXPECT_NEAR(bp.get(dSSS), -0.00034829031543903103, 1e-14);
    EXPECT_NEAR(bp.get(dSSK), 0.00022070863517584553, 1e-14);
    EXPECT_NEAR(bp.get(dSSv), -0.06346246194775129, 2e-14);
    EXPECT_NEAR(bp.get(dSST), -0.005584696651402094, 1e-14);
    EXPECT_NEAR(bp.get(dSKK), -7.619916364388522e-05, 1e-14);
    EXPECT_NEAR(bp.get(dSKv), 0.06734791880169519, 1e-14);
    EXPECT_NEAR(bp.get(dSKT), 0.00592661685454916, 1e-14);
    EXPECT_NEAR(bp.get(dSvv), 3.6870433621286907, 1e-13);
    EXPECT_NEAR(bp.get(dSvT), 0.24395997142737658, 1e-13);
    EXPECT_NEAR(bp.get(dSTT), 0.035636450107039974, 1e-14);
    EXPECT_NEAR(bp.get(dKKK), -8.683265788681426e-05, 1e-14);
    EXPECT_NEAR(bp.get(dKKv), -0.07147126076914595, 1e-14);
    EXPECT_NEAR(bp.get(dKKT), -0.006289470947684829, 2e-14);
    EXPECT_NEAR(bp.get(dKvv), -3.8257277858046095, 5e-14);
    EXPECT_NEAR(bp.get(dKvT), -0.07407373936398809, 1e-14);
    EXPECT_NEAR(bp.get(dKTT), -0.05273438288251053, 1e-14);
    EXPECT_NEAR(bp.get(dvvv), -168.913390475117, 1e-11);
    EXPECT_NEAR(bp.get(dvvT), -8.039429039783727, 1e-13);
    EXPECT_NEAR(bp.get(dvTT), -7.952513983811457, 1e-13);
    EXPECT_NEAR(bp.get(dTTT), 1.6390247075907936, 1e-14);

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
    // print(lambdify([S, K, v, T], diff(result, (S, S)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (S, K)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (S, v)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (S, T)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (K, K)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (K, v)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (K, T)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (v, v)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (v, T)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (T, T)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (S, S, S)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (S, S, K)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (S, S, v)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (S, S, T)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (S, K, K)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (S, K, v)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (S, K, T)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (S, v, v)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (S, v, T)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (S, T, T)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (K, K, K)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (K, K, v)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (K, K, T)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (K, v, v)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (K, v, T)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (K, T, T)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (v, v, v)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (v, v, T)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (v, T, T)))(S0, K0, v0, T0))
    // print(lambdify([S, K, v, T], diff(result, (T, T, T)))(S0, K0, v0, T0))

    TEST_END;
}
