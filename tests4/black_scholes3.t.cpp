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

TEST(BlackScholes, BS3) {
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

    // order 3
    auto dSSS = d<3>(S);
    auto dKKK = d<3>(K);
    auto dvvv = d<3>(v);
    auto dTTT = d<3>(T);

    auto dSSK = d<2>(S) * d(K);
    auto dSSv = d<2>(S) * d(v);
    auto dSST = d<2>(S) * d(T);
    auto dSKK = d(S) * d<2>(K);
    auto dSvv = d(S) * d<2>(v);
    auto dSTT = d(S) * d<2>(T);
    auto dKKv = d<2>(K) * d(v);
    auto dKKT = d<2>(K) * d(T);
    auto dKvv = d(K) * d<2>(v);
    auto dKTT = d(K) * d<2>(T);
    auto dvvT = d<2>(v) * d(T);
    auto dvTT = d(v) * d<2>(T);

    auto dSKv = d(S) * d(K) * d(v);
    auto dSKT = d(S) * d(K) * d(T);
    auto dSvT = d(S) * d(v) * d(T);
    auto dKvT = d(K) * d(v) * d(T);

    auto dres = d(res);

    BackPropagator t(dS, dK, dv, dT, dSS, dKK, dvv, dTT, dSK, dSv, dST, dKv,
                     dKT, dvT, dSSS, dKKK, dvvv, dTTT, dSSK, dSSv, dSST, dSKK,
                     dSvv, dSTT, dKKv, dKKT, dKvv, dKTT, dvvT, dvTT, dSKv, dSKT,
                     dSvT, dKvT, dres);

    t.set(dres) = 1.;
    t.backpropagate(ct);

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    std::cout << ct.get(res) << std::endl;
    std::cout << t.get(dS) << std::endl;
}

} // namespace adhoc4
