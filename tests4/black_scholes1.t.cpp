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

TEST(BlackScholes, BS_FD) {
    double S = 1.01;
    auto res = call_price2(S, 1.02, 0.15, 0.5);

    double epsilon = 1e-8;
    double Sup = S + epsilon;
    auto resup = call_price2(Sup, 1.02, 0.15, 0.5);

    auto dS = (resup - res) / epsilon;

    std::cout.precision(std::numeric_limits<double>::max_digits10);

    std::cout << res << std::endl;
    std::cout << dS << std::endl;
}

TEST(BlackScholes, BS1) {
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

    auto dS = d(S);

    auto dres = d(res);

    BackPropagator t(dS, dres);

    t.set(dres) = 1.;
    t.backpropagate2(ct);

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    std::cout << ct.get(res) << std::endl;
    std::cout << t.get(dS) << std::endl;
}

} // namespace adhoc4
