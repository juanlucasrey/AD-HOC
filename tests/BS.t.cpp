#include <adhoc.hpp>
#include <constants_constexpr.hpp>
#include <constants_type.hpp>
#include <erf_inv.hpp>
#include <evaluate_bwd.hpp>
#include <evaluate_fwd.hpp>
#include <tape_static.hpp>

#include "type_name.hpp"

#include <gtest/gtest.h>

namespace adhoc {

template <class D> inline auto cdf_n(D const &x) {
    using constants::CD;
    using constants::encode;
    using std::erfc;
    constexpr double minus_one_over_root_two =
        -1.0 / constexpression::sqrt(2.0);
    return CD<encode(0.5)>() * erfc(x * CD<encode(minus_one_over_root_two)>());
    // return CD<0.5>() * erfc(x * CD<minus_one_over_root_two>());
}

template <class I1, class I2, class I3, class I4>
auto call_price(const I1 &S, const I2 &K, const I3 &v, const I4 &T) {
    using constants::CD;
    using constants::encode;
    using std::log;
    using std::sqrt;
    auto totalvol = v * sqrt(T);
    auto d1 = log(S / K) / totalvol + totalvol * CD<encode(0.5)>();
    auto d2 = d1 + totalvol;
    return S * cdf_n(d1) - K * cdf_n(d2);
}

TEST(BlackScholes, erfcinverse) {
    double x = 0.5;
    double erfc_x = std::erfc(x);
    double x_2 = erfc_inv<true, algorithm::acklam>(erfc_x);
    EXPECT_NEAR(x, x_2, 1e-16);

    double erf_x = std::erf(x);
    double x_3 = erf_inv<true, algorithm::acklam>(erf_x);
    EXPECT_NEAR(x, x_3, 1e-16);

    std::vector<double> original = {0.1, 0.2, 0.3, 0.5, 0.75};
    std::vector<double> ins = {0.1, 0.2, 0.3, 0.5, 0.75};
    for (auto &val : ins) {
        val = std::erf(val);
    }

    erf_inv<algorithm::acklam>(ins.begin(), ins.end());
    for (std::size_t i = 0; i < ins.size(); i++) {
        EXPECT_NEAR(original[i], ins[i], 1e-9);
    }

    std::vector<double> ins2 = {0.1, 0.2, 0.3, 0.5, 0.75};
    for (auto &val : ins2) {
        val = std::erfc(val);
    }

    erfc_inv<algorithm::acklam>(ins2.begin(), ins2.end());
    for (std::size_t i = 0; i < ins2.size(); i++) {
        EXPECT_NEAR(original[i], ins2[i], 1e-9);
    }
}

TEST(BlackScholes, Inverse) {
    double S = 100.0;
    double K = 102.0;
    double v = 0.15;
    double T = 0.5;

    auto totalvol = v * std::sqrt(T);
    auto d1 = std::log(S / K) / totalvol + totalvol * 0.5;
    auto d2 = d1 + totalvol;
    double cdfd1 = cdf_n(d1);
    double cdfd2 = cdf_n(d2);
    double r1 = S * cdfd1;
    double r2 = K * cdfd2;
    double result = r1 - r2;

    double r1_inv = result + r2;
    double r2_inv = r1 - result;
    double S_inv = r1_inv / cdfd1;
    double cdfd1_inv = r1_inv / S;
    double K_inv = r2_inv / cdfd2;
    double cdfd2_inv = r2_inv / K;

    constexpr double minus_root_two = -constexpression::sqrt(2.0);
    double cdfd1_2_inv = cdfd1_inv * 2.;
    // double cdfd1_3_inv =
    //     quantile</* standard_normal, */ false,
    //     algorithm::acklam>(cdfd1_2_inv);
    // double cdfd1_3_inv = erf_inv(cdfd1_2_inv);
    // double d1_inv = cdfd1_3_inv * minus_root_two;

    constexpr double minus_one_over_root_two =
        -1.0 / constexpression::sqrt(2.0);
    double d1_2 = d1 * minus_one_over_root_two;
    double d1_3 = std::erfc(d1_2);
    double cdfd1_2 = 0.5 * d1_3;

    double cdfd2_2_inv = cdfd2_inv * 2.;
    // double cdfd2_3_inv =
    //     quantile</* standard_normal, */ false,
    //     algorithm::acklam>(cdfd2_2_inv);
    // double d2_inv = cdfd2_3_inv * minus_root_two;
}

TEST(BlackScholes, SinglePrice) {
    double result = 0.;

    {
        double S = 100.0;
        double K = 102.0;
        double v = 0.15;
        double T = 0.5;
        result = call_price(S, K, v, T);
    }

    auto tape = CreateTapeInitial<4>();
    auto [S, K, v, T] = tape.getalias();
    tape.set(S, 100.0);
    tape.set(K, 102.0);
    tape.set(v, 0.15);
    tape.set(T, 0.5);
    auto result_adhoc = call_price(S, K, v, T);

    auto intermediate_tape = CreateTapeIntermediate(result_adhoc);
    evaluate_fwd(tape, intermediate_tape);

    double result2 = intermediate_tape.get(result_adhoc);
    EXPECT_EQ(result2, result);

    auto derivatives =
        evaluate_bwd(tape, intermediate_tape, result_adhoc, 1.0, S, K, v, T);
}

} // namespace adhoc