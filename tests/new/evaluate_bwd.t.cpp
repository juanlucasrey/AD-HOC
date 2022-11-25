#include <new/evaluate_bwd.hpp>
#include <new/evaluate_fwd.hpp>
#include <new/init.hpp>

#include <gtest/gtest.h>

namespace adhoc2 {

TEST(EvaluateBwd, Univariate) {
    auto [v0] = Init<1>();
    // auto r = exp(cos(v0) * (v1 * v1));
    auto r = exp(cos(log(v0)));

    auto leaves_and_roots = TapeRootsAndLeafs(r);
    leaves_and_roots.set(v0, 0.5);

    auto intermediate_tape = evaluate_fwd_return_vals(leaves_and_roots);
    auto tape_d = TapeDerivatives(v0, r);
    tape_d.get(r) = 1.0;

    evaluate_bwd(leaves_and_roots, intermediate_tape, tape_d);
    EXPECT_NEAR(tape_d.get(v0), 2.757914160416556, 1e-13);
}

TEST(EvaluateBwd, Univariate2) {
    auto [v0] = Init<1>();
    // auto r = exp(cos(v0) * (v1 * v1));
    auto r = exp(v0);

    auto leaves_and_roots = TapeRootsAndLeafs(r);
    leaves_and_roots.set(v0, 0.5);

    auto intermediate_tape = evaluate_fwd_return_vals(leaves_and_roots);
    auto tape_d = TapeDerivatives(v0, r);
    tape_d.get(r) = 1.0;

    evaluate_bwd(leaves_and_roots, intermediate_tape, tape_d);
    EXPECT_NEAR(tape_d.get(v0), 1.6487212707001282, 1e-13);
}

TEST(EvaluateBwd, BivariateCutLeaf1) {
    auto [v0, v1] = Init<2>();
    // auto r = exp(cos(v0) * (v1 * v1));
    auto r = cos(v0) * cos(v1);

    auto leaves_and_roots = TapeRootsAndLeafs(r);
    leaves_and_roots.set(v0, 0.5);
    leaves_and_roots.set(v1, 2.);

    auto intermediate_tape = evaluate_fwd_return_vals(leaves_and_roots);
    auto tape_d = TapeDerivatives(v0, r);
    tape_d.get(r) = 1.0;

    evaluate_bwd(leaves_and_roots, intermediate_tape, tape_d);
    EXPECT_NEAR(tape_d.get(v0), -std::sin(0.5) * std::cos(2.), 1e-13);
}

TEST(EvaluateBwd, BivariateCutLeaf2) {
    auto [v0, v1] = Init<2>();
    // auto r = exp(cos(v0) * (v1 * v1));
    auto r = cos(v0) * cos(v1);

    auto leaves_and_roots = TapeRootsAndLeafs(r);
    leaves_and_roots.set(v0, 0.5);
    leaves_and_roots.set(v1, 2.);

    auto intermediate_tape = evaluate_fwd_return_vals(leaves_and_roots);
    auto tape_d = TapeDerivatives(v1, r);
    tape_d.get(r) = 1.0;

    evaluate_bwd(leaves_and_roots, intermediate_tape, tape_d);
    EXPECT_NEAR(tape_d.get(v1), std::cos(0.5) * -std::sin(2.), 1e-13);
}

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

TEST(EvaluateBwd, BlackScholes) {
    double result = 0.;

    {
        double S = 100.0;
        double K = 102.0;
        double v = 0.15;
        double T = 0.5;
        result = call_price(S, K, v, T);
    }

    auto [S, K, v, T] = Init<4>();
    auto result_adhoc = call_price(S, K, v, T);

    auto leaves_and_roots = TapeRootsAndLeafs(result_adhoc);
    leaves_and_roots.set(S, 100.0);
    leaves_and_roots.set(K, 102.0);
    leaves_and_roots.set(v, 0.15);
    leaves_and_roots.set(T, 0.5);

    auto intermediate_tape = evaluate_fwd_return_vals(leaves_and_roots);
    auto tape_d = TapeDerivatives(S, K, v, T, result_adhoc);

    tape_d.get(result_adhoc) = 1.0;

    evaluate_bwd(leaves_and_roots, intermediate_tape, tape_d);

    std::cout << tape_d.get(S) << std::endl;
    std::cout << tape_d.get(K) << std::endl;
    std::cout << tape_d.get(v) << std::endl;
    std::cout << tape_d.get(T) << std::endl;
}

} // namespace adhoc2