#include <new/evaluate_fwd.hpp>

#include <gtest/gtest.h>

namespace adhoc {

TEST(adhoc2, EvaluateFwdUni) {
    auto tape = CreateTapeInitial<1>();
    auto [val0] = tape.getalias();
    auto temp = exp(val0);

    tape.set(val0, 1.0);
    auto intermediate_tape = CreateTapeIntermediate(temp);
    evaluate_fwd(tape, intermediate_tape);

    double result = intermediate_tape.get(temp);
    EXPECT_EQ(result, std::exp(1.0));
}

TEST(adhoc2, EvaluateFwdBi) {
    auto tape = CreateTapeInitial<2>();
    auto [val0, val1] = tape.getalias();

    auto temp = val0 * val1;

    tape.set(val0, 2.0);
    tape.set(val1, 3.0);
    auto intermediate_tape = CreateTapeIntermediate(temp);
    evaluate_fwd(tape, intermediate_tape);

    double result = intermediate_tape.get(temp);
    EXPECT_EQ(result, 6.0);
}

TEST(adhoc2, EvaluateFwdConst) {
    auto tape = CreateTapeInitial<1>();
    auto [val0] = tape.getalias();
    using namespace constants;
    auto temp = val0 * CD<encode(0.5)>();

    tape.set(val0, 1.0);
    auto intermediate_tape = CreateTapeIntermediate(temp);
    evaluate_fwd(tape, intermediate_tape);

    double result = intermediate_tape.get(temp);
    EXPECT_EQ(result, 0.5);
}

template <class D> inline auto cdf_n(D const &x) {
    using std::erfc;
    using namespace constants;
    constexpr double minus_one_over_root_two =
        -1.0 / adhoc::constexpression::sqrt(2.0);
    return CD<encode(0.5)>() * erfc(x * CD<encode(minus_one_over_root_two)>());
    // return CD<0.5>() * erfc(x * CD<minus_one_over_root_two>());
}

template <class I1, class I2, class I3, class I4>
auto call_price(const I1 &S, const I2 &K, const I3 &v, const I4 &T) {
    using std::log;
    using std::sqrt;
    using namespace constants;
    auto totalvol = v * sqrt(T);
    auto d1 = log(S / K) / totalvol + totalvol * CD<encode(0.5)>();
    auto d2 = d1 + totalvol;
    return S * cdf_n(d1) - K * cdf_n(d2);
}

TEST(adhoc2, BSAdhoc) {
    double result = 0;

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
}

} // namespace adhoc