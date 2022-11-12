#include <adhoc.hpp>
#include <constants_constexpr.hpp>
#include <constants_type.hpp>
#include <evaluate_bwd.hpp>
#include <evaluate_fwd.hpp>
#include <tape_static.hpp>

#include "../type_name.hpp"

#include <gtest/gtest.h>

namespace adhoc {

template <class D> inline auto cdf_n(D const &x) {
    using std::erfc;
    using namespace constants;
    constexpr double minus_one_over_root_two =
        -1.0 / constexpression::sqrt(2.0);
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

TEST(BlackScholes, SinglePrice) {
    double result;

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