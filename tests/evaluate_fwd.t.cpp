#include <constants_type.hpp>
#include <evaluate_fwd.hpp>
#include <init.hpp>
#include <tape.hpp>

#include "call_price.hpp"

#include <gtest/gtest.h>

namespace adhoc {

TEST(EvaluateFwd, EvaluateFwdUni) {
    auto [val0] = Init<1>();
    auto temp = exp(val0);

    auto t = Tape3(temp);

    t.val(val0) = 1.0;
    evaluate_fwd(t);

    double result = t.val(temp);
    EXPECT_EQ(result, std::exp(1.0));
}

TEST(EvaluateFwd, EvaluateFwdBi) {
    auto [val0, val1] = Init<2>();
    auto temp = val0 * val1;

    auto t = Tape3(temp);

    t.val(val0) = 2.0;
    t.val(val1) = 3.0;
    evaluate_fwd(t);

    double result = t.val(temp);
    EXPECT_EQ(result, 6.0);
}

TEST(EvaluateFwd, EvaluateFwdConst) {
    auto [val0] = Init<1>();
    using namespace constants;
    auto temp = val0 * CD<encode(0.5)>();

    auto t = Tape3(temp);
    t.val(val0) = 1.0;
    evaluate_fwd(t);

    double result = t.val(temp);
    EXPECT_EQ(result, 0.5);
}

TEST(EvaluateFwd, BSAdhoc) {
    double result = 0;

    {
        double S = 100.0;
        double K = 102.0;
        double v = 0.15;
        double T = 0.5;
        result = call_price(S, K, v, T);
    }

    auto [S, K, v, T] = Init<4>();
    auto result_adhoc = call_price(S, K, v, T);

    auto t = Tape3(result_adhoc);

    t.val(S) = 100.0;
    t.val(K) = 102.0;
    t.val(v) = 0.15;
    t.val(T) = 0.5;
    evaluate_fwd(t);

    double result2 = t.val(result_adhoc);
    EXPECT_EQ(result2, result);
}

} // namespace adhoc