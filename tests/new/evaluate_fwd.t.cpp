#include <new/constants_type.hpp>
#include <new/evaluate_fwd.hpp>
#include <new/init.hpp>

#include "call_price.hpp"

#include <gtest/gtest.h>

namespace adhoc2 {

TEST(EvaluateFwd, EvaluateFwdUni) {
    auto [val0] = Init<1>();
    auto temp = exp(val0);

    auto leaves_and_roots = TapeRootsAndLeafs(temp);

    leaves_and_roots.get(val0) = 1.0;
    auto intermediate_tape = evaluate_fwd_return_vals(leaves_and_roots);

    double result = leaves_and_roots.get(temp);
    EXPECT_EQ(result, std::exp(1.0));
}

TEST(EvaluateFwd, EvaluateFwdBi) {
    auto [val0, val1] = Init<2>();
    auto temp = val0 * val1;

    auto leaves_and_roots = TapeRootsAndLeafs(temp);

    leaves_and_roots.get(val0) = 2.0;
    leaves_and_roots.get(val1) = 3.0;
    auto intermediate_tape = evaluate_fwd_return_vals(leaves_and_roots);

    double result = leaves_and_roots.get(temp);
    EXPECT_EQ(result, 6.0);
}

TEST(EvaluateFwd, EvaluateFwdConst) {
    auto [val0] = Init<1>();
    using namespace constants;
    auto temp = val0 * CD<encode(0.5)>();

    auto leaves_and_roots = TapeRootsAndLeafs(temp);
    leaves_and_roots.get(val0) = 1.0;
    auto intermediate_tape = evaluate_fwd_return_vals(leaves_and_roots);

    double result = leaves_and_roots.get(temp);
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

    auto leaves_and_roots = TapeRootsAndLeafs(result_adhoc);

    leaves_and_roots.get(S) = 100.0;
    leaves_and_roots.get(K) = 102.0;
    leaves_and_roots.get(v) = 0.15;
    leaves_and_roots.get(T) = 0.5;
    auto intermediate_tape = evaluate_fwd_return_vals(leaves_and_roots);

    double result2 = leaves_and_roots.get(result_adhoc);
    EXPECT_EQ(result2, result);
}

} // namespace adhoc2