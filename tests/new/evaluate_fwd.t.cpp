#include <new/evaluate_fwd.hpp>
#include <new/init.hpp>

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

} // namespace adhoc2