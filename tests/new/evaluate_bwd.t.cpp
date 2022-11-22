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
    // std::cout << tape_d.get(v0) << std::endl;
    // std::cout << -std::sin(0.5) * std::cos(2.) << std::endl;
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
    // std::cout << tape_d.get(v1) << std::endl;
    // std::cout << std::cos(0.5) * -std::sin(2.) << std::endl;
}

} // namespace adhoc2