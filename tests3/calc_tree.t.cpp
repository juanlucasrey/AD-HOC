#include <calc_tree.hpp>
#include <init.hpp>

#include "type_name.hpp"

#include <gtest/gtest.h>

namespace adhoc3 {

TEST(CalcTree, EvaluateFwdUni) {
    auto [val0] = Init<1>();
    auto temp = exp(val0);

    std::cout << type_name2<decltype(temp)>() << std::endl;

    CalcTree t(temp);

    std::cout << type_name2<decltype(val0)>() << std::endl;

    t.set(val0) = 1.0;
    // t.set(temp) = 1.0; // fails because it's not a leaf
    t.evaluate();

    double result = t.val(temp);
    EXPECT_EQ(result, std::exp(1.0));
}

TEST(CalcTree, EvaluateFwdBi) {
    auto [val0, val1] = Init<2>();
    auto temp = val0 * val1;

    CalcTree t(temp);

    t.set(val0) = 2.0;
    t.set(val1) = 3.0;
    // t.set(temp) = 3.0; // fails because it's not a leaf
    t.evaluate();

    double result = t.val(temp);
    EXPECT_EQ(result, 6.0);
}

TEST(CalcTree, EvaluateFwdConst) {
    auto [val0] = Init<1>();
    using constants::CD;
    using constants::encode;
    auto temp = val0 * CD<encode(0.5)>();

    CalcTree t(temp);
    t.set(val0) = 1.0;
    // t.set(temp) = 1.0; // fails because it's not a leaf
    t.evaluate();

    double result = t.val(temp);
    EXPECT_EQ(result, 0.5);
}

TEST(CalcTree, EvaluateFwdTwoVals) {
    auto [val0, val1] = Init<2>();
    auto temp = val0 * val1;
    auto temp2 = cos(val0) * val1;

    CalcTree t(temp, temp2);
    t.set(val0) = 1.0;
    t.set(val1) = 2.0;
    t.evaluate();

    double result = t.val(temp);
    double result2 = t.val(temp2);
    EXPECT_NEAR(result, 2., 1e-10);
    EXPECT_NEAR(result2, 1.0806046117362795, 1e-10);
}

TEST(CalcTree, EvaluateFwdDivision) {
    auto [val0, val1] = Init<2>();
    auto temp = val0 / val1;

    CalcTree t(temp);
    t.set(val0) = 1.0;
    t.set(val1) = 2.0;
    t.evaluate();

    double result = t.val(temp);
    EXPECT_NEAR(result, 1.0 / 2.0, 1e-10);
}

// TEST(CalcTree, BSAdhoc) {
//     double result = 0;

//     {
//         double S = 100.0;
//         double K = 102.0;
//         double v = 0.15;
//         double T = 0.5;
//         result = call_price(S, K, v, T);
//     }

//     auto [S, K, v, T] = Init<4>();
//     auto result_adhoc = call_price(S, K, v, T);

//     CalcTree t(result_adhoc);

//     t.set(S) = 100.0;
//     t.set(K) = 102.0;
//     t.set(v) = 0.15;
//     t.set(T) = 0.5;
//     t.evaluate_fwd();

//     double result2 = t.val(result_adhoc);
//     EXPECT_EQ(result2, result);
// }

// TEST(CalcTree, CallAndPut) {
//     double result_call = 0.;
//     double result_put = 0.;

//     {
//         double S = 100.0;
//         double K = 102.0;
//         double v = 0.15;
//         double T = 0.5;
//         result_call = call_price(S, K, v, T);
//         result_put = put_price(S, K, v, T);
//     }

//     auto [S, K, v, T] = Init<4>();
//     auto rc = call_price(S, K, v, T);
//     auto rp = put_price(S, K, v, T);

//     CalcTree t(rc, rp);

//     t.set(S) = 100.0;
//     t.set(K) = 102.0;
//     t.set(v) = 0.15;
//     t.set(T) = 0.5;
//     t.evaluate_fwd();

//     double result_call2 = t.val(rc);
//     EXPECT_EQ(result_call2, result_call);
//     double result_put2 = t.val(rp);
//     EXPECT_EQ(result_put2, result_put);
// }

} // namespace adhoc3
