#include <calc_tree.hpp>

#include <adhoc.hpp>
#include <dependency.hpp>

#include "call_price.hpp"
#include "type_name.hpp"

#include <gtest/gtest.h>

namespace adhoc4 {

TEST(CalcTree, IsInput) {
    ADHOC(vol);
    ADHOC(strike);
    ADHOC(spot);

    auto res = vol * strike + exp(spot);
    static_assert(is_input(vol));
    static_assert(is_input(strike));
    static_assert(is_input(spot));
    static_assert(!is_input(res));

    constexpr auto tree = detail::calc_tree_t(res);
    std::cout << type_name2<decltype(tree)>() << std::endl;

    constexpr auto resb1 = (vol == strike);
    static_assert(!resb1);

    constexpr auto resb2 = (vol == vol);
    static_assert(resb2);

    constexpr auto resb3 = (res == res);
    static_assert(resb3);
}

TEST(CalcTree, ConstType) {
    ADHOC(val0);

    using constants::CD;
    auto res = val0 * CD<0.5>();

    constexpr auto tree = detail::calc_tree_t(res);
    std::cout << type_name2<decltype(tree)>() << std::endl;

    // constexpr auto resb1 = (vol == strike);
    // static_assert(!resb1);

    // constexpr auto resb2 = (vol == vol);
    // static_assert(resb2);

    // constexpr auto resb3 = (res == res);
    // static_assert(resb3);
}

TEST(CalcTree, Repeat) {
    ADHOC(var1);
    ADHOC(var2);
    ADHOC(var3);
    ADHOC(var4);

    auto res1 = (var1 + var2) * (var3 * var4);
    auto res2 = (var4 * var3) * (var2 + var1);
    // auto res2 = (var3 * var4) * (var1 + var2);
    constexpr auto resb1 = (res1 == res2);
    static_assert(resb1);

    auto res3 = (var1 * var2) + (var3 + var4);
    auto res4 = (var4 + var3) + (var2 * var1);
    constexpr auto resb2 = (res3 == res4);
    static_assert(resb2);
}

TEST(CalcTree, EvaluateFwdUni2) {
    ADHOC(val0);
    ADHOC(val1);
    auto temp = exp(val0);
    CalcTree t(temp);

    t.set(val0) = 1.0;
    // t.set(val1) = 1.0; // fails because it's not on the calculation tree
    // t.set(temp) = 1.0; // fails because it's not an input
    t.evaluate();

    double const result = t.get(temp);
    EXPECT_EQ(result, std::exp(1.0));
}

TEST(CalcTree, EvaluateFwdBi2) {
    ADHOC(val0);
    ADHOC(val1);
    auto temp = val0 * val1;

    CalcTree t(temp);

    t.set(val0) = 2.0;
    t.set(val1) = 3.0;
    // t.set(temp) = 3.0; // fails because it's not an input
    t.evaluate();

    double const result = t.get(temp);
    EXPECT_EQ(result, 6.0);
}

TEST(CalcTree, EvaluateFwdConst2) {
    ADHOC(val0);
    using constants::CD;
    auto temp = val0 * CD<0.5>();

    CalcTree t(temp);
    t.set(val0) = 1.0;
    // t.set(temp) = 1.0; // fails because it's not an input
    t.evaluate();

    double const result = t.get(temp);
    EXPECT_EQ(result, 0.5);
}

TEST(CalcTree, EvaluateFwdTwoVals2) {
    ADHOC(val0);
    ADHOC(val1);
    auto temp = val0 * val1;
    auto temp2 = cos(val0) * val1;

    CalcTree t(temp, temp2);
    t.set(val0) = 1.0;
    t.set(val1) = 2.0;
    t.evaluate();

    double const result = t.get(temp);
    double const result2 = t.get(temp2);
    EXPECT_NEAR(result, 2., 1e-10);
    EXPECT_NEAR(result2, 1.0806046117362795, 1e-10);
}

TEST(CalcTree, EvaluateFwdDivision2) {
    ADHOC(val0);
    ADHOC(val1);
    auto temp = val0 / val1;

    CalcTree t(temp);
    t.set(val0) = 1.0;
    t.set(val1) = 2.0;
    t.evaluate();

    double const result = t.get(temp);
    EXPECT_NEAR(result, 1.0 / 2.0, 1e-10);
}

TEST(CalcTree, CallPrice) {
    ADHOC(S);
    ADHOC(K);
    ADHOC(v);
    ADHOC(T);

    auto res = call_price(S, K, v, T);

    CalcTree ct(res);
    ct.set(S) = 100.0;
    ct.set(K) = 102.0;
    ct.set(v) = 0.15;
    ct.set(T) = 0.5;
    ct.evaluate();

    double const result_calctree = ct.get(res);

    double const result_double = call_price(100.0, 102.0, 0.15, 0.5);
    EXPECT_EQ(result_calctree, result_double);
}

} // namespace adhoc4
