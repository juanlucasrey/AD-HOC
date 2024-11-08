#include <calc_tree.hpp>

#include <adhoc.hpp>
#include <dependency.hpp>

#include "call_price.hpp"
#include "test_tools.hpp"

int main() {
    TEST_START;
    using namespace adhoc4;

    // EvaluateFwdUni2
    {
        ADHOC(val0);
        auto temp = exp(val0);
        CalcTree t(temp);

        t.set(val0) = 1.0;
        t.evaluate();

        double const result = t.get(temp);
        EXPECT_EQUAL(result, std::exp(1.0));
    }

    // EvaluateFwdBi2
    {
        ADHOC(val0);
        ADHOC(val1);
        auto temp = val0 * val1;

        CalcTree t(temp);

        t.set(val0) = 2.0;
        t.set(val1) = 3.0;
        t.evaluate();

        double const result = t.get(temp);
        EXPECT_EQUAL(result, 6.0);
    }

    // EvaluateFwdConst2
    {
        ADHOC(val0);
        using constants::CD;
        auto temp = val0 * CD<0.5>();

        CalcTree t(temp);
        t.set(val0) = 1.0;
        t.evaluate();

        double const result = t.get(temp);
        EXPECT_EQUAL(result, 0.5);
    }

    // EvaluateFwdTwoVals2
    {
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
        EXPECT_NEAR_ABS(result, 2., 1e-10);
        EXPECT_NEAR_ABS(result2, 1.0806046117362795, 1e-10);
    }

    // EvaluateFwdDivision2
    {
        ADHOC(val0);
        ADHOC(val1);
        auto temp = val0 / val1;

        CalcTree t(temp);
        t.set(val0) = 1.0;
        t.set(val1) = 2.0;
        t.evaluate();

        double const result = t.get(temp);
        EXPECT_NEAR_ABS(result, 1.0 / 2.0, 1e-10);
    }

    // CallPrice
    {
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
        // inlining and fma CAN create differences
        EXPECT_NEAR_ABS(result_calctree, result_double, 1e-14);
    }
    TEST_END;
}
