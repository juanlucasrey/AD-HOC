#include <adhoc.hpp>
#include <backpropagator.hpp>
#include <calc_tree.hpp>
#include <differential_operator.hpp>

#include <utils/tuple.hpp>

#include <sort.hpp>

#include "call_price.hpp"

#include "type_name.hpp"

#include <gtest/gtest.h>

#include <random>
#include <utility>

namespace adhoc4 {

TEST(BackPropagator, UnivariateLogExp) {
    ADHOC(x);
    auto res = log(exp(x));

    CalcTree ct(res);

    const std::size_t evals = 100;
    std::mt19937 gen;
    std::uniform_real_distribution<> dis(-3., 3.);

    for (std::size_t i = 0; i < evals; i++) {
        ct.set(x) = dis(gen);
        ct.evaluate();

        BackPropagator bp(d(res), d(x), d<2>(x), d<3>(x), d<4>(x), d<5>(x));
        bp.set(d(res)) = 1.;
        bp.backpropagate(ct);

        EXPECT_NEAR(bp.get(d(x)), 1., 1e-14);
        EXPECT_NEAR(bp.get(d<2>(x)), 0., 1e-14);
        EXPECT_NEAR(bp.get(d<3>(x)), 0., 1e-14);
        EXPECT_NEAR(bp.get(d<4>(x)), 0., 1e-14);
        EXPECT_NEAR(bp.get(d<5>(x)), 0., 1e-13);
    }
}

TEST(BackPropagator, UnivariateExpLog) {
    ADHOC(x);
    auto res = exp(log(x));

    CalcTree ct(res);

    const std::size_t evals = 100;
    std::mt19937 gen;
    std::uniform_real_distribution<> dis(0.1, 3.);

    for (std::size_t i = 0; i < evals; i++) {
        ct.set(x) = dis(gen);
        ct.evaluate();

        BackPropagator bp(d(res), d(x), d<2>(x), d<3>(x), d<4>(x), d<5>(x));
        bp.set(d(res)) = 1.;
        bp.backpropagate(ct);

        EXPECT_NEAR(bp.get(d(x)), 1., 1e-14);
        EXPECT_NEAR(bp.get(d<2>(x)), 0., 1e-14);
        EXPECT_NEAR(bp.get(d<3>(x)), 0., 1e-13);
        EXPECT_NEAR(bp.get(d<4>(x)), 0., 1e-12);
        EXPECT_NEAR(bp.get(d<5>(x)), 0., 1e-10);
    }
}

TEST(BackPropagator2, UnivariateSingle) {
    ADHOC(x);
    auto res = erfc(x);

    CalcTree ct(res);
    ct.set(x) = 1.2;
    ct.evaluate();

    BackPropagator bp(d(res), d(x), d<2>(x), d<3>(x));

    constexpr auto buffer_size = bp.buffer_size<decltype(ct)>();
    static_assert(buffer_size == 0U);
    // std::cout << buffer_size << std::endl;

    bp.set(d(res)) = 1.;
    bp.backpropagate(ct);

    // from sympy import *
    // x = Symbol('x')
    // f = erfc(x)

    // def taylor_coeff(f,x,x0, n):
    //     fdiff = f
    //     for i in range(0,n):
    //         fdiff = fdiff.diff(x)

    //     val = lambdify([x], fdiff)(x0)
    //     val = val / factorial(n)
    //     return val

    // print(taylor_coeff(f, x, 1.2, 1))
    // print(taylor_coeff(f, x, 1.2, 2))
    // print(taylor_coeff(f, x, 1.2, 3))

    EXPECT_NEAR(bp.get(d(x)), -0.267344347003539, 1e-14);
    EXPECT_NEAR(bp.get(d<2>(x)), 0.64162643280849396, 1e-14);
    EXPECT_NEAR(bp.get(d<3>(x)), -1.0052147447333071, 1e-14);
}

TEST(BackPropagator2, UnivariateDouble) {
    ADHOC(x);
    auto res = log(erfc(x));

    CalcTree ct(res);
    ct.set(x) = 1.2;
    ct.evaluate();

    BackPropagator bp(d(res), d(x), d<2>(x), d<3>(x));

    constexpr auto buffer_size = bp.buffer_size<decltype(ct)>();
    static_assert(buffer_size == 3U);
    // std::cout << buffer_size << std::endl;

    bp.set(d(res)) = 1.;
    bp.backpropagate(ct);

    // from sympy import *
    // x = Symbol('x')
    // f = log(erfc(x))

    // def taylor_coeff(f,x,x0, n):
    //     fdiff = f
    //     for i in range(0,n):
    //         fdiff = fdiff.diff(x)

    //     val = lambdify([x], fdiff)(x0)
    //     val = val / factorial(n)
    //     return val

    // print(taylor_coeff(f, x, 1.2, 1))
    // print(taylor_coeff(f, x, 1.2, 2))
    // print(taylor_coeff(f, x, 1.2, 3))

    EXPECT_NEAR(bp.get(d(x)), -2.98089202449025, 1e-14);
    EXPECT_NEAR(bp.get(d<2>(x)), -1.731576402892995, 1e-14);
    EXPECT_NEAR(bp.get(d<3>(x)), -0.20571716243479798, 1e-14);
}

TEST(BackPropagator2, SumSingle) {
    ADHOC(x);
    ADHOC(y);
    // auto res = x + y;
    auto res = log(x + y);

    CalcTree ct(res);
    ct.set(x) = 1.2;
    ct.set(y) = 0.4;
    ct.evaluate();

    BackPropagator bp(d(res), d(x), d(y), d(x) * d(y), d<5>(x));

    constexpr auto buffer_size = bp.buffer_size<decltype(ct)>();
    static_assert(buffer_size == 3U);
    // std::cout << buffer_size << std::endl;

    bp.set(d(res)) = 1.;
    bp.backpropagate(ct);

    EXPECT_NEAR(bp.get(d(x)), 0.625, 1e-13);
    EXPECT_NEAR(bp.get(d(y)), 0.625, 1e-13);
    EXPECT_NEAR(bp.get(d(x) * d(y)), -0.390625, 1e-13);
}

TEST(BackPropagator2, MulSingle) {
    ADHOC(x);
    ADHOC(y);
    auto res = log(x * y);

    CalcTree ct(res);
    ct.set(x) = 1.2;
    ct.set(y) = 0.4;
    ct.evaluate();

    BackPropagator bp(d(res), d(x), d(y), d(x) * d(y), d<5>(x));
    bp.set(d(res)) = 1.;
    bp.backpropagate(ct);

    EXPECT_NEAR(bp.get(d(x)), 0.833333333333333, 1e-14);
    EXPECT_NEAR(bp.get(d(y)), 2.5, 1e-14);
    EXPECT_NEAR(bp.get(d(x) * d(y)), 0., 1e-14);
    EXPECT_NEAR(bp.get(d<5>(x)), 9.645061728395067, 1e-14);
}

TEST(BackPropagator2, MulSingle2) {
    ADHOC(x);
    ADHOC(y);
    auto res = erfc(x * y);

    CalcTree ct(res);
    ct.set(x) = 1.2;
    ct.set(y) = 0.4;
    ct.evaluate();

    BackPropagator bp(d(res), d(x), d(y), d(x) * d(y), d<5>(x));

    constexpr auto buffer_size = bp.buffer_size<decltype(ct)>();
    static_assert(buffer_size == 3U);
    // std::cout << buffer_size << std::endl;

    bp.set(d(res)) = 1.;
    bp.backpropagate(ct);

    // std::cout.precision(std::numeric_limits<double>::max_digits10);
    // std::cout << bp.get(d(x)) << std::endl;
    // std::cout << bp.get(d(y)) << std::endl;
    // std::cout << bp.get(d(x) * d(y)) << std::endl;
    // std::cout << bp.get(d<5>(x)) << std::endl;

    EXPECT_NEAR(bp.get(d(x)), -0.358470648907805, 1e-14);
    EXPECT_NEAR(bp.get(d(y)), -1.07541194672341, 1e-14);
    EXPECT_NEAR(bp.get(d(x) * d(y)), -0.483218434727721, 1e-14);
    EXPECT_NEAR(bp.get(d<5>(x)), -0.016427903974483831, 1e-14);
}

TEST(BackPropagator2, DivSingle) {
    ADHOC(x);
    ADHOC(y);
    auto res = erfc(x / y);

    CalcTree ct(res);
    ct.set(x) = 1.2;
    ct.set(y) = 0.4;
    ct.evaluate();

    BackPropagator bp(d(res), d(x), d(y), d(x) * d(y), d<5>(x));
    bp.set(d(res)) = 1.;
    bp.backpropagate(ct);

    EXPECT_NEAR(bp.get(d(x)), -0.000348132629866870, 1e-14);
    EXPECT_NEAR(bp.get(d(y)), 0.00104439788960061, 1e-14);
    EXPECT_NEAR(bp.get(d(x) * d(y)), -0.0147956367693420, 1e-14);
    EXPECT_NEAR(bp.get(d<5>(x)), -11.912663428256959, 1e-14);
}

TEST(BackPropagator2, MulConstant) {
    ADHOC(x);
    using constants::CD;
    auto res = log(x * CD<0.4>());

    CalcTree ct(res);
    ct.set(x) = 1.2;
    ct.evaluate();

    BackPropagator bp(d(res), d(x), d<5>(x));
    bp.set(d(res)) = 1.;
    bp.backpropagate(ct);

    EXPECT_NEAR(bp.get(d(x)), 0.833333333333333, 1e-14);
    // EXPECT_NEAR(bp.get(d(y)), 2.5, 1e-14);
    // EXPECT_NEAR(bp.get(d(x) * d(y)), 0., 1e-14);
    EXPECT_NEAR(bp.get(d<5>(x)), 9.645061728395067, 1e-14);
}

namespace detail {

template <class I1, class I2>
constexpr auto equal(I1 /* lhs */, I2 /* rhs */) -> bool {
    return std::is_same_v<I1, I2>;
}

} // namespace detail

TEST(BackPropagator2, MultiplyOrdered) {
    using constants::CD;
    ADHOC(S);
    ADHOC(K);
    ADHOC(v);
    ADHOC(T);

    auto res = call_price(S, K, v, T);

    CalcTree ct(res);
    using NodesValue = decltype(ct)::ValuesTupleInverse;

    {
        constexpr auto dS = d(S);
        constexpr auto diffop1mul =
            detail::multiply_ordered(dS, dS, NodesValue{});

        constexpr auto result = d<2>(S);
        static_assert(detail::equal(diffop1mul, result));
    }

    {
        constexpr auto dS = d(S);
        constexpr auto dK = d(K);
        constexpr auto diffop1mul =
            detail::multiply_ordered(dS, dK, NodesValue{});

        constexpr auto result = d(S) * d(K);
        static_assert(detail::equal(diffop1mul, result));
    }

    {
        constexpr auto dS = d(S);
        constexpr auto dK = d(K) * d(T);
        constexpr auto diffop1mul =
            detail::multiply_ordered(dS, dK, NodesValue{});

        constexpr auto result = d(S) * d(K) * d(T);
        static_assert(detail::equal(diffop1mul, result));
    }

    {
        constexpr auto dS = d(S) * d(v);
        constexpr auto dK = d(K) * d(T);
        constexpr auto diffop1mul =
            detail::multiply_ordered(dS, dK, NodesValue{});

        constexpr auto result = d(S) * d(K) * d(v) * d(T);
        static_assert(detail::equal(diffop1mul, result));
    }

    {
        constexpr auto dS = d(S) * d<2>(K) * d(v);
        constexpr auto dK = d(K) * d(T);
        constexpr auto diffop1mul =
            detail::multiply_ordered(dS, dK, NodesValue{});

        constexpr auto result = d(S) * d<3>(K) * d(v) * d(T);
        static_assert(detail::equal(diffop1mul, result));
    }
}

TEST(BackPropagator, IsInputDer) {
    ADHOC(var);
    ADHOC(var2);
    auto calc = exp(var);
    auto res = d(var) * d<2>(var2);

    constexpr auto resb1 = detail::is_derivative_input(res);
    static_assert(resb1);

    auto res2 = d(calc) * d<2>(var2);

    constexpr auto resb2 = detail::is_derivative_input(res2);
    static_assert(!resb2);
}

TEST(BackPropagator, SecondDerCheck) {

    double x0 = 1.02;
    ADHOC(x);
    auto res = exp(x) + x * exp(x);
    CalcTree ct(res);
    ct.set(x) = x0;
    ct.evaluate();

    auto dXX = d<2>(x);
    auto dres = d(res);
    BackPropagator t(dXX, dres);

    t.set(dres) = 1.;
    t.backpropagate(ct);

    // from sympy import *
    // x = Symbol('x')
    // result = exp(x) + x * exp(x)
    // def diff_function(f, diffs):
    //     fdiff = f
    //     for i in diffs:
    //         fdiff = fdiff.diff(i)
    //     return fdiff
    // x0 = 1.02
    // print(lambdify([x], diff_function(result, (x, x)))(x0))

    EXPECT_NEAR(t.get(dXX), 11.148242951136478, 1e-14);
}

} // namespace adhoc4
