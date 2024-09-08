
#include <adhoc.hpp>
#include <constants_type.hpp>
#include <order.hpp>

#include <monomial_included.hpp>

#include <differential_operator.hpp>

#include <gtest/gtest.h>

#include "type_name.hpp"

namespace adhoc4 {

TEST(SaturationArithmetic, Inverse) {

    constexpr auto max = std::numeric_limits<std::size_t>::max();
    constexpr auto min = std::numeric_limits<std::size_t>::min();

    constexpr auto a = detail::add_sat(static_cast<std::size_t>(3),
                                       static_cast<std::size_t>(4));
    static_assert(a == 7);

    constexpr auto b = detail::add_sat(max, static_cast<std::size_t>(4));
    static_assert(b == max);

    constexpr auto c = detail::add_sat(max, max, max);
    static_assert(c == max);

    constexpr auto d = detail::add_sat(static_cast<std::size_t>(3),
                                       static_cast<std::size_t>(3),
                                       static_cast<std::size_t>(3));
    static_assert(d == 9);
}

TEST(Order, Const) {
    using constants::CD;

    ADHOC(var);
    auto constant = CD<-0.5>();

    // static_assert(order<decltype(constant), decltype(var)>::call() == 0);
    static_assert(order(constant, var) == 0);
}

TEST(Order, VarType) {
    ADHOC(var);
    ADHOC(var2);

    // static_assert(order<decltype(var), decltype(var)>::call() == 1);
    // static_assert(order<decltype(var), decltype(var2)>::call() == 0);
    static_assert(order(var, var) == 1);
    static_assert(order(var, var2) == 0);
}

TEST(Order, UnivariateType) {
    ADHOC(var);

    auto var2 = exp(var);
    // static_assert(order<decltype(var2), decltype(var)>::call() ==
    //               std::numeric_limits<std::size_t>::max());
    static_assert(order(var2, var) == std::numeric_limits<std::size_t>::max());
}

TEST(Order, Add) {
    ADHOC(var1);
    ADHOC(var2);

    auto var3 = var1 + var2;
    // static_assert(order<decltype(var3), decltype(var1)>::call() == 1);
    // static_assert(order<decltype(var3), decltype(var2)>::call() == 1);
    static_assert(order(var3, var1) == 1);
    static_assert(order(var3, var1) == 1);
}

TEST(Order, Sub) {
    ADHOC(var1);
    ADHOC(var2);

    auto var3 = var1 - var2;
    // static_assert(order<decltype(var3), decltype(var1)>::call() == 1);
    // static_assert(order<decltype(var3), decltype(var2)>::call() == 1);
    static_assert(order(var3, var1) == 1);
    static_assert(order(var3, var1) == 1);
}

TEST(Order, Mul) {
    ADHOC(var1);
    ADHOC(var2);

    auto var3 = var1 * var2;
    // static_assert(order<decltype(var3), decltype(var1)>::call() == 1);
    // static_assert(order<decltype(var3), decltype(var2)>::call() == 1);
    static_assert(order(var3, var1) == 1);
    static_assert(order(var3, var1) == 1);

    auto var4 = var1 * var1;
    // static_assert(order<decltype(var4), decltype(var1)>::call() == 2);
    static_assert(order(var4, var1) == 2);
}

TEST(Order, Div) {
    ADHOC(var1);
    ADHOC(var2);

    auto var3 = var1 / var2;
    // static_assert(order<decltype(var3), decltype(var1)>::call() == 1);
    // static_assert(order<decltype(var3), decltype(var2)>::call() ==
    //               std::numeric_limits<std::size_t>::max());
    static_assert(order(var3, var1) == 1);
    static_assert(order(var3, var2) == std::numeric_limits<std::size_t>::max());
}

TEST(Order, First) {
    ADHOC(S);
    ADHOC(K);
    ADHOC(v);
    ADHOC(T);

    auto var = exp(T);
    // static_assert(order<decltype(var), decltype(T)>::call() ==
    //               std::numeric_limits<std::size_t>::max());
    static_assert(order(var, T) == std::numeric_limits<std::size_t>::max());

    auto sumv = K + v;
    // static_assert(order<decltype(sumv), decltype(S)>::call() == 0);
    // static_assert(order<decltype(sumv), decltype(K)>::call() == 1);
    // static_assert(order<decltype(sumv), decltype(v)>::call() == 1);
    // static_assert(order<decltype(sumv), decltype(T)>::call() == 0);
    static_assert(order(sumv, S) == 0);
    static_assert(order(sumv, K) == 1);
    static_assert(order(sumv, v) == 1);
    static_assert(order(sumv, T) == 0);

    auto mulv = K * v * K;
    // static_assert(order<decltype(mulv), decltype(S)>::call() == 0);
    // static_assert(order<decltype(mulv), decltype(K)>::call() == 2);
    // static_assert(order<decltype(mulv), decltype(v)>::call() == 1);
    // static_assert(order<decltype(mulv), decltype(T)>::call() == 0);
    static_assert(order(mulv, S) == 0);
    static_assert(order(mulv, K) == 2);
    static_assert(order(mulv, v) == 1);
    static_assert(order(mulv, T) == 0);

    auto mulexpv = K * v * exp(K);
    // static_assert(order<decltype(mulexpv), decltype(S)>::call() == 0);
    // static_assert(order<decltype(mulexpv), decltype(K)>::call() ==
    //               std::numeric_limits<std::size_t>::max());
    // static_assert(order<decltype(mulexpv), decltype(v)>::call() == 1);
    // static_assert(order<decltype(mulexpv), decltype(T)>::call() == 0);
    static_assert(order(mulexpv, S) == 0);
    static_assert(order(mulexpv, K) == std::numeric_limits<std::size_t>::max());
    static_assert(order(mulexpv, v) == 1);
    static_assert(order(mulexpv, T) == 0);

    auto divv = v / K;
    // static_assert(order<decltype(divv), decltype(S)>::call() == 0);
    // static_assert(order<decltype(divv), decltype(K)>::call() ==
    //               std::numeric_limits<std::size_t>::max());
    // static_assert(order<decltype(divv), decltype(v)>::call() == 1);
    // static_assert(order<decltype(divv), decltype(T)>::call() == 0);
    static_assert(order(divv, S) == 0);
    static_assert(order(divv, K) == std::numeric_limits<std::size_t>::max());
    static_assert(order(divv, v) == 1);
    static_assert(order(divv, T) == 0);

    using constants::CD;

    auto poly = K * (K + CD<3.>()) +
                (K + CD<0.5>()) * (K + CD<1.5>()) * (K + CD<2.5>());
    // static_assert(order<decltype(poly), decltype(K)>::call() == 3);
    static_assert(order(poly, K) == 3);
}

TEST(Order, AddSub) {
    ADHOC(var1);
    ADHOC(var2);

    auto var3 = var1 + var2;
    auto var4 = var3 - var2;
    // static_assert(order<decltype(var3), decltype(var1)>::call() == 1);
    // static_assert(order<decltype(var3), decltype(var2)>::call() ==
    //               std::numeric_limits<std::size_t>::max());
    static_assert(order(var4, var1) == 1);
    static_assert(order(var4, var2) == 1);
}

TEST(Order, Vars) {
    ADHOC(var1);
    ADHOC(var2);
    ADHOC(var33);

    static_assert(var1 == var1);
    static_assert(var1 != var2);
    static_assert(var1 != var33);

    static_assert(depends(var1, var1));
    static_assert(!depends(var1, var2));
    static_assert(!depends(var1, var33));

    auto ret1 = exp(var1);
    static_assert(depends(ret1, var1));
    static_assert(!depends(ret1, var2));

    auto ret2 = var1 + var2;
    static_assert(depends(ret2, var1));
    static_assert(depends(ret2, var2));
    static_assert(!depends(ret2, var33));

    auto ret3 = exp(ret2);

    auto result = order_vars(ret3, var1, var2, var33);
    std::cout << type_name2<decltype(result)>() << std::endl;

    auto ret4 = var1 * var2;
    auto result2 = order_vars(var1 * var2, var1, var2, var33);
    // auto result2 = order_vars(ret4, var1);
    std::cout << type_name2<decltype(result2)>() << std::endl;

    auto ret5 = var1 + var2;
    auto result3 = expand(ret5);
    std::cout << type_name2<decltype(result3)>() << std::endl;
}

TEST(Order, Expand) {
    ADHOC(var1);
    ADHOC(var2);
    ADHOC(var3);
    ADHOC(var4);

    // {
    //     auto var = var3 * (var1 + var2);
    //     auto result = expand(var);
    //     std::cout << type_name2<decltype(result)>() << std::endl;
    // }

    // {
    //     auto var = (var1 + var2) * var3;
    //     auto result = expand(var);
    //     std::cout << type_name2<decltype(result)>() << std::endl;
    // }

    // {
    //     auto var = (var1 + var2) * (var3 + var4);
    //     auto result = expand(var);
    //     std::cout << type_name2<decltype(result)>() << std::endl;
    // }

    // {
    //     auto var = (var1 + var2) * (var3 + var4) * exp(var1);

    //     auto result2 = order_vars2(var, var3, var2);
    //     std::cout << type_name2<decltype(result2)>() << std::endl;
    // }

    {
        auto var = (var1 + (var3 * var4)) * (var3 + (var3 * var4)) * exp(var1);
        // auto var = var1 + (var1 + var3) * (var3);

        // auto var = var1 + var1 * (var3);

        auto vars = std::tuple(var1, var2, var3, var4);
        auto result2 = convex_hull_orders(var, vars);
        std::cout << type_name2<decltype(result2)>() << std::endl;
    }

    auto v1 = std::index_sequence<18446744073709551615UL, 0, 1, 0>{};
    auto v2 = std::index_sequence<18446744073709551615UL, 0, 1, 1>{};
    auto v3 = std::index_sequence<18446744073709551615UL, 0, 0, 1>{};

    // bool res = v1 <= v2;
    std::cout << (v1 <= v2) << std::endl;
    std::cout << (v2 <= v1) << std::endl;
    std::cout << (v1 <= v1) << std::endl;
    std::cout << (v2 <= v2) << std::endl;

    std::cout << (v1 <= v3) << std::endl;
    std::cout << (v3 <= v1) << std::endl;

    auto vartoreduce = std::make_tuple(v1, v2, v3);
    auto varreduced = filter_smaller(vartoreduce);
    std::cout << type_name2<decltype(vartoreduce)>() << std::endl;
    std::cout << type_name2<decltype(varreduced)>() << std::endl;
}

TEST(Order, GetVars) {
    ADHOC(var1);
    ADHOC(var2);
    ADHOC(var3);
    auto res = d<2>(var1) * d(var2);
    auto res2 = d<2>(var1) * d(var3);

    auto vars = get_vars_single(res);
    auto vars2 = get_vars_single(res2);
    std::cout << type_name2<decltype(vars)>() << std::endl;

    auto res3 = std::tuple(res, res2);
    auto vars3 = get_vars(res3);
    std::cout << type_name2<decltype(vars3)>() << std::endl;
}

TEST(Order, MinimumPolys) {
    std::tuple<std::index_sequence<1, 2, 0>, std::index_sequence<2, 0, 0>,
               std::index_sequence<0, 1, 1>>
        val;
    auto result_exp = get_minimum_polys(val);
    std::cout << type_name2<decltype(result_exp)>() << std::endl;
}

TEST(Order, MonomialIncluded) {
    ADHOC(var1);
    ADHOC(var2);
    ADHOC(var3);
    ADHOC(var4);
    auto res = d<2>(var1) * d(var2);
    auto res2 = d<2>(var1) * d(var3);

    auto res3 = std::tuple(res, res2);

    // auto var = d((var1 + (var3 * var4)) * (var3 + (var3 * var4)) *
    // exp(var1));
    // auto var = d(var1 + (var3 * var4)) * d<2>(var1);
    auto var = d<2>(var1) * d(var1 + (var3 * var4));

    constexpr auto resbool = monomial_included(var, res3);
    static_assert(resbool);
}

TEST(Order, MonomialNotIncluded) {
    ADHOC(var1);
    ADHOC(var2);

    auto res = d<2>(var1);
    auto res2 = d<2>(var2);
    auto res_t = std::tuple(res, res2);
    auto out1 = d(exp(var1)) * d(exp(var2));
    constexpr auto resbool = monomial_included(out1, res_t);
    static_assert(!resbool);

    ADHOC(var3);

    auto res3 = d<2>(var3) * d<2>(var1);
    auto res2_t = std::tuple(res2, res3);
    auto out2 = d(exp(var1)) * d(exp(var3));
    constexpr auto resbool2 = monomial_included(out2, res2_t);
    static_assert(resbool2);
}

} // namespace adhoc4
