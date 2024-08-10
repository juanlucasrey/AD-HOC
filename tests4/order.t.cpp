#include <adhoc.hpp>
#include <constants_type.hpp>
#include <order.hpp>

#include <gtest/gtest.h>

namespace adhoc4 {

TEST(Order, Const) {
    using constants::CD;

    ADHOC(var);
    auto constant = CD<-0.5>();

    static_assert(order<decltype(constant), decltype(var)>::call() == 0);
}

TEST(Order, VarType) {
    ADHOC(var);
    ADHOC(var2);

    static_assert(order<decltype(var), decltype(var)>::call() == 1);
    static_assert(order<decltype(var), decltype(var2)>::call() == 0);
}

TEST(Order, UnivariateType) {
    ADHOC(var);

    auto var2 = exp(var);
    static_assert(order<decltype(var2), decltype(var)>::call() ==
                  std::numeric_limits<std::size_t>::max());
}

TEST(Order, Add) {
    ADHOC(var1);
    ADHOC(var2);

    auto var3 = var1 + var2;
    static_assert(order<decltype(var3), decltype(var1)>::call() == 1);
    static_assert(order<decltype(var3), decltype(var2)>::call() == 1);
}

TEST(Order, Sub) {
    ADHOC(var1);
    ADHOC(var2);

    auto var3 = var1 - var2;
    static_assert(order<decltype(var3), decltype(var1)>::call() == 1);
    static_assert(order<decltype(var3), decltype(var2)>::call() == 1);
}

TEST(Order, Mul) {
    ADHOC(var1);
    ADHOC(var2);

    auto var3 = var1 * var2;
    static_assert(order<decltype(var3), decltype(var1)>::call() == 1);
    static_assert(order<decltype(var3), decltype(var2)>::call() == 1);

    auto var4 = var1 * var1;
    static_assert(order<decltype(var4), decltype(var1)>::call() == 2);
}

TEST(Order, Div) {
    ADHOC(var1);
    ADHOC(var2);

    auto var3 = var1 / var2;
    static_assert(order<decltype(var3), decltype(var1)>::call() == 1);
    static_assert(order<decltype(var3), decltype(var2)>::call() ==
                  std::numeric_limits<std::size_t>::max());
}

TEST(Order, First) {
    ADHOC(S);
    ADHOC(K);
    ADHOC(v);
    ADHOC(T);

    auto var = exp(T);
    static_assert(order<decltype(var), decltype(T)>::call() ==
                  std::numeric_limits<std::size_t>::max());

    auto sumv = K + v;
    static_assert(order<decltype(sumv), decltype(S)>::call() == 0);
    static_assert(order<decltype(sumv), decltype(K)>::call() == 1);
    static_assert(order<decltype(sumv), decltype(v)>::call() == 1);
    static_assert(order<decltype(sumv), decltype(T)>::call() == 0);

    auto mulv = K * v * K;
    static_assert(order<decltype(sumv), decltype(S)>::call() == 0);
    static_assert(order<decltype(mulv), decltype(K)>::call() == 2);
    static_assert(order<decltype(mulv), decltype(v)>::call() == 1);
    static_assert(order<decltype(sumv), decltype(T)>::call() == 0);

    auto mulexpv = K * v * exp(K);
    static_assert(order<decltype(mulexpv), decltype(S)>::call() == 0);
    static_assert(order<decltype(mulexpv), decltype(K)>::call() ==
                  std::numeric_limits<std::size_t>::max());
    static_assert(order<decltype(mulexpv), decltype(v)>::call() == 1);
    static_assert(order<decltype(mulexpv), decltype(T)>::call() == 0);

    auto divv = v / K;
    static_assert(order<decltype(divv), decltype(S)>::call() == 0);
    static_assert(order<decltype(divv), decltype(K)>::call() ==
                  std::numeric_limits<std::size_t>::max());
    static_assert(order<decltype(divv), decltype(v)>::call() == 1);
    static_assert(order<decltype(divv), decltype(T)>::call() == 0);

    using constants::CD;

    auto poly = K * (K + CD<3.>()) +
                (K + CD<0.5>()) * (K + CD<1.5>()) * (K + CD<2.5>());
    static_assert(order<decltype(poly), decltype(K)>::call() == 3);
}

} // namespace adhoc4
