
#include <adhoc.hpp>
#include <constants_type.hpp>
#include <differential_operator.hpp>
#include <monomial_included.hpp>
#include <order.hpp>

int main() {
    using namespace adhoc4;

    // Inverse
    {

        constexpr auto max = std::numeric_limits<std::size_t>::max();

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

    // Const
    {
        using constants::CD;

        ADHOC(var);
        auto constant = CD<-0.5>();

        static_assert(order(constant, var) == 0);
    }

    // VarType
    {
        ADHOC(var);
        ADHOC(var2);

        static_assert(order(var, var) == 1);
        static_assert(order(var, var2) == 0);
    }

    // UnivariateType
    {
        ADHOC(var);

        auto var2 = exp(var);
        static_assert(order(var2, var) ==
                      std::numeric_limits<std::size_t>::max());
    }

    // Add
    {
        ADHOC(var1);
        ADHOC(var2);

        auto var3 = var1 + var2;
        static_assert(order(var3, var1) == 1);
        static_assert(order(var3, var1) == 1);
    }

    // Sub
    {
        ADHOC(var1);
        ADHOC(var2);

        auto var3 = var1 - var2;
        static_assert(order(var3, var1) == 1);
        static_assert(order(var3, var1) == 1);
    }

    // Mul
    {
        ADHOC(var1);
        ADHOC(var2);

        auto var3 = var1 * var2;
        static_assert(order(var3, var1) == 1);
        static_assert(order(var3, var1) == 1);

        auto var4 = var1 * var1;
        static_assert(order(var4, var1) == 2);
    }

    // Div
    {
        ADHOC(var1);
        ADHOC(var2);

        auto var3 = var1 / var2;
        static_assert(order(var3, var1) == 1);
        static_assert(order(var3, var2) ==
                      std::numeric_limits<std::size_t>::max());
    }

    // First
    {
        ADHOC(S);
        ADHOC(K);
        ADHOC(v);
        ADHOC(T);

        auto var = exp(T);
        static_assert(order(var, T) == std::numeric_limits<std::size_t>::max());

        auto sumv = K + v;
        static_assert(order(sumv, S) == 0);
        static_assert(order(sumv, K) == 1);
        static_assert(order(sumv, v) == 1);
        static_assert(order(sumv, T) == 0);

        auto mulv = K * v * K;
        static_assert(order(mulv, S) == 0);
        static_assert(order(mulv, K) == 2);
        static_assert(order(mulv, v) == 1);
        static_assert(order(mulv, T) == 0);

        auto mulexpv = K * v * exp(K);
        static_assert(order(mulexpv, S) == 0);
        static_assert(order(mulexpv, K) ==
                      std::numeric_limits<std::size_t>::max());
        static_assert(order(mulexpv, v) == 1);
        static_assert(order(mulexpv, T) == 0);

        auto divv = v / K;
        static_assert(order(divv, S) == 0);
        static_assert(order(divv, K) ==
                      std::numeric_limits<std::size_t>::max());
        static_assert(order(divv, v) == 1);
        static_assert(order(divv, T) == 0);

        using constants::CD;

        auto poly = K * (K + CD<3.>()) +
                    (K + CD<0.5>()) * (K + CD<1.5>()) * (K + CD<2.5>());
        static_assert(order(poly, K) == 3);
    }

    // AddSub
    {
        ADHOC(var1);
        ADHOC(var2);

        auto var3 = var1 + var2;
        auto var4 = var3 - var2;
        static_assert(order(var4, var1) == 1);
        static_assert(order(var4, var2) == 1);
    }

    // Vars
    {
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

        {
            auto ret3 = exp(ret2);
            auto result = order_vars(ret3, var1, var2, var33);
            std::index_sequence<std::numeric_limits<std::size_t>::max(),
                                std::numeric_limits<std::size_t>::max(), 0>
                res_check;
            static_assert(
                std::is_same_v<decltype(result), decltype(res_check)>);
        }

        {
            auto result = order_vars(var1 * var2, var1, var2, var33);
            std::index_sequence<1, 1, 0> res_check;
            static_assert(
                std::is_same_v<decltype(result), decltype(res_check)>);
        }

        {
            auto ret5 = var1 + var2;
            auto result = expand(ret5);
            auto res_check = std::make_tuple(var1, var2);
            static_assert(
                std::is_same_v<decltype(result), decltype(res_check)>);
        }
    }

    // Expand
    {
        ADHOC(var1);
        ADHOC(var2);
        ADHOC(var3);
        ADHOC(var4);

        {
            auto var =
                (var1 + (var3 * var4)) * (var3 + (var3 * var4)) * exp(var1);
            auto vars = std::tuple(var1, var2, var3, var4);
            auto result2 = convex_hull_orders(var, vars);

            std::tuple<std::index_sequence<
                std::numeric_limits<std::size_t>::max(), 0, 2, 2>>
                res_check;

            static_assert(
                std::is_same_v<decltype(result2), decltype(res_check)>);
        }

        {
            auto v1 = std::index_sequence<18446744073709551615UL, 0, 1, 0>{};
            auto v2 = std::index_sequence<18446744073709551615UL, 0, 1, 1>{};
            auto v3 = std::index_sequence<18446744073709551615UL, 0, 0, 1>{};

            static_assert(v1 <= v2);
            static_assert(!(v2 <= v1));
            static_assert(v1 <= v1);
            static_assert(v2 <= v2);

            static_assert(!(v1 <= v3));
            static_assert(!(v3 <= v1));

            auto vartoreduce = std::make_tuple(v1, v2, v3);
            auto varreduced = filter_smaller(vartoreduce);
            std::tuple<std::index_sequence<
                std::numeric_limits<std::size_t>::max(), 0, 1, 1>>
                res_check;
            static_assert(
                std::is_same_v<decltype(varreduced), decltype(res_check)>);
        }
    }

    // MinimumPolys
    {
        std::tuple<std::index_sequence<1, 2, 0>, std::index_sequence<2, 0, 0>,
                   std::index_sequence<0, 1, 1>>
            val;
        auto result_exp = get_minimum_polys(val);

        std::tuple<std::index_sequence<2, 1, 0>, std::index_sequence<2, 0, 1>,
                   std::index_sequence<1, 2, 0>, std::index_sequence<1, 1, 1>>
            res_check;

        static_assert(
            std::is_same_v<decltype(result_exp), decltype(res_check)>);
    }
    return 0;
}
