#include <adhoc.hpp>
#include <monomial_included.hpp>
#include <utils/flatten.hpp>

int main() {
    using namespace adhoc;

    // FlattenTest
    {
        constexpr auto in = std::make_tuple(
            std::make_tuple(std::index_sequence<1, 0, 0, 0, 0, 0>{},
                            std::index_sequence<0, 1, 0, 0, 0, 0>{}),
            std::make_tuple(std::index_sequence<0, 0, 1, 0, 0, 0>{},
                            std::index_sequence<0, 0, 0, 1, 0, 0>{}),
            std::make_tuple(std::index_sequence<0, 0, 0, 0, 1, 0>{},
                            std::index_sequence<0, 0, 0, 0, 0, 1>{}));

        constexpr auto out = flatten(in);

        constexpr auto out2 = std::make_tuple(
            std::make_tuple(std::index_sequence<1, 0, 0, 0, 0, 0>{},
                            std::index_sequence<0, 0, 1, 0, 0, 0>{},
                            std::index_sequence<0, 0, 0, 0, 1, 0>{}),
            std::make_tuple(std::index_sequence<1, 0, 0, 0, 0, 0>{},
                            std::index_sequence<0, 0, 1, 0, 0, 0>{},
                            std::index_sequence<0, 0, 0, 0, 0, 1>{}),
            std::make_tuple(std::index_sequence<1, 0, 0, 0, 0, 0>{},
                            std::index_sequence<0, 0, 0, 1, 0, 0>{},
                            std::index_sequence<0, 0, 0, 0, 1, 0>{}),
            std::make_tuple(std::index_sequence<1, 0, 0, 0, 0, 0>{},
                            std::index_sequence<0, 0, 0, 1, 0, 0>{},
                            std::index_sequence<0, 0, 0, 0, 0, 1>{}),
            std::make_tuple(std::index_sequence<0, 1, 0, 0, 0, 0>{},
                            std::index_sequence<0, 0, 1, 0, 0, 0>{},
                            std::index_sequence<0, 0, 0, 0, 1, 0>{}),
            std::make_tuple(std::index_sequence<0, 1, 0, 0, 0, 0>{},
                            std::index_sequence<0, 0, 1, 0, 0, 0>{},
                            std::index_sequence<0, 0, 0, 0, 0, 1>{}),
            std::make_tuple(std::index_sequence<0, 1, 0, 0, 0, 0>{},
                            std::index_sequence<0, 0, 0, 1, 0, 0>{},
                            std::index_sequence<0, 0, 0, 0, 1, 0>{}),
            std::make_tuple(std::index_sequence<0, 1, 0, 0, 0, 0>{},
                            std::index_sequence<0, 0, 0, 1, 0, 0>{},
                            std::index_sequence<0, 0, 0, 0, 0, 1>{}));

        static_assert(std::is_same_v<decltype(out), decltype(out2)>);
    }

    // ErrorFix
    {
        ADHOC(v1);
        auto res = erfc(v1);
        auto out1 = d<1>(res) * d(v1);
        {
            auto res = d(v1);
            auto res_t = std::make_tuple(
                std::make_pair(res, std::integral_constant<std::size_t, 2>{}));
            constexpr auto resbool = monomial_included(out1, res_t);
            static_assert(!resbool);
        }

        {
            auto res = d<2>(v1);
            auto res_t = std::make_tuple(
                std::make_pair(res, std::integral_constant<std::size_t, 2>{}));
            constexpr auto resbool = monomial_included(out1, res_t);
            static_assert(resbool);
        }

        {
            auto res3 = d<3>(v1);
            auto res_t = std::make_tuple(
                std::make_pair(res3, std::integral_constant<std::size_t, 2>{}));
            constexpr auto resbool = monomial_included(out1, res_t);
            static_assert(resbool);
        }
    }

    // SumOrder2
    {
        ADHOC(x);
        ADHOC(y);
        constexpr auto resbool = monomial_included(
            d<2>(x + y),
            std::make_tuple(
                std::make_pair(d(x), std::integral_constant<std::size_t, 1>{}),
                std::make_pair(d(y), std::integral_constant<std::size_t, 2>{}),
                std::make_pair(d(x) * d(y),
                               std::integral_constant<std::size_t, 3>{})));

        static_assert(resbool);
    }

    // ErrorFix2
    {
        ADHOC(S);
        ADHOC(K);
        ADHOC(v);

        {
            constexpr auto d1 = d(K * S) * d(sqrt(v));
            constexpr auto in = d<2>(v);

            constexpr auto resbool = monomial_included(
                d1, std::make_tuple(std::make_pair(
                        in, std::integral_constant<std::size_t, 2>{})));
            static_assert(!resbool);
        }

        {
            constexpr auto d1 = d(K /* * S */) * d(v * v);
            constexpr auto in = d<2>(v);

            constexpr auto resbool = monomial_included(
                d1, std::make_tuple(std::make_pair(
                        in, std::integral_constant<std::size_t, 2>{})));
            static_assert(!resbool);
        }

        {
            constexpr auto d1 = d(K * K) * d(v * v);
            constexpr auto in = d<2>(v);
            constexpr auto in2 = d<2>(K);

            constexpr auto resbool = monomial_included(
                d1, std::make_tuple(
                        std::make_pair(
                            in, std::integral_constant<std::size_t, 2>{}),
                        std::make_pair(
                            in2, std::integral_constant<std::size_t, 3>{})));
            static_assert(!resbool);
        }
    }

    // ErrorFix3
    {

        ADHOC(x);
        ADHOC(y);

        auto inputs_do = std::make_tuple(
            std::make_pair(d(x), std::integral_constant<std::size_t, 1>{}),
            std::make_pair(d(y), std::integral_constant<std::size_t, 2>{}));

        auto d1 = d(x / y);

        constexpr auto resbool = monomial_included(d1, inputs_do);
        static_assert(resbool);
    }

    // OldTest
    {
        ADHOC(x);
        ADHOC(y);
        constexpr auto dnin = std::make_tuple(
            std::make_pair(d(x), std::integral_constant<std::size_t, 1>{}),
            std::make_pair(d(y), std::integral_constant<std::size_t, 2>{}),
            std::make_pair(d(x) * d(y),
                           std::integral_constant<std::size_t, 3>{}));
        constexpr auto node = d(x) * d<2>(y);
        constexpr auto resbool = monomial_included(node, dnin);
        static_assert(!resbool);
    }

    // MonomialIncluded
    {
        ADHOC(var1);
        ADHOC(var2);
        ADHOC(var3);
        ADHOC(var4);
        auto res = d<2>(var1) * d(var2);
        auto res2 = d<2>(var1) * d(var3);

        auto res3 = std::tuple(
            std::make_pair(res, std::integral_constant<std::size_t, 1>{}),
            std::make_pair(res2, std::integral_constant<std::size_t, 2>{}));

        auto var = d<2>(var1) * d(var1 + (var3 * var4));

        constexpr auto resbool = monomial_included(var, res3);
        static_assert(resbool);
    }

    // MonomialNotIncluded
    {
        ADHOC(var1);
        ADHOC(var2);

        auto res = d<2>(var1);
        auto res2 = d<2>(var2);
        auto res_t = std::tuple(
            std::make_pair(res, std::integral_constant<std::size_t, 1>{}),
            std::make_pair(res2, std::integral_constant<std::size_t, 2>{}));
        auto out1 = d(exp(var1)) * d(exp(var2));
        constexpr auto resbool = monomial_included(out1, res_t);
        static_assert(!resbool);

        ADHOC(var3);

        auto res3 = d<2>(var3) * d<2>(var1);
        auto res2_t = std::tuple(
            std::make_pair(res2, std::integral_constant<std::size_t, 2>{}),
            std::make_pair(res3, std::integral_constant<std::size_t, 3>{}));
        auto out2 = d(exp(var1)) * d(exp(var3));
        constexpr auto resbool2 = monomial_included(out2, res2_t);
        static_assert(resbool2);
    }

    return 0;
}
