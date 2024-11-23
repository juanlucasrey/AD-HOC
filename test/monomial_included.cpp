#include <adhoc.hpp>
#include <backpropagator_tools.hpp>
#include <calc_tree.hpp>
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

        // dependency on calc_tree should be removed
        CalcTree ct(res);
        using NodesValue = decltype(ct)::ValuesTupleInverse;

        auto out1 = d<1>(res) * d(v1);
        auto out1_i = detail::convert_to_index2<NodesValue>(out1);
        {
            auto res = d(v1);
            auto res_t = std::make_tuple(
                std::make_pair(res, std::integral_constant<std::size_t, 2>{}));
            auto res_t_i = detail::convert_to_index_many<NodesValue>(res_t);

            constexpr auto resbool =
                monomial_included2<NodesValue>(out1_i, res_t_i);
            static_assert(!resbool);
        }

        {
            auto res = d<2>(v1);
            auto res_t = std::make_tuple(
                std::make_pair(res, std::integral_constant<std::size_t, 2>{}));
            auto res_t_i = detail::convert_to_index_many<NodesValue>(res_t);

            constexpr auto resbool =
                monomial_included2<NodesValue>(out1_i, res_t_i);
            static_assert(resbool);
        }

        {
            auto res3 = d<3>(v1);
            auto res_t = std::make_tuple(
                std::make_pair(res3, std::integral_constant<std::size_t, 2>{}));
            auto res_t_i = detail::convert_to_index_many<NodesValue>(res_t);

            constexpr auto resbool =
                monomial_included2<NodesValue>(out1_i, res_t_i);
            static_assert(resbool);
        }
    }

    // SumOrder2
    {
        ADHOC(x);
        ADHOC(y);
        auto res = x + y;
        CalcTree ct(res);
        using NodesValue = decltype(ct)::ValuesTupleInverse;

        auto out1 = d<2>(res);
        auto out1_i = detail::convert_to_index2<NodesValue>(out1);

        auto inputs = std::make_tuple(
            std::make_pair(d(x), std::integral_constant<std::size_t, 1>{}),
            std::make_pair(d(y), std::integral_constant<std::size_t, 2>{}),
            std::make_pair(d(x) * d(y),
                           std::integral_constant<std::size_t, 3>{}));
        auto inputs_i = detail::convert_to_index_many<NodesValue>(inputs);

        constexpr auto resbool =
            monomial_included2<NodesValue>(out1_i, inputs_i);
        static_assert(resbool);
    }

    // ErrorFix2
    {
        ADHOC(S);
        ADHOC(K);
        ADHOC(v);

        {
            auto res = (K * S) * sqrt(v);
            CalcTree ct(res);
            using NodesValue = decltype(ct)::ValuesTupleInverse;

            constexpr auto d1 = d(K * S) * d(sqrt(v));
            constexpr auto in = d<2>(v);

            constexpr auto inputs = std::make_tuple(
                std::make_pair(in, std::integral_constant<std::size_t, 2>{}));

            auto d1_i = detail::convert_to_index2<NodesValue>(d1);
            auto inputs_i = detail::convert_to_index_many<NodesValue>(inputs);

            constexpr auto resbool =
                monomial_included2<NodesValue>(d1_i, inputs_i);
            static_assert(!resbool);
        }

        {
            auto res = K * (v * v);
            CalcTree ct(res);
            using NodesValue = decltype(ct)::ValuesTupleInverse;

            constexpr auto d1 = d(K /* * S */) * d(v * v);
            constexpr auto in = d<2>(v);

            constexpr auto inputs = std::make_tuple(
                std::make_pair(in, std::integral_constant<std::size_t, 2>{}));

            auto d1_i = detail::convert_to_index2<NodesValue>(d1);
            auto inputs_i = detail::convert_to_index_many<NodesValue>(inputs);

            constexpr auto resbool =
                monomial_included2<NodesValue>(d1_i, inputs_i);
            static_assert(!resbool);
        }

        {
            auto res = (K * K) * (v * v);
            CalcTree ct(res);
            using NodesValue = decltype(ct)::ValuesTupleInverse;

            constexpr auto d1 = d(K * K) * d(v * v);
            constexpr auto in = d<2>(v);
            constexpr auto in2 = d<2>(K);

            constexpr auto inputs = std::make_tuple(
                std::make_pair(in, std::integral_constant<std::size_t, 2>{}),
                std::make_pair(in2, std::integral_constant<std::size_t, 3>{}));

            auto d1_i = detail::convert_to_index2<NodesValue>(d1);
            auto inputs_i = detail::convert_to_index_many<NodesValue>(inputs);

            constexpr auto resbool =
                monomial_included2<NodesValue>(d1_i, inputs_i);
            static_assert(!resbool);
        }
    }

    // ErrorFix3
    {

        ADHOC(x);
        ADHOC(y);

        auto res = x / y;
        CalcTree ct(res);
        using NodesValue = decltype(ct)::ValuesTupleInverse;

        auto inputs_do = std::make_tuple(
            std::make_pair(d(x), std::integral_constant<std::size_t, 1>{}),
            std::make_pair(d(y), std::integral_constant<std::size_t, 2>{}));

        auto d1 = d(x / y);

        auto d1_i = detail::convert_to_index2<NodesValue>(d1);
        auto inputs_i = detail::convert_to_index_many<NodesValue>(inputs_do);

        constexpr auto resbool = monomial_included2<NodesValue>(d1_i, inputs_i);
        static_assert(resbool);
    }

    // OldTest
    {
        ADHOC(x);
        ADHOC(y);

        auto res = x * y;
        CalcTree ct(res);
        using NodesValue = decltype(ct)::ValuesTupleInverse;

        constexpr auto dnin = std::make_tuple(
            std::make_pair(d(x), std::integral_constant<std::size_t, 1>{}),
            std::make_pair(d(y), std::integral_constant<std::size_t, 2>{}),
            std::make_pair(d(x) * d(y),
                           std::integral_constant<std::size_t, 3>{}));
        constexpr auto node = d(x) * d<2>(y);

        auto node_i = detail::convert_to_index2<NodesValue>(node);
        auto dnin_i = detail::convert_to_index_many<NodesValue>(dnin);

        constexpr auto resbool = monomial_included2<NodesValue>(node_i, dnin_i);
        static_assert(!resbool);
    }

    // MonomialIncluded
    {
        ADHOC(var1);
        ADHOC(var2);
        ADHOC(var3);
        ADHOC(var4);

        auto resct = var1 + (var3 * var4) + var2;
        CalcTree ct(resct);
        using NodesValue = decltype(ct)::ValuesTupleInverse;

        auto res = d<2>(var1) * d(var2);
        auto res2 = d<2>(var1) * d(var3);

        auto res3 = std::tuple(
            std::make_pair(res, std::integral_constant<std::size_t, 1>{}),
            std::make_pair(res2, std::integral_constant<std::size_t, 2>{}));

        auto var = d<2>(var1) * d(var1 + (var3 * var4));

        auto var_i = detail::convert_to_index2<NodesValue>(var);
        auto res3_i = detail::convert_to_index_many<NodesValue>(res3);

        constexpr auto resbool = monomial_included2<NodesValue>(var_i, res3_i);
        static_assert(resbool);
    }

    // MonomialNotIncluded
    {
        ADHOC(var1);
        ADHOC(var2);

        auto resct = exp(var1) * exp(var2);
        CalcTree ct(resct);
        using NodesValue = decltype(ct)::ValuesTupleInverse;

        auto res = d<2>(var1);
        auto res2 = d<2>(var2);
        auto res_t = std::tuple(
            std::make_pair(res, std::integral_constant<std::size_t, 1>{}),
            std::make_pair(res2, std::integral_constant<std::size_t, 2>{}));
        auto out1 = d(exp(var1)) * d(exp(var2));

        auto out1_i = detail::convert_to_index2<NodesValue>(out1);
        auto res_t_i = detail::convert_to_index_many<NodesValue>(res_t);

        constexpr auto resbool =
            monomial_included2<NodesValue>(out1_i, res_t_i);
        static_assert(!resbool);

        ADHOC(var3);

        auto resct2 = exp(var1) * exp(var2) * exp(var3);
        CalcTree ct2(resct2);
        using NodesValue2 = decltype(ct2)::ValuesTupleInverse;

        auto res3 = d<2>(var3) * d<2>(var1);
        auto res2_t = std::tuple(
            std::make_pair(res2, std::integral_constant<std::size_t, 2>{}),
            std::make_pair(res3, std::integral_constant<std::size_t, 3>{}));
        auto out2 = d(exp(var1)) * d(exp(var3));

        auto out2_i = detail::convert_to_index2<NodesValue2>(out2);
        auto res2_t_i = detail::convert_to_index_many<NodesValue2>(res2_t);

        constexpr auto resbool2 =
            monomial_included2<NodesValue2>(out2_i, res2_t_i);
        static_assert(resbool2);
    }

    return 0;
}
