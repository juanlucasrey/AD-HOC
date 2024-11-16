#include <adhoc.hpp>
#include <backpropagator_tools.hpp>
#include <calc_tree.hpp>
#include <differential_operator.hpp>
#include <sort.hpp>
#include <utils/index_sequence.hpp>
#include <utils/tuple.hpp>

#include "call_price.hpp"

#include <tuple>
#include <type_traits>
#include <utility>

int main() {
    using namespace adhoc;

    // MultiplyOrdered
    {
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
            static_assert(
                std::is_same_v<decltype(diffop1mul), decltype(result)>);
        }

        {
            constexpr auto dS = d(S);
            constexpr auto dK = d(K);
            constexpr auto diffop1mul =
                detail::multiply_ordered(dS, dK, NodesValue{});

            constexpr auto result = d(S) * d(K);
            static_assert(
                std::is_same_v<decltype(diffop1mul), decltype(result)>);
        }

        {
            constexpr auto dS = d(S);
            constexpr auto dK = d(K) * d(T);
            constexpr auto diffop1mul =
                detail::multiply_ordered(dS, dK, NodesValue{});

            constexpr auto result = d(S) * d(K) * d(T);
            static_assert(
                std::is_same_v<decltype(diffop1mul), decltype(result)>);
        }

        {
            constexpr auto dS = d(S) * d(v);
            constexpr auto dK = d(K) * d(T);
            constexpr auto diffop1mul =
                detail::multiply_ordered(dS, dK, NodesValue{});

            constexpr auto result = d(S) * d(K) * d(v) * d(T);
            static_assert(
                std::is_same_v<decltype(diffop1mul), decltype(result)>);
        }

        {
            constexpr auto dS = d(S) * d<2>(K) * d(v);
            constexpr auto dK = d(K) * d(T);
            constexpr auto diffop1mul =
                detail::multiply_ordered(dS, dK, NodesValue{});

            constexpr auto result = d(S) * d<3>(K) * d(v) * d(T);
            static_assert(
                std::is_same_v<decltype(diffop1mul), decltype(result)>);
        }
    }

    // IsInputDer
    {
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

    // MergeSorted2
    {

        const std::tuple<
            std::pair<std::tuple<adhoc::der::d<
                          1, std::integral_constant<unsigned long, 2>>>,
                      std::pair<adhoc::detail::on_buffer_t,
                                std::integral_constant<unsigned long, 1>>>>
            next_derivatives_new_with_pos_idx;

        std::tuple<std::pair<
            std::tuple<
                adhoc::der::d<1, std::integral_constant<unsigned long, 2>>,
                adhoc::der::d<1, std::integral_constant<unsigned long, 3>>>,
            std::pair<adhoc::detail::on_buffer_t,
                      std::integral_constant<unsigned long, 0>>>>
            dnn_i;

        std::tuple<
            adhoc::erfc_t<adhoc::mul_t<
                adhoc::double_t<detail::StringLiteral<2>{"x"}>,
                adhoc::inv_t<adhoc::double_t<detail::StringLiteral<2>{"y"}>>>>,
            adhoc::mul_t<
                adhoc::double_t<detail::StringLiteral<2>{"x"}>,
                adhoc::inv_t<adhoc::double_t<detail::StringLiteral<2>{"y"}>>>,
            adhoc::inv_t<adhoc::double_t<detail::StringLiteral<2>{"y"}>>,
            adhoc::double_t<detail::StringLiteral<2>{"x"}>,
            adhoc::double_t<detail::StringLiteral<2>{"y"}>>
            nodes;

        constexpr auto dnn_i_new = detail::merge_sorted(
            next_derivatives_new_with_pos_idx, dnn_i, nodes);

        constexpr std::tuple<
            std::pair<
                std::tuple<
                    adhoc::der::d<1, std::integral_constant<unsigned long, 2>>,
                    adhoc::der::d<1, std::integral_constant<unsigned long, 3>>>,
                std::pair<adhoc::detail::on_buffer_t,
                          std::integral_constant<unsigned long, 0>>>,
            std::pair<std::tuple<adhoc::der::d<
                          1, std::integral_constant<unsigned long, 2>>>,
                      std::pair<adhoc::detail::on_buffer_t,
                                std::integral_constant<unsigned long, 1>>>>
            res_check;

        static_assert(std::is_same_v<decltype(dnn_i_new), decltype(res_check)>);
    }

    // SortPair
    {
        ADHOC(x);
        ADHOC(y);
        auto res = log(x + y);

        CalcTree ct(res);
        using NodesValue = decltype(ct)::ValuesTupleInverse;

        constexpr auto result = std::make_tuple(x, y);
        constexpr auto pair = std::make_tuple(x, y);
        constexpr auto ordered_pair = detail::sort_pair(pair, NodesValue{});
        static_assert(ordered_pair == result);

        constexpr auto inverted_pair = std::make_tuple(y, x);
        constexpr auto ordered_pair_2 =
            detail::sort_pair(inverted_pair, NodesValue{});
        static_assert(ordered_pair_2 == result);
    }

    // PowDiffOp
    {
        ADHOC(x);
        ADHOC(y);

        {
            constexpr auto res = d<4>(x) * d(y);
            constexpr auto res2 = detail::elevate_by_power<0>(res);
            static_assert(std::is_same_v<const std::tuple<>, decltype(res2)>);
        }

        {
            constexpr auto vars = std::make_tuple(d(x) * d(y), d(x), d(y));
            constexpr auto is = std::index_sequence<1, 2, 3>{};
            constexpr auto res3 = create_differential_operator(vars, is);

            constexpr auto res_check = d<3>(x) * d<4>(y);
            static_assert(std::is_same_v<decltype(res_check), decltype(res3)>);
        }

        {
            constexpr auto var1 = detail::elevate_by_power<1>(d(x) * d(y));
            constexpr auto var2 = detail::elevate_by_power<2>(d(x));
            constexpr auto var3 = var1 * var2;

            constexpr auto res_check = d<3>(x) * d(y);
            static_assert(std::is_same_v<decltype(res_check), decltype(var3)>);
        }

        {
            constexpr auto var1 = d(x) * d(y);
            constexpr auto var2 = d<2>(x);
            constexpr auto var3 = var1 * var2;

            constexpr auto res_check = d<3>(x) * d(y);
            static_assert(std::is_same_v<decltype(res_check), decltype(var3)>);
        }

        {
            constexpr auto var4 = (d(x) * d(y)) * d(x);
            constexpr auto res_check = d<2>(x) * d(y);
            static_assert(std::is_same_v<decltype(res_check), decltype(var4)>);
        }
        {
            constexpr auto vars = std::make_tuple(d(x), d(y));
            constexpr auto is1 = std::index_sequence<4, 1>{};
            constexpr auto diff_op = create_differential_operator(vars, is1);

            constexpr auto res = d<4>(x) * d(y);
            static_assert(diff_op == res);
        }
    }

    return 0;
}
