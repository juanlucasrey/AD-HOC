#include <adhoc.hpp>
#include <backpropagator_tools.hpp>
#include <calc_tree.hpp>
#include <differential_operator.hpp>
#include <utility>
#include <utils/bivariate.hpp>

#include <utils/tuple.hpp>

#include <sort.hpp>

#include "call_price.hpp"
#include "type_name.hpp"

#include <gtest/gtest.h>

#include <random>
#include <tuple>

namespace adhoc4 {

TEST(BackPropagatorTools, FilteredCopy) {

    constexpr std::tuple<std::true_type, std::true_type, std::false_type,
                         std::true_type>
        temp;
    std::array<double, 4> arr{1, 2, 3, 4};
    std::array<double, 3> arr_to;

    detail::copy_filtered(arr, arr_to, temp);
    EXPECT_EQ(arr_to[0], 1);
    EXPECT_EQ(arr_to[1], 2);
    EXPECT_EQ(arr_to[2], 4);
}

TEST(BackPropagatorTools, FilteredCopyInverse) {

    {
        constexpr std::tuple<std::true_type, std::true_type, std::false_type,
                             std::true_type>
            temp;
        std::array<double, 4> arr{1, 2, 3, 4};
        std::array<double, 3> arr_to;

        detail::copy_filtered_inverted(arr, arr_to, temp);
        EXPECT_EQ(arr_to[0], 4);
        EXPECT_EQ(arr_to[1], 3);
        EXPECT_EQ(arr_to[2], 1);
    }

    {
        constexpr std::tuple<std::true_type, std::true_type> temp;
        std::array<double, 3> arr{1, 2, 3};
        std::array<double, 2> arr_to;

        detail::copy_filtered_inverted(arr, arr_to, temp);

        EXPECT_EQ(arr_to[0], 3);
        EXPECT_EQ(arr_to[1], 2);
    }

    {
        constexpr std::tuple<std::true_type, std::false_type> temp;
        std::array<double, 3> arr{1, 2, 3};
        std::array<double, 2> arr_to;
        arr_to.fill(0);

        detail::copy_filtered_inverted(arr, arr_to, temp);

        EXPECT_EQ(arr_to[0], 3);
        EXPECT_EQ(arr_to[1], 0);
    }
}

TEST(BackPropagatorTools, SortPair) {
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

TEST(BackPropagatorTools, CreateDiffOpFromIS2) {
    ADHOC(x);
    ADHOC(y);
    constexpr auto vars = std::make_tuple(d(x), d(y));

    constexpr auto is1 = std::index_sequence<4, 1>{};

    constexpr auto diff_op = create_differential_operator(vars, is1);

    constexpr auto res = d<4>(x) * d(y);
    static_assert(diff_op == res);
}

TEST(BackPropagatorTools, PowDiffOp) {
    ADHOC(x);
    ADHOC(y);
    // constexpr auto vars = std::make_tuple(x, y);

    // constexpr auto is1 = std::index_sequence<4, 1>{};

    // constexpr auto diff_op = create_differential_operator(vars, is1);

    constexpr auto res = d<4>(x) * d(y);
    constexpr auto res2 = power<0>(res);

    std::cout << "res2" << std::endl;
    std::cout << type_name2<decltype(res2)>() << std::endl;

    constexpr auto vars = std::make_tuple(d(x) * d(y), d(x), d(y));
    constexpr auto is = std::index_sequence<1, 2, 3>{};
    constexpr auto res3 = create_differential_operator(vars, is);

    std::cout << "res3" << std::endl;
    std::cout << type_name2<decltype(res3)>() << std::endl;

    {
        constexpr auto var1 = power<1>(d(x) * d(y));
        constexpr auto var2 = power<2>(d(x));
        constexpr auto var3 = var1 * var2;
        std::cout << "var3" << std::endl;
        std::cout << type_name2<decltype(var3)>() << std::endl;
    }

    {
        constexpr auto var1 = d(x) * d(y);
        constexpr auto var2 = d<2>(x);
        constexpr auto var3 = var1 * var2;
        // constexpr auto var3 = var1 + var2;
        std::cout << "var3" << std::endl;
        std::cout << type_name2<decltype(var3)>() << std::endl;
    }

    // constexpr auto var4 = (d(x) * d(y)) * d(x);
    // std::cout << "var4" << std::endl;
    // std::cout << type_name2<decltype(var4)>() << std::endl;

    // static_assert(diff_op == res);
}

} // namespace adhoc4
