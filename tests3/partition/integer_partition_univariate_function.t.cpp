#include <differential_operator/differential_operator.hpp>
#include <differential_operator/integer_partition_univariate_function.hpp>
#include <init.hpp>

#include <gtest/gtest.h>

namespace adhoc3 {

TEST(DifferentialOperator, ExpandUnivariateDifferentialOperator1) {
    auto [x1] = Init<1>();
    constexpr auto expansion =
        integer_partition_univariate_function<decltype(x1), 1>();

    constexpr auto result =
        std::tuple<std::tuple<der2::p<1, der2::d<1, decltype(x1)>>>>{};

    static_assert(std::is_same<decltype(expansion), decltype(result)>::value);
}

TEST(DifferentialOperator, ExpandUnivariateDifferentialOperator2) {
    auto [x1] = Init<1>();
    constexpr auto expansion =
        integer_partition_univariate_function<decltype(x1), 2>();

    constexpr auto result =
        std::tuple<std::tuple<der2::p<1, der2::d<2, decltype(x1)>>>,
                   std::tuple<der2::p<2, der2::d<1, decltype(x1)>>>>{};

    static_assert(std::is_same<decltype(expansion), decltype(result)>::value);
}

TEST(DifferentialOperator, ExpandUnivariateDifferentialOperator3) {
    auto [x1] = Init<1>();
    constexpr auto expansion =
        integer_partition_univariate_function<decltype(x1), 3>();

    constexpr auto result =
        std::tuple<std::tuple<der2::p<1, der2::d<3, decltype(x1)>>>,
                   std::tuple<der2::p<1, der2::d<2, decltype(x1)>>,
                              der2::p<1, der2::d<1, decltype(x1)>>>,
                   std::tuple<der2::p<3, der2::d<1, decltype(x1)>>>>{};

    static_assert(std::is_same<decltype(expansion), decltype(result)>::value);
}
TEST(DifferentialOperator, ExpandUnivariateDifferentialOperator7) {
    auto [x1] = Init<1>();

    //  {0, 0, 0, 0, 0, 0, 1},
    //  {1, 0, 0, 0, 0, 1, 0},
    constexpr std::tuple<der2::p<1, der2::d<7, decltype(x1)>>> in1;
    constexpr auto out1 = detail::get_next_partition(in1);
    constexpr std::tuple<der2::p<1, der2::d<6, decltype(x1)>>,
                         der2::p<1, der2::d<1, decltype(x1)>>>
        result1;
    static_assert(std::is_same<decltype(result1), decltype(out1)>::value);

    //  {1, 0, 0, 0, 0, 1, 0},
    //  {0, 1, 0, 0, 1, 0, 0},
    std::tuple<der2::p<1, der2::d<6, decltype(x1)>>,
               der2::p<1, der2::d<1, decltype(x1)>>>
        in2;
    constexpr auto out2 = detail::get_next_partition(in2);
    constexpr std::tuple<der2::p<1, der2::d<5, decltype(x1)>>,
                         der2::p<1, der2::d<2, decltype(x1)>>>
        result2;
    static_assert(std::is_same<decltype(result2), decltype(out2)>::value);

    //  {0, 1, 0, 0, 1, 0, 0},
    //  {2, 0, 0, 0, 1, 0, 0},
    std::tuple<der2::p<1, der2::d<5, decltype(x1)>>,
               der2::p<1, der2::d<2, decltype(x1)>>>
        in3;
    constexpr auto out3 = detail::get_next_partition(in3);
    constexpr std::tuple<der2::p<1, der2::d<5, decltype(x1)>>,
                         der2::p<2, der2::d<1, decltype(x1)>>>
        result3;

    static_assert(std::is_same<decltype(result3), decltype(out3)>::value);

    //  {2, 0, 0, 0, 1, 0, 0},
    //  {0, 0, 1, 1, 0, 0, 0},
    std::tuple<der2::p<1, der2::d<5, decltype(x1)>>,
               der2::p<2, der2::d<1, decltype(x1)>>>
        in4;
    constexpr auto out4 = detail::get_next_partition(in4);
    constexpr std::tuple<der2::p<1, der2::d<4, decltype(x1)>>,
                         der2::p<1, der2::d<3, decltype(x1)>>>
        result4;
    static_assert(std::is_same<decltype(result4), decltype(out4)>::value);

    //  {0, 0, 1, 1, 0, 0, 0},
    //  {1, 1, 0, 1, 0, 0, 0},
    std::tuple<der2::p<1, der2::d<4, decltype(x1)>>,
               der2::p<1, der2::d<3, decltype(x1)>>>
        in5;
    constexpr auto out5 = detail::get_next_partition(in5);
    constexpr std::tuple<der2::p<1, der2::d<4, decltype(x1)>>,
                         der2::p<1, der2::d<2, decltype(x1)>>,
                         der2::p<1, der2::d<1, decltype(x1)>>>
        result5;
    static_assert(std::is_same<decltype(result5), decltype(out5)>::value);

    //  {1, 1, 0, 1, 0, 0, 0},
    //  {3, 0, 0, 1, 0, 0, 0},
    std::tuple<der2::p<1, der2::d<4, decltype(x1)>>,
               der2::p<1, der2::d<2, decltype(x1)>>,
               der2::p<1, der2::d<1, decltype(x1)>>>
        in6;
    constexpr auto out6 = detail::get_next_partition(in6);
    constexpr std::tuple<der2::p<1, der2::d<4, decltype(x1)>>,
                         der2::p<3, der2::d<1, decltype(x1)>>>
        result6;

    static_assert(std::is_same<decltype(result6), decltype(out6)>::value);

    //  {3, 0, 0, 1, 0, 0, 0},
    //  {1, 0, 2, 0, 0, 0, 0},
    std::tuple<der2::p<1, der2::d<4, decltype(x1)>>,
               der2::p<3, der2::d<1, decltype(x1)>>>
        in7;
    constexpr auto out7 = detail::get_next_partition(in7);
    constexpr std::tuple<der2::p<2, der2::d<3, decltype(x1)>>,
                         der2::p<1, der2::d<1, decltype(x1)>>>
        result7;
    static_assert(std::is_same<decltype(result7), decltype(out7)>::value);

    //  {1, 0, 2, 0, 0, 0, 0},
    //  {0, 2, 1, 0, 0, 0, 0},
    std::tuple<der2::p<2, der2::d<3, decltype(x1)>>,
               der2::p<1, der2::d<1, decltype(x1)>>>
        in8;
    constexpr auto out8 = detail::get_next_partition(in8);
    constexpr std::tuple<der2::p<1, der2::d<3, decltype(x1)>>,
                         der2::p<2, der2::d<2, decltype(x1)>>>
        result8;
    static_assert(std::is_same<decltype(result8), decltype(out8)>::value);

    //  {0, 2, 1, 0, 0, 0, 0},
    //  {2, 1, 1, 0, 0, 0, 0},
    std::tuple<der2::p<1, der2::d<3, decltype(x1)>>,
               der2::p<2, der2::d<2, decltype(x1)>>>
        in9;
    constexpr auto out9 = detail::get_next_partition(in9);
    constexpr std::tuple<der2::p<1, der2::d<3, decltype(x1)>>,
                         der2::p<1, der2::d<2, decltype(x1)>>,
                         der2::p<2, der2::d<1, decltype(x1)>>>
        result9;
    static_assert(std::is_same<decltype(result9), decltype(out9)>::value);

    //  {2, 1, 1, 0, 0, 0, 0},
    //  {4, 0, 1, 0, 0, 0, 0},
    std::tuple<der2::p<1, der2::d<3, decltype(x1)>>,
               der2::p<1, der2::d<2, decltype(x1)>>,
               der2::p<2, der2::d<1, decltype(x1)>>>
        in10;
    constexpr auto out10 = detail::get_next_partition(in10);
    constexpr std::tuple<der2::p<1, der2::d<3, decltype(x1)>>,
                         der2::p<4, der2::d<1, decltype(x1)>>>
        result10;
    static_assert(std::is_same<decltype(result10), decltype(out10)>::value);

    //  {4, 0, 1, 0, 0, 0, 0},
    //  {1, 3, 0, 0, 0, 0, 0},
    std::tuple<der2::p<1, der2::d<3, decltype(x1)>>,
               der2::p<4, der2::d<1, decltype(x1)>>>
        in11;
    constexpr auto out11 = detail::get_next_partition(in11);
    constexpr std::tuple<der2::p<3, der2::d<2, decltype(x1)>>,
                         der2::p<1, der2::d<1, decltype(x1)>>>
        result11;
    static_assert(std::is_same<decltype(result11), decltype(out11)>::value);

    //  {1, 3, 0, 0, 0, 0, 0},
    //  {3, 2, 0, 0, 0, 0, 0},
    std::tuple<der2::p<3, der2::d<2, decltype(x1)>>,
               der2::p<1, der2::d<1, decltype(x1)>>>
        in12;
    constexpr auto out12 = detail::get_next_partition(in12);
    constexpr std::tuple<der2::p<2, der2::d<2, decltype(x1)>>,
                         der2::p<3, der2::d<1, decltype(x1)>>>
        result12;
    static_assert(std::is_same<decltype(result12), decltype(out12)>::value);

    //  {3, 2, 0, 0, 0, 0, 0},
    //  {5, 1, 0, 0, 0, 0, 0},
    std::tuple<der2::p<2, der2::d<2, decltype(x1)>>,
               der2::p<3, der2::d<1, decltype(x1)>>>
        in13;
    constexpr auto out13 = detail::get_next_partition(in13);
    constexpr std::tuple<der2::p<1, der2::d<2, decltype(x1)>>,
                         der2::p<5, der2::d<1, decltype(x1)>>>
        result13;
    static_assert(std::is_same<decltype(result13), decltype(out13)>::value);

    //  {5, 1, 0, 0, 0, 0, 0},
    //  {7, 0, 0, 0, 0, 0, 0}
    std::tuple<der2::p<1, der2::d<2, decltype(x1)>>,
               der2::p<5, der2::d<1, decltype(x1)>>>
        in14;
    constexpr auto out14 = detail::get_next_partition(in14);
    constexpr std::tuple<der2::p<7, der2::d<1, decltype(x1)>>> result14;
    static_assert(std::is_same<decltype(result14), decltype(out14)>::value);

    constexpr auto expansion =
        integer_partition_univariate_function<decltype(x1), 7>();

    static_assert(
        std::is_same<decltype(in1), std::remove_reference<decltype(std::get<0>(
                                        expansion))>::type>::value);

    static_assert(std::is_same<decltype(result1),
                               std::remove_reference<decltype(std::get<1>(
                                   expansion))>::type>::value);

    static_assert(std::is_same<decltype(result2),
                               std::remove_reference<decltype(std::get<2>(
                                   expansion))>::type>::value);

    static_assert(std::is_same<decltype(result3),
                               std::remove_reference<decltype(std::get<3>(
                                   expansion))>::type>::value);

    static_assert(std::is_same<decltype(result4),
                               std::remove_reference<decltype(std::get<4>(
                                   expansion))>::type>::value);

    static_assert(std::is_same<decltype(result5),
                               std::remove_reference<decltype(std::get<5>(
                                   expansion))>::type>::value);

    static_assert(std::is_same<decltype(result6),
                               std::remove_reference<decltype(std::get<6>(
                                   expansion))>::type>::value);

    static_assert(std::is_same<decltype(result7),
                               std::remove_reference<decltype(std::get<7>(
                                   expansion))>::type>::value);

    static_assert(std::is_same<decltype(result8),
                               std::remove_reference<decltype(std::get<8>(
                                   expansion))>::type>::value);

    static_assert(std::is_same<decltype(result9),
                               std::remove_reference<decltype(std::get<9>(
                                   expansion))>::type>::value);

    static_assert(std::is_same<decltype(result10),
                               std::remove_reference<decltype(std::get<10>(
                                   expansion))>::type>::value);

    static_assert(std::is_same<decltype(result11),
                               std::remove_reference<decltype(std::get<11>(
                                   expansion))>::type>::value);

    static_assert(std::is_same<decltype(result12),
                               std::remove_reference<decltype(std::get<12>(
                                   expansion))>::type>::value);

    static_assert(std::is_same<decltype(result13),
                               std::remove_reference<decltype(std::get<13>(
                                   expansion))>::type>::value);
}

} // namespace adhoc3
