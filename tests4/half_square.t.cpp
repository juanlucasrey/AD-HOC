#include <adhoc.hpp>
#include <backpropagator.hpp>
#include <calc_tree.hpp>
#include <differential_operator.hpp>

#include <utils/tuple.hpp>

#include <sort.hpp>

#include "call_price.hpp"
#include "type_name.hpp"

#include <gtest/gtest.h>

#include <tuple>

namespace adhoc4 {

namespace {

constexpr auto half_square_size(std::size_t max) -> std::size_t {
    return (max * (max + 1)) / 2;
}

constexpr auto half_square_index(std::size_t y, std::size_t x, std::size_t max)
    -> std::size_t {
    return x + max * y - (y * (y - 1)) / 2;
}

template <std::size_t, class Type> using return_type = Type;

template <std::size_t Size, std::size_t... I>
constexpr auto init_flags_aux(std::index_sequence<I...> /* i */) {
    return std::make_tuple(return_type<I, std::false_type>{}...);
}

template <std::size_t Size> constexpr auto init_flags() {
    return init_flags_aux<Size>(std::make_index_sequence<Size>{});
}

template <std::size_t Loc, class TupleBool, std::size_t... First,
          std::size_t... Last>
constexpr auto switch_flag_aux(TupleBool tb,
                               std::index_sequence<First...> /* i */,
                               std::index_sequence<Last...> /* i */) {
    return std::tuple_cat(
        std::make_tuple(std::get<First>(tb))...,
        std::make_tuple(
            std::conditional_t<std::tuple_element_t<Loc, TupleBool>::value,
                               std::false_type, std::true_type>{}),
        std::make_tuple(std::get<Last + Loc + 1>(tb))...);
}

template <std::size_t Loc, class TupleBool>
constexpr auto switch_flag(TupleBool tb) {
    constexpr auto total_size = std::tuple_size_v<TupleBool>;
    static_assert(total_size > Loc);
    return switch_flag_aux<Loc>(
        tb, std::make_index_sequence<Loc>{},
        std::make_index_sequence<total_size - Loc - 1>{});
}

template <std::size_t Row, std::size_t Col, std::size_t Side, std::size_t Size,
          class TupleBool>
auto set(std::array<double, Size> & /* vals */, TupleBool tb) {
    if constexpr (std::tuple_element_t<half_square_index(Row, Col, Side),
                                       TupleBool>::value) {
        return tb;
    } else {
        return tb;
    }
}

template <std::size_t Loc, std::size_t End, std::size_t Size>
auto fill_side1(std::array<double, Size> &vals, double val) {
    if constexpr (Loc < End) {
        vals[half_square_index(0, Loc, End)] =
            vals[half_square_index(0, Loc - 1, End)] * val;
        fill_side1<Loc + 1, End>(vals, val);
    }
}

template <std::size_t Loc, std::size_t End, std::size_t Size>
auto fill_side2(std::array<double, Size> &vals, double val) {
    if constexpr (Loc < End) {
        vals[half_square_index(Loc, 0, End)] =
            vals[half_square_index(Loc - 1, 0, End)] * val;
        fill_side2<Loc + 1, End>(vals, val);
    }
}

template <std::size_t Loc, std::size_t Size>
auto fill_inside(std::array<double, Size> & /* vals */) {}

template <std::size_t Order> auto mult_coeffs(double a, double b) {
    static_assert(Order > 0);

    std::array<double, half_square_size(Order + 1)> vals;
    vals[0] = 1.;

    vals[half_square_index(0, 1, Order + 1)] = a;
    fill_side1<2, Order + 1>(vals, a);

    vals[half_square_index(1, 0, Order + 1)] = b;
    fill_side2<2, Order + 1>(vals, b);

    fill_inside<Order + 1>(vals);

    return vals;
}

// template <std::size_t Loc, std::size_t Size>
// auto powers_fill(std::array<double, Size> &vals, double val) {
//     if constexpr (Loc < Size) {
//         vals[Loc] = vals[Loc - 1] * val;
//         powers_fill<Loc + 1>(vals, val);
//     }
// }

// template <std::size_t Order> auto powers(double a) {
//     static_assert(Order > 0);
//     std::array<double, Order> vals;
//     vals[0] = a;
//     powers_fill<1>(vals, a);
//     return vals;
// }

} // namespace

TEST(BackPropagator, TriangleIndex) {
    static_assert(half_square_size(5) == 15);
    static_assert(half_square_index(0, 0, 5) == 0);
    static_assert(half_square_index(0, 1, 5) == 1);
    static_assert(half_square_index(0, 2, 5) == 2);
    static_assert(half_square_index(0, 3, 5) == 3);
    static_assert(half_square_index(0, 4, 5) == 4);
    static_assert(half_square_index(1, 0, 5) == 5);
    static_assert(half_square_index(1, 1, 5) == 6);
    static_assert(half_square_index(1, 2, 5) == 7);
    static_assert(half_square_index(1, 3, 5) == 8);
    static_assert(half_square_index(2, 0, 5) == 9);
    static_assert(half_square_index(2, 1, 5) == 10);
    static_assert(half_square_index(2, 2, 5) == 11);
    static_assert(half_square_index(3, 0, 5) == 12);
    static_assert(half_square_index(3, 1, 5) == 13);
    static_assert(half_square_index(4, 0, 5) == 14);

    double a = 0.1;
    auto const res = detail::powers<5>(a);
    for (std::size_t i = 0; i < res.size(); i++) {
        std::cout << res[i] << ", ";
    }

    // double a = 0.1;
    // double b = 3.14;
    // auto res = mult_coeffs<4>(a, b);

    // for (std::size_t i = 0; i < res.size(); i++) {
    //     std::cout << res[i] << ", ";
    // }
    // std::cout << std::endl;

    // constexpr auto size = half_square_size(5);
    // constexpr auto flags = init_flags<size>();
    // constexpr auto flags1 = switch_flag<half_square_index(0, 0, 5)>(flags);
    // constexpr auto flags2 = switch_flag<half_square_index(1, 0, 5)>(flags1);
    // constexpr auto flags3 = switch_flag<half_square_index(1, 0, 5)>(flags2);

    // std::array<double, size> vals;
    // vals[half_square_index(0, 0, 5)] = 1.;
    // vals[half_square_index(1, 0, 5)] = a;
    // vals[half_square_index(0, 1, 5)] = b;

    // auto flags4 = set<0, 0, 5>(vals, flags3);

    // auto new_flags = switch_flag<half_square_index(0, 0, size)>(flags);
    // {
    //     auto flags = init_flags<5>();
    //     std::cout << type_name2<decltype(flags)>() << std::endl;

    //     auto new_flags = switch_flag<3>(flags);
    //     std::cout << type_name2<decltype(new_flags)>() << std::endl;

    //     auto new_flags2 = switch_flag<2>(new_flags);
    //     std::cout << type_name2<decltype(new_flags2)>() << std::endl;

    //     auto new_flags3 = switch_flag<0>(new_flags2);
    //     std::cout << type_name2<decltype(new_flags3)>() << std::endl;

    //     auto new_flags4 = switch_flag<4>(new_flags3);
    //     std::cout << type_name2<decltype(new_flags4)>() << std::endl;
    // }
}

} // namespace adhoc4
