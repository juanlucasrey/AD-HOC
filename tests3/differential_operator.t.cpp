#include <calc_order.hpp>
#include <init.hpp>
#include <tuple_utils.hpp>
#include <utils.hpp>

#include "type_name.hpp"

#include <gtest/gtest.h>

namespace adhoc3 {

template <class Id, std::size_t Order = 1> struct d {};
template <class Id, std::size_t Power = 1> struct p {};

template <class In, class... Ids, std::size_t... Powers, std::size_t... Orders>
constexpr auto create_single_tuple(In /* in */,
                                   std::tuple<p<d<Ids, Orders>, Powers>...> o) {
    if constexpr (has_type<In, Ids...>()) {
        constexpr auto idx = idx_type2<In, Ids...>();
        constexpr auto v = std::get<idx>(o);
        return std::tuple<decltype(v)>{};
    } else {
        return std::tuple<>{};
    }
};

template <class... Ids, std::size_t... Powers, std::size_t... Orders,
          class... CalculationNodes>
constexpr auto order_differential_operator(
    std::tuple<p<d<Ids, Orders>, Powers>...> diff_operator,
    std::tuple<CalculationNodes...> nodes) {
    return std::apply(
        [diff_operator](auto &&...args) {
            return std::tuple_cat(create_single_tuple(args, diff_operator)...);
        },
        nodes);
}

TEST(DifferentialOperator, dID) {
    auto [x, y] = Init<2>();
    auto res = (x * y) * (x * cos(y));
    constexpr auto co = calc_order_t<true>(res);

    std::tuple<p<d<decltype(x), 2>, 3>, p<d<decltype(y), 4>, 5>> der3{};
    std::tuple<p<d<decltype(y), 4>, 5>, p<d<decltype(x), 2>, 3>> der3_inv{};

    constexpr auto der3ordered = order_differential_operator(der3, co);
    constexpr auto der3_inv_ordered = order_differential_operator(der3_inv, co);

    static_assert(
        std::is_same<decltype(der3ordered), decltype(der3_inv_ordered)>::value);
}

} // namespace adhoc3