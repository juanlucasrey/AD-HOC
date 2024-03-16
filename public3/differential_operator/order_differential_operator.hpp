#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_ORDER_DIFFERENTIAL_OPERATOR_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_ORDER_DIFFERENTIAL_OPERATOR_HPP

#include "../adhoc.hpp"
#include "../dependency.hpp"
#include "derivative_non_null.hpp"
#include "differential_operator.hpp"

namespace adhoc3 {

namespace detail {

template <class In, class Id, std::size_t Power, std::size_t Order>
constexpr auto
create_single_tuple2(In /* in */, der2::p<Power, der2::d<Order, Id>> /* o */) {
    if constexpr (std::is_same_v<In, Id>) {
        return std::tuple<der2::p<Power, der2::d<Order, Id>>>{};
    } else {
        return std::tuple<>{};
    }
};

template <std::size_t I, class Id, std::size_t Power, std::size_t Order>
constexpr auto order_by_order4(der2::p<Power, der2::d<Order, Id>> /* o */) {
    if constexpr (I == Order) {
        return std::tuple<der2::p<Power, der2::d<Order, Id>>>{};
    } else {
        return std::tuple<>{};
    }
}

template <std::size_t I, class Id, std::size_t... Powers, std::size_t... Orders>
constexpr auto
order_by_order3(std::tuple<der2::p<Powers, der2::d<Orders, Id>>...> o) {

    return std::apply(
        [](auto &&...args) {
            return std::tuple_cat(order_by_order4<I>(args)...);
        },
        o);
}

template <class Id, std::size_t... Powers, std::size_t... Orders,
          std::size_t... I>
constexpr auto
order_by_order2(std::tuple<der2::p<Powers, der2::d<Orders, Id>>...> o,
                std::index_sequence<I...> /* i */) {
    return std::tuple_cat(order_by_order3<I + 1>(o)...);
}

template <class In> constexpr auto order_by_order(In o);

template <class Id, std::size_t... Powers, std::size_t... Orders>
constexpr auto
order_by_order(std::tuple<der2::p<Powers, der2::d<Orders, Id>>...> o) {
    auto constexpr numer_max_d = detail::max<Orders...>();
    auto seq = std::make_integer_sequence<std::size_t, numer_max_d>{};
    return order_by_order2(o, seq);
};

template <> constexpr auto order_by_order(std::tuple<> /* o */) {
    return std::tuple<>{};
};

template <class In, class... Ids, std::size_t... Powers, std::size_t... Orders>
constexpr auto
create_single_tuple(In in,
                    std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...> o) {
    auto tuple_to_order = std::apply(
        [in](auto &&...args) {
            return std::tuple_cat(detail::create_single_tuple2(in, args)...);
        },
        o);

    return order_by_order(tuple_to_order);
};

} // namespace detail

template <class... Ids, std::size_t... Orders, std::size_t... Powers,
          class... CalculationNodes>
constexpr auto order_differential_operator(
    std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...> diff_operator,
    std::tuple<CalculationNodes...> nodes) {
    return std::apply(
        [diff_operator](auto &&...args) {
            return std::tuple_cat(
                detail::create_single_tuple(args, diff_operator)...);
        },
        nodes);
}

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_ORDER_DIFFERENTIAL_OPERATOR_HPP
