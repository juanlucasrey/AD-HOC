#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_INTEGER_PARTITION_UNIVARIATE_FUNCTION_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_INTEGER_PARTITION_UNIVARIATE_FUNCTION_HPP

#include "../partition/partition_function.hpp"
#include "differential_operator.hpp"

namespace adhoc3 {

namespace detail {

template <std::size_t... Ns, class... Ts>
constexpr auto pop_back_impl(std::index_sequence<Ns...>, std::tuple<Ts...> t) {
    return std::make_tuple(std::get<Ns>(t)...);
}

template <class... Ts> constexpr auto pop_back(std::tuple<Ts...> t) {
    if constexpr (sizeof...(Ts) < 2) {
        return std::tuple<>{};
    } else {
        return pop_back_impl(std::make_index_sequence<sizeof...(Ts) - 1>(), t);
    }
}

template <class... Ts> constexpr auto back(std::tuple<Ts...> t) {
    if constexpr (sizeof...(Ts) == 0) {
        return std::tuple<>{};
    } else {
        return std::get<sizeof...(Ts) - 1u>(t);
    }
}

template <std::size_t Rem, class Id, std::size_t... Orders,
          std::size_t... Powers, std::size_t Order, std::size_t Power>
constexpr auto get_next_partition_sep(
    std::tuple<der2::p<Powers, der2::d<Orders, Id>>...> /* id */,
    der2::p<Power, der2::d<Order, Id>> /* id */) {
    constexpr std::size_t remainder = Rem + Order;
    constexpr std::size_t lower_order = Order - 1;
    constexpr std::size_t div = remainder / lower_order;
    constexpr std::size_t new_remainder = remainder - lower_order * div;
    if constexpr (Power > 1) {
        if constexpr (new_remainder) {
            return std::tuple<der2::p<Powers, der2::d<Orders, Id>>...,
                              der2::p<Power - 1, der2::d<Order, Id>>,
                              der2::p<div, der2::d<lower_order, Id>>,
                              der2::p<1, der2::d<new_remainder, Id>>>{};
        } else {
            return std::tuple<der2::p<Powers, der2::d<Orders, Id>>...,
                              der2::p<Power - 1, der2::d<Order, Id>>,
                              der2::p<div, der2::d<lower_order, Id>>>{};
        }

    } else {
        if constexpr (new_remainder) {
            return std::tuple<der2::p<Powers, der2::d<Orders, Id>>...,
                              der2::p<div, der2::d<lower_order, Id>>,
                              der2::p<1, der2::d<new_remainder, Id>>>{};
        } else {
            return std::tuple<der2::p<Powers, der2::d<Orders, Id>>...,
                              der2::p<div, der2::d<lower_order, Id>>>{};
        }
    }
}

template <std::size_t Rem, class Id, std::size_t... Orders,
          std::size_t... Powers>
constexpr auto
get_next_partition_aux(std::tuple<der2::p<Powers, der2::d<Orders, Id>>...> id) {
    static_assert(sizeof...(Orders));
    return get_next_partition_sep<Rem>(pop_back(id), back(id));
}

template <class Id, std::size_t... Orders, std::size_t... Powers>
constexpr auto
get_next_partition(std::tuple<der2::p<Powers, der2::d<Orders, Id>>...> id) {

    constexpr auto last_order =
        std::get<sizeof...(Orders) - 1>(std::forward_as_tuple(Orders...));

    constexpr auto rem =
        last_order == 1
            ? std::get<sizeof...(Powers) - 1>(std::forward_as_tuple(Powers...))
            : 0;

    if constexpr (last_order == 1) {
        constexpr auto id_no_1 = pop_back(id);
        return get_next_partition_aux<rem>(pop_back(id));
    } else {
        return get_next_partition_aux<rem>(id);
    }
}

template <std::size_t Idx, std::size_t End, class Last, class... Out>
constexpr auto
integer_partition_univariate_function(Last last, std::tuple<Out...> /* out */) {
    if constexpr (Idx == End) {
        return std::tuple<Out..., Last>{};
    } else {
        auto constexpr next = get_next_partition(last);
        return integer_partition_univariate_function<Idx + 1, End>(
            next, std::tuple<Out..., Last>{});
    }
};

} // namespace detail

template <class Id, std::size_t Order>
constexpr auto integer_partition_univariate_function() {
    return detail::integer_partition_univariate_function<1, partition_function(
                                                                Order)>(
        std::tuple<der2::p<1, der2::d<Order, Id>>>{}, std::tuple<>{});
}

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_INTEGER_PARTITION_UNIVARIATE_FUNCTION_HPP
