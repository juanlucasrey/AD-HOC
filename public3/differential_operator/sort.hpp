#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_SORT_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_SORT_HPP

#include "../adhoc.hpp"
#include "../dependency.hpp"
#include "../tuple_utils.hpp"
#include "derivative_non_null.hpp"
#include "differential_operator.hpp"
#include <differential_operator/is_sorted.hpp>

#include "less_than.hpp"
#include <tuple>

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
    auto seq = std::make_index_sequence<numer_max_d>{};
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
    if constexpr (is_ordered(nodes, diff_operator)) {
        return diff_operator;
    } else {
        return std::apply(
            [diff_operator](auto &&...args) {
                return std::tuple_cat(
                    detail::create_single_tuple(args, diff_operator)...);
            },
            nodes);
    }
}

namespace detail {

template <class ResultOp, class CandidateOp, class Nodes>
constexpr auto return_largest(ResultOp res, CandidateOp candidate,
                              Nodes nodes) {
    if constexpr (less_than(nodes, candidate, res)) {
        return res;
    } else {
        return candidate;
    }
}

template <class ResultOp, class Nodes>
constexpr auto find_largest_aux(ResultOp res, std::tuple<> /* diff_operators */,
                                Nodes /* nodes */) {
    return res;
}

template <class ResultOp, class Op1, class... Ops, class Nodes>
constexpr auto find_largest_aux(ResultOp res,
                                std::tuple<Op1, Ops...> /* diff_operators */,
                                Nodes nodes) {
    return find_largest_aux(return_largest(res, Op1{}, nodes),
                            std::tuple<Ops...>{}, nodes);
}

template <class Op1, class... Ops, class Nodes>
constexpr auto find_largest(std::tuple<Op1, Ops...> /* diff_operators */,
                            Nodes nodes) {
    return find_largest_aux(Op1{}, std::tuple<Ops...>{}, nodes);
}

template <class Nodes, class Result>
constexpr auto
order_differential_operators_aux(std::tuple<> /* diff_operators */,
                                 Nodes /* nodes */, Result res) {
    return res;
}

template <class Ops, class Nodes, class Result>
constexpr auto order_differential_operators_aux(Ops diff_operators, Nodes nodes,
                                                Result res) {
    constexpr auto largest = find_largest(diff_operators, nodes);
    constexpr auto rem = remove(diff_operators, largest);
    return order_differential_operators_aux(
        rem, nodes, std::tuple_cat(res, std::make_tuple(largest)));
}

} // namespace detail

template <class Ops, class Nodes>
constexpr auto order_differential_operators(Ops diff_operators, Nodes nodes) {
    if constexpr (is_ordered(nodes, diff_operators)) {
        return diff_operators;
    } else {
        return detail::order_differential_operators_aux(diff_operators, nodes,
                                                        std::tuple<>{});
    }
}

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_SORT_HPP
