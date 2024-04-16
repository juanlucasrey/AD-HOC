#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_ORDER_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_ORDER_HPP

#include <differential_operator/differential_operator.hpp>

namespace adhoc3 {

constexpr auto max(auto const &value, auto const &...args) {
    if constexpr (sizeof...(args) == 0U) {
        return value;
    } else {
        const auto max_val = max(args...);
        return value > max_val ? value : max_val;
    }
}

template <class... Ids, std::size_t... Orders, std::size_t... Powers>
constexpr auto
max_order(std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...> /* id1 */) {
    return max(Orders...);
}

template <class... Ops>
constexpr auto max_orders(std::tuple<Ops...> /* id1 */) {
    return max(max_order(Ops{})...);
}

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_ORDER_HPP
