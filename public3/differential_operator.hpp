#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_HPP

// #include "adhoc.hpp"
#include "dependency.hpp"

namespace adhoc3 {

namespace der2 {

template <class Id, std::size_t Order = 1> struct d {};
template <class Id, std::size_t Power = 1> struct p {};

} // namespace der2

namespace detail {

template <class IdInput, std::size_t... PowersNode, std::size_t... OrdersNode,
          class... IdsNode>
constexpr auto der_non_null_aux(
    std::tuple<der2::p<der2::d<IdsNode, OrdersNode>, PowersNode>...>)
    -> std::size_t {
    return (detail::sum_no_overflow<order<IdsNode, IdInput>::call()...>());
}

} // namespace detail

template <std::size_t... PowersNode, std::size_t... OrdersNode,
          class... IdsNode, std::size_t... PowersInput, class... IdsInput>
constexpr auto
der_non_null(std::tuple<der2::p<der2::d<IdsNode, OrdersNode>, PowersNode>...> n,
             std::tuple<der2::p<der2::d<IdsInput, 1>, PowersInput>...>)
    -> bool {
    return ((PowersInput <= detail::der_non_null_aux<IdsInput>(n)) && ...);
}

namespace detail {

template <class In, class... Ids, std::size_t... Powers, std::size_t... Orders>
constexpr auto
create_single_tuple(In /* in */,
                    std::tuple<der2::p<der2::d<Ids, Orders>, Powers>...> o) {
    if constexpr (has_type<In, Ids...>()) {
        constexpr auto idx = idx_type2<In, Ids...>();
        constexpr auto v = std::get<idx>(o);
        return std::tuple<decltype(v)>{};
    } else {
        return std::tuple<>{};
    }
};

} // namespace detail

template <class... Ids, std::size_t... Powers, std::size_t... Orders,
          class... CalculationNodes>
constexpr auto order_differential_operator(
    std::tuple<der2::p<der2::d<Ids, Orders>, Powers>...> diff_operator,
    std::tuple<CalculationNodes...> nodes) {
    return std::apply(
        [diff_operator](auto &&...args) {
            return std::tuple_cat(
                detail::create_single_tuple(args, diff_operator)...);
        },
        nodes);
}

template <std::size_t Order, class Id> constexpr auto diff_op(Id /* id */) {
    return std::tuple<der2::p<der2::d<Id, 1>, Order>>{};
}

template <std::size_t... Orders, class... Ids>
constexpr auto diff_op(std::tuple<der2::p<der2::d<Ids, 1>, Orders>>... id) {
    return std::tuple_cat(id...);
}

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_HPP