#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_EXPAND_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_EXPAND_HPP

#include "../adhoc.hpp"
#include "../constants_type.hpp"
#include "../dependency.hpp"
#include "differential_operator.hpp"
#include "is_ordered.hpp"
#include "merge_ordered.hpp"
#include "order_differential_operator.hpp"
#include "select_root_derivatives.hpp"

namespace adhoc3 {

namespace detail {

template <typename T> struct is_constant_class : std::false_type {};

template <constants::ArgType N>
struct is_constant_class<constants::CD<N>> : std::true_type {};

template <class T>
constexpr bool is_constant_class_v = is_constant_class<T>::value;

template <class... Ids, std::size_t... Orders, std::size_t... Powers,
          class... CalculationNodes>
constexpr auto
expand_single(std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...> id,
              std::tuple<CalculationNodes...> nodes);

template <class... CalculationNodes>
constexpr auto expand_single(std::tuple<> /* id */,
                             std::tuple<CalculationNodes...> /* nodes */) {
    return std::tuple<std::tuple<>>{};
}

template <std::size_t N, class... Ids, std::size_t Order, std::size_t... Orders,
          std::size_t Power, std::size_t... Powers, class... CalculationNodes>
constexpr auto
expand_single(std::tuple<der2::p<Power, der2::d<Order, double_t<N>>>,
                         der2::p<Powers, der2::d<Orders, Ids>>...> /* id */,
              std::tuple<CalculationNodes...> /* nodes */) {
    return std::tuple<std::tuple<der2::p<Power, der2::d<Order, double_t<N>>>,
                                 der2::p<Powers, der2::d<Orders, Ids>>...>>{};
}

template <class Id1, class Id2, class... Ids, std::size_t Order,
          std::size_t... Orders, std::size_t Power, std::size_t... Powers,
          class... CalculationNodes>
constexpr auto
expand_single(std::tuple<der2::p<Power, der2::d<Order, add_t<Id1, Id2>>>,
                         der2::p<Powers, der2::d<Orders, Ids>>...> /* id */,
              std::tuple<CalculationNodes...> nodes) {
    if constexpr (is_constant_class_v<Id1>) {
        return expand_single(
            std::tuple<der2::p<Power, der2::d<Order, Id2>>,
                       der2::p<Powers, der2::d<Orders, Ids>>...>{});
    } else if constexpr (is_constant_class_v<Id2>) {
        return expand_single(
            std::tuple<der2::p<Power, der2::d<Order, Id1>>,
                       der2::p<Powers, der2::d<Orders, Ids>>...>{});
    } else {
        auto constexpr unordered_differential_operator_1 =
            std::tuple<der2::p<Power, der2::d<Order, Id1>>,
                       der2::p<Powers, der2::d<Orders, Ids>>...>{};

        auto constexpr ordered_differential_operator_1 =
            order_differential_operator(unordered_differential_operator_1,
                                        nodes);

        auto constexpr expanded_differential_operator_1 =
            expand_single(ordered_differential_operator_1, nodes);

        auto constexpr unordered_differential_operator_2 =
            std::tuple<der2::p<Power, der2::d<Order, Id2>>,
                       der2::p<Powers, der2::d<Orders, Ids>>...>{};

        auto constexpr ordered_differential_operator_2 =
            order_differential_operator(unordered_differential_operator_2,
                                        nodes);

        auto constexpr expanded_differential_operator_2 =
            expand_single(ordered_differential_operator_2, nodes);

        return merge_ordered(expanded_differential_operator_1,
                             expanded_differential_operator_2, nodes);
    }
}

template <class Id1, class Id2, class... Ids, std::size_t Order,
          std::size_t... Orders, std::size_t Power, std::size_t... Powers,
          class... CalculationNodes>
constexpr auto
expand_single(std::tuple<der2::p<Power, der2::d<Order, sub_t<Id1, Id2>>>,
                         der2::p<Powers, der2::d<Orders, Ids>>...> /* id */,
              std::tuple<CalculationNodes...> nodes) {
    // substraction and addition have the same expansion
    return expand_single(
        std::tuple<der2::p<Power, der2::d<Order, add_t<Id1, Id2>>>,
                   der2::p<Powers, der2::d<Orders, Ids>>...>{},
        nodes);
}

} // namespace detail

template <class... Ids, class... CalculationNodes>
constexpr auto expand(std::tuple<Ids...> in1,
                      std::tuple<CalculationNodes...> nodes);

template <class... CalculationNodes>
constexpr auto expand(std::tuple<> /* in1 */,
                      std::tuple<CalculationNodes...> /* nodes */) {
    return std::tuple<>{};
}

template <class Id, class... Ids, class... CalculationNodes>
constexpr auto expand(std::tuple<Id, Ids...> in,
                      std::tuple<CalculationNodes...> nodes) {
    static_assert(is_ordered(in, nodes));
    return merge_ordered(detail::expand_single(Id{}, nodes),
                         expand(std::tuple<Ids...>{}, nodes), nodes);
}

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_EXPAND_HPP
