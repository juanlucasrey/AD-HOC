#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_ORDERED_MULT_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_ORDERED_MULT_HPP

#include "../adhoc.hpp"
#include "../dependency.hpp"
#include "differential_operator.hpp"
#include "is_ordered.hpp"
#include "less_than.hpp"

namespace adhoc3 {

namespace detail {

template <class Id1, class... Ids1, std::size_t Order1, std::size_t... Orders1,
          std::size_t Power1, std::size_t... Powers1, class Id2, class... Ids2,
          std::size_t Order2, std::size_t... Orders2, std::size_t Power2,
          std::size_t... Powers2, class... Out, class... CalculationNodes>
constexpr auto
multiply_ordered_compare(std::tuple<der2::p<Power1, der2::d<Order1, Id1>>,
                                    der2::p<Powers1, der2::d<Orders1, Ids1>>...>
                             in1,
                         std::tuple<der2::p<Power2, der2::d<Order2, Id2>>,
                                    der2::p<Powers2, der2::d<Orders2, Ids2>>...>
                             in2,
                         std::tuple<Out...> out,
                         std::tuple<CalculationNodes...> nodes);

template <class... Ids1, std::size_t... Orders1, std::size_t... Powers1,
          class... Ids2, std::size_t... Orders2, std::size_t... Powers2,
          class... Out, class... CalculationNodes>
constexpr auto multiply_ordered_empty(
    std::tuple<der2::p<Powers1, der2::d<Orders1, Ids1>>...> in1,
    std::tuple<der2::p<Powers2, der2::d<Orders2, Ids2>>...> in2,
    std::tuple<Out...> out, std::tuple<CalculationNodes...> nodes) {

    if constexpr (sizeof...(Ids2) == 0) {
        return std::tuple_cat(out, in1);
    } else if constexpr (sizeof...(Ids1) == 0) {
        return std::tuple_cat(out, in2);
    } else {
        return multiply_ordered_compare(in1, in2, out, nodes);
    }
}

template <class Id1, class... Ids1, std::size_t Order1, std::size_t... Orders1,
          std::size_t Power1, std::size_t... Powers1, class Id2, class... Ids2,
          std::size_t Order2, std::size_t... Orders2, std::size_t Power2,
          std::size_t... Powers2, class... Out, class... CalculationNodes>
constexpr auto
multiply_ordered_compare(std::tuple<der2::p<Power1, der2::d<Order1, Id1>>,
                                    der2::p<Powers1, der2::d<Orders1, Ids1>>...>
                             in1,
                         std::tuple<der2::p<Power2, der2::d<Order2, Id2>>,
                                    der2::p<Powers2, der2::d<Orders2, Ids2>>...>
                             in2,
                         std::tuple<Out...> /* out */,
                         std::tuple<CalculationNodes...> nodes) {

    if constexpr (less_than(der2::p<Power1, der2::d<Order1, Id1>>{},
                            der2::p<Power2, der2::d<Order2, Id2>>{}, nodes)) {

        return multiply_ordered_empty(
            in1, std::tuple<der2::p<Powers2, der2::d<Orders2, Ids2>>...>{},
            std::tuple<Out..., der2::p<Power2, der2::d<Order2, Id2>>>{}, nodes);

    } else if constexpr (less_than(der2::p<Power2, der2::d<Order2, Id2>>{},
                                   der2::p<Power1, der2::d<Order1, Id1>>{},
                                   nodes)) {

        return multiply_ordered_empty(
            std::tuple<der2::p<Powers1, der2::d<Orders1, Ids1>>...>{}, in2,
            std::tuple<Out..., der2::p<Power1, der2::d<Order1, Id1>>>{}, nodes);

    } else {

        static_assert(
            std::is_same_v<der2::d<Order1, Id1>, der2::d<Order2, Id2>>);
        return multiply_ordered_empty(
            std::tuple<der2::p<Powers1, der2::d<Orders1, Ids1>>...>{},
            std::tuple<der2::p<Powers2, der2::d<Orders2, Ids2>>...>{},
            std::tuple<Out...,
                       der2::p<Power1 + Power2, der2::d<Order1, Id1>>>{},
            nodes);
    }
}

} // namespace detail

template <class... Ids1, std::size_t... Orders1, std::size_t... Powers1,
          class... Ids2, std::size_t... Orders2, std::size_t... Powers2,
          class... CalculationNodes>
constexpr auto
multiply_ordered(std::tuple<der2::p<Powers1, der2::d<Orders1, Ids1>>...> in1,
                 std::tuple<der2::p<Powers2, der2::d<Orders2, Ids2>>...> in2,
                 std::tuple<CalculationNodes...> nodes) {
    static_assert(is_ordered(in1, nodes));
    static_assert(is_ordered(in2, nodes));
    return detail::multiply_ordered_empty(in1, in2, std::tuple<>{}, nodes);
};

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_ORDERED_MULT_HPP
