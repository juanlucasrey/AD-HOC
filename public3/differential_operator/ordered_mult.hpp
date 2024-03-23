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
          std::size_t... Powers2, class... Out, class Nodes>
constexpr auto
multiply_ordered_compare(std::tuple<der2::p<Power1, der2::d<Order1, Id1>>,
                                    der2::p<Powers1, der2::d<Orders1, Ids1>>...>
                             in1,
                         std::tuple<der2::p<Power2, der2::d<Order2, Id2>>,
                                    der2::p<Powers2, der2::d<Orders2, Ids2>>...>
                             in2,
                         std::tuple<Out...> out, Nodes nodes);

template <class DerOp1, class DerOp2, class Out, class Nodes>
constexpr auto multiply_ordered_empty(DerOp1 in1, DerOp2 in2, Out out,
                                      Nodes nodes) {
    if constexpr (std::tuple_size_v<DerOp2> == 0) {
        return std::tuple_cat(out, in1);
    } else if constexpr (std::tuple_size_v<DerOp1> == 0) {
        return std::tuple_cat(out, in2);
    } else {
        return multiply_ordered_compare(in1, in2, out, nodes);
    }
}

template <class Id1, class... Ids1, std::size_t Order1, std::size_t... Orders1,
          std::size_t Power1, std::size_t... Powers1, class Id2, class... Ids2,
          std::size_t Order2, std::size_t... Orders2, std::size_t Power2,
          std::size_t... Powers2, class... Out, class Nodes>
constexpr auto
multiply_ordered_compare(std::tuple<der2::p<Power1, der2::d<Order1, Id1>>,
                                    der2::p<Powers1, der2::d<Orders1, Ids1>>...>
                             in1,
                         std::tuple<der2::p<Power2, der2::d<Order2, Id2>>,
                                    der2::p<Powers2, der2::d<Orders2, Ids2>>...>
                             in2,
                         std::tuple<Out...> /* out */, Nodes nodes) {

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

template <class Nodes> constexpr auto multiply_ordered(Nodes /* nodes */) {
    return std::tuple<>{};
};

template <class Nodes, class DerOp>
constexpr auto multiply_ordered(Nodes nodes, DerOp in1) {
    static_assert(is_ordered(in1, nodes));
    return in1;
};

template <class Nodes, class DerOp1, class DerOp2>
constexpr auto multiply_ordered(Nodes nodes, DerOp1 in1, DerOp2 in2) {
    static_assert(is_ordered(in1, nodes));
    static_assert(is_ordered(in2, nodes));
    return detail::multiply_ordered_empty(in1, in2, std::tuple<>{}, nodes);
};

template <class Nodes, class DerOp1, class DerOp2, class... DerOps>
constexpr auto multiply_ordered(Nodes nodes, DerOp1 in1, DerOp2 in2,
                                DerOps... ins) {
    static_assert(is_ordered(in1, nodes));
    static_assert(is_ordered(in2, nodes));
    constexpr auto mult = multiply_ordered(nodes, in1, in2);
    return multiply_ordered(nodes, mult, ins...);
};

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_ORDERED_MULT_HPP
