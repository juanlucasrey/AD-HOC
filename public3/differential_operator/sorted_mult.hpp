#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_SORTED_MULT_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_SORTED_MULT_HPP

#include "../adhoc.hpp"
#include "../dependency.hpp"
#include "differential_operator.hpp"
#include "less_than.hpp"
#include <differential_operator/is_sorted.hpp>

namespace adhoc3 {

namespace detail {

template <class Nodes, class Id1, class... Ids1, std::size_t Order1,
          std::size_t... Orders1, std::size_t Power1, std::size_t... Powers1,
          class Id2, class... Ids2, std::size_t Order2, std::size_t... Orders2,
          std::size_t Power2, std::size_t... Powers2, class... Out>
constexpr auto
multiply_ordered_compare(Nodes nodes,
                         std::tuple<der2::p<Power1, der2::d<Order1, Id1>>,
                                    der2::p<Powers1, der2::d<Orders1, Ids1>>...>
                             in1,
                         std::tuple<der2::p<Power2, der2::d<Order2, Id2>>,
                                    der2::p<Powers2, der2::d<Orders2, Ids2>>...>
                             in2,
                         std::tuple<Out...> out);

template <class Nodes, class DerOp1, class DerOp2, class Out>
constexpr auto multiply_ordered_empty(Nodes nodes, DerOp1 in1, DerOp2 in2,
                                      Out out) {
    if constexpr (std::tuple_size_v<DerOp2> == 0) {
        return std::tuple_cat(out, in1);
    } else if constexpr (std::tuple_size_v<DerOp1> == 0) {
        return std::tuple_cat(out, in2);
    } else {
        return multiply_ordered_compare(nodes, in1, in2, out);
    }
}

template <class Nodes, class Id1, class... Ids1, std::size_t Order1,
          std::size_t... Orders1, std::size_t Power1, std::size_t... Powers1,
          class Id2, class... Ids2, std::size_t Order2, std::size_t... Orders2,
          std::size_t Power2, std::size_t... Powers2, class... Out>
constexpr auto
multiply_ordered_compare(Nodes nodes,
                         std::tuple<der2::p<Power1, der2::d<Order1, Id1>>,
                                    der2::p<Powers1, der2::d<Orders1, Ids1>>...>
                             in1,
                         std::tuple<der2::p<Power2, der2::d<Order2, Id2>>,
                                    der2::p<Powers2, der2::d<Orders2, Ids2>>...>
                             in2,
                         std::tuple<Out...> /* out */) {

    if constexpr (less_than(nodes, der2::p<Power1, der2::d<Order1, Id1>>{},
                            der2::p<Power2, der2::d<Order2, Id2>>{})) {

        return multiply_ordered_empty(
            nodes, in1,
            std::tuple<der2::p<Powers2, der2::d<Orders2, Ids2>>...>{},
            std::tuple<Out..., der2::p<Power2, der2::d<Order2, Id2>>>{});

    } else if constexpr (less_than(nodes,
                                   der2::p<Power2, der2::d<Order2, Id2>>{},
                                   der2::p<Power1, der2::d<Order1, Id1>>{})) {

        return multiply_ordered_empty(
            nodes, std::tuple<der2::p<Powers1, der2::d<Orders1, Ids1>>...>{},
            in2, std::tuple<Out..., der2::p<Power1, der2::d<Order1, Id1>>>{});

    } else {

        static_assert(
            std::is_same_v<der2::d<Order1, Id1>, der2::d<Order2, Id2>>);
        return multiply_ordered_empty(
            nodes, std::tuple<der2::p<Powers1, der2::d<Orders1, Ids1>>...>{},
            std::tuple<der2::p<Powers2, der2::d<Orders2, Ids2>>...>{},
            std::tuple<Out...,
                       der2::p<Power1 + Power2, der2::d<Order1, Id1>>>{});
    }
}

} // namespace detail

template <class Nodes> constexpr auto multiply_ordered(Nodes /* nodes */) {
    return std::tuple<>{};
};

template <class Nodes, class DerOp>
constexpr auto multiply_ordered(Nodes nodes, DerOp in1) {
    static_assert(is_ordered(nodes, in1));
    return in1;
};

template <class Nodes, class DerOp1, class DerOp2>
constexpr auto multiply_ordered(Nodes nodes, DerOp1 in1, DerOp2 in2) {
    static_assert(is_ordered(nodes, in1));
    static_assert(is_ordered(nodes, in2));
    return detail::multiply_ordered_empty(nodes, in1, in2, std::tuple<>{});
};

template <class Nodes, class DerOp1, class DerOp2, class... DerOps>
constexpr auto multiply_ordered(Nodes nodes, DerOp1 in1, DerOp2 in2,
                                DerOps... ins) {
    static_assert(is_ordered(nodes, in1));
    static_assert(is_ordered(nodes, in2));
    constexpr auto mult = multiply_ordered(nodes, in1, in2);
    return multiply_ordered(nodes, mult, ins...);
};

template <class Nodes, class... DerOps1, class ToMult>
constexpr auto multiply_ordered_tuple(Nodes nodes,
                                      std::tuple<DerOps1...> /* in */,
                                      ToMult in2) {
    return std::make_tuple(multiply_ordered(nodes, DerOps1{}, in2)...);
};

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_SORTED_MULT_HPP
