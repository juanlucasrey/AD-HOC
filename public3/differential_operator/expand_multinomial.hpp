#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_EXPAND_MULTINOMIAL_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_EXPAND_MULTINOMIAL_HPP

#include "../partition/combinations.hpp"
#include "../partition/multinomial_coefficient_index_sequence.hpp"
#include "ordered_mult.hpp"

#include <cstddef>
#include <tuple>
#include <utility>

namespace adhoc3 {

namespace detail {

template <class Nodes, class... Ids, std::size_t... I>
constexpr auto expand_single(Nodes nodes, std::tuple<Ids...> /* in */,
                             std::index_sequence<I...> /* i */) {
    static_assert(sizeof...(Ids) == sizeof...(I));
    return multiply_ordered(nodes, pow<I>(Ids{})...);
}

template <std::size_t Idx, std::size_t End, class Nodes, class... Ids,
          class Prev, class... Out>
constexpr auto expand_multinomial(Nodes nodes, std::tuple<Ids...> in, Prev prev,
                                  std::tuple<Out...> out) {
    auto constexpr curent = expand_single(nodes, in, prev);
    if constexpr (Idx == End) {
        return std::tuple_cat(out, std::make_tuple(curent));
        // return std::tuple<Out..., Prev>{};
    } else {
        auto constexpr next = NextMultinomialIS(prev);
        // return expand_multinomial<Idx + 1, End>(nodes, in, next,
        //                                         std::tuple<Out..., Prev>{});
        return expand_multinomial<Idx + 1, End>(
            nodes, in, next, std::tuple_cat(out, std::make_tuple(curent)));
    }
};

} // namespace detail

template <std::size_t Power, class Nodes, class... Ids>
constexpr auto expand_multinomial(Nodes nodes, std::tuple<Ids...> in) {
    constexpr auto first = FirstMultinomialIS<sizeof...(Ids), Power>();
    return detail::expand_multinomial<1, combinations(sizeof...(Ids), Power)>(
        nodes, in, first, std::tuple<>{});
}

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_EXPAND_MULTINOMIAL_HPP
