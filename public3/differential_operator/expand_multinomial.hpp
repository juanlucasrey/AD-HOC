#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_EXPAND_MULTINOMIAL_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_EXPAND_MULTINOMIAL_HPP

#include "../partition/combinations.hpp"
#include "../partition/multinomial_coefficient_index_sequence.hpp"
#include "is_ordered.hpp"
#include "ordered_mult.hpp"

#include <cstddef>
#include <tuple>
#include <utility>

namespace adhoc3 {

namespace detail {

template <std::size_t Idx, std::size_t End, class Nodes, class... Ids,
          std::size_t... I, class... Out>
constexpr auto expand_multinomial_aux(Nodes nodes, std::tuple<Ids...> in,
                                      std::index_sequence<I...> prev,
                                      std::tuple<Out...> out) {
    auto constexpr curent = multiply_ordered(nodes, pow<I>(Ids{})...);
    if constexpr (Idx == End) {
        return std::tuple_cat(out, std::make_tuple(curent));
    } else {
        auto constexpr next = NextMultinomialIS(prev);
        return expand_multinomial_aux<Idx + 1, End>(
            nodes, in, next, std::tuple_cat(out, std::make_tuple(curent)));
    }
};

} // namespace detail

template <std::size_t Power, class Nodes, class... Ids>
constexpr auto expand_multinomial(Nodes nodes, std::tuple<Ids...> in) {
    static_assert(is_ordered(nodes, in));
    constexpr auto first = FirstMultinomialIS<sizeof...(Ids), Power>();
    return detail::expand_multinomial_aux<1,
                                          combinations(sizeof...(Ids), Power)>(
        nodes, in, first, std::tuple<>{});
}

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_EXPAND_MULTINOMIAL_HPP
