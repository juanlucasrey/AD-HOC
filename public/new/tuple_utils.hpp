#ifndef ADHOC_TUPLE_UTILS_HPP
#define ADHOC_TUPLE_UTILS_HPP

#include "utils.hpp"

#include <iostream>
#include <tuple>

namespace adhoc2 {

namespace detail {

template <std::size_t IdxToReplace, class Replacement, class... TupleValues,
          std::size_t... IndicesStart, std::size_t... IndicesEnd>
auto constexpr replace_first_aux2(
    std::tuple<TupleValues...> /* in */,
    std::index_sequence<IndicesStart...> /* idx1 */,
    std::index_sequence<IndicesEnd...> /* idx2 */) {
    return std::make_tuple(
        std::tuple_element_t<IndicesStart, std::tuple<TupleValues...>>{}...,
        Replacement{},
        std::tuple_element_t<IndicesEnd + IdxToReplace + 1,
                             std::tuple<TupleValues...>>{}...);
}

template <std::size_t IdxToReplace, class Replacement, class... TupleValues,
          class IndicesStart = std::make_index_sequence<IdxToReplace>,
          class IndicesEnd = std::make_index_sequence<sizeof...(TupleValues) -
                                                      IdxToReplace - 1>>
auto constexpr replace_first_aux(std::tuple<TupleValues...> /* in */) {
    return replace_first_aux2<IdxToReplace, Replacement>(
        std::tuple<TupleValues...>{}, IndicesStart{}, IndicesEnd{});
}

} // namespace detail

template <class TypeToReplace, class Replacement, class... TupleValues>
auto constexpr replace_first(std::tuple<TupleValues...> /* in */) {
    static_assert(has_type2<TypeToReplace, TupleValues...>(),
                  "type to replace not on this tuple");
    constexpr auto idx =
        detail::get_index<TypeToReplace, TupleValues...>::value;
    return detail::replace_first_aux<idx, Replacement>(
        std::tuple<TupleValues...>{});
}

} // namespace adhoc2
#endif // ADHOC_TUPLE_UTILS_HPP