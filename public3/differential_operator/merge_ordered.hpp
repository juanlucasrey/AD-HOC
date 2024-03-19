#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_MERGE_ORDERED_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_MERGE_ORDERED_HPP

#include "../adhoc.hpp"
#include "../dependency.hpp"
#include "differential_operator.hpp"
#include "is_ordered.hpp"
#include "less_than.hpp"

namespace adhoc3 {

namespace detail {

template <class Id1, class... Ids1, class Id2, class... Ids2, class... Out,
          class... CalculationNodes>
constexpr auto merge_ordered_compare(std::tuple<Id1, Ids1...> in1,
                                     std::tuple<Id2, Ids2...> in2,
                                     std::tuple<Out...> out,
                                     std::tuple<CalculationNodes...> nodes);

template <class... Ids1, class... Ids2, class... Out, class... CalculationNodes>
constexpr auto merge_ordered_empty(std::tuple<Ids1...> in1,
                                   std::tuple<Ids2...> in2,
                                   std::tuple<Out...> out,
                                   std::tuple<CalculationNodes...> nodes) {

    if constexpr (sizeof...(Ids2) == 0) {
        return std::tuple_cat(out, in1);
    } else if constexpr (sizeof...(Ids1) == 0) {
        return std::tuple_cat(out, in2);
    } else {
        return merge_ordered_compare(in1, in2, out, nodes);
    }
}

template <class Id1, class... Ids1, class Id2, class... Ids2, class... Out,
          class... CalculationNodes>
constexpr auto merge_ordered_compare(std::tuple<Id1, Ids1...> in1,
                                     std::tuple<Id2, Ids2...> in2,
                                     std::tuple<Out...> /* out */,
                                     std::tuple<CalculationNodes...> nodes) {
    if constexpr (std::is_same<Id1, Id2>::value) {
        return merge_ordered_empty(std::tuple<Ids1...>{}, std::tuple<Ids2...>{},
                                   std::tuple<Out..., Id1>{}, nodes);

    } else if constexpr (less_than(Id1{}, Id2{}, nodes)) {
        return merge_ordered_empty(in1, std::tuple<Ids2...>{},
                                   std::tuple<Out..., Id2>{}, nodes);

    } else {
        return merge_ordered_empty(std::tuple<Ids1...>{}, in2,
                                   std::tuple<Out..., Id1>{}, nodes);
    }
}

} // namespace detail

template <class... Ids1, class... Ids2, class... CalculationNodes>
constexpr auto merge_ordered(std::tuple<Ids1...> in1, std::tuple<Ids2...> in2,
                             std::tuple<CalculationNodes...> nodes) {

    static_assert(is_ordered(in1, nodes));
    static_assert(is_ordered(in2, nodes));
    return detail::merge_ordered_empty(in1, in2, std::tuple<>{}, nodes);
}

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_MERGE_ORDERED_HPP
