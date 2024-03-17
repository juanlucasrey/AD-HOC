#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_IS_ORDERED_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_IS_ORDERED_HPP

#include "../adhoc.hpp"
#include "../dependency.hpp"
#include "differential_operator.hpp"
#include "less_than.hpp"

namespace adhoc3 {

namespace detail {

template <class... CalculationNodes>
constexpr auto is_ordered_aux(std::tuple<> /* in */,
                              std::tuple<CalculationNodes...> /* nodes */) {

    return true;
}

template <class Id, class... CalculationNodes>
constexpr auto is_ordered_aux(std::tuple<Id> /* in */,
                              std::tuple<CalculationNodes...> /* nodes */) {

    return true;
}

template <class Id1, class Id2, class... Ids, class... CalculationNodes>
constexpr auto is_ordered_aux(std::tuple<Id1, Id2, Ids...> /* in */,
                              std::tuple<CalculationNodes...> nodes) {

    if constexpr (less_than(Id2{}, Id1{}, nodes)) {
        return is_ordered_aux(std::tuple<Id2, Ids...>{}, nodes);
    } else {
        return false;
    }
}

} // namespace detail

template <class... Ids, class... CalculationNodes>
constexpr auto is_ordered(std::tuple<Ids...> in,
                          std::tuple<CalculationNodes...> nodes) {

    return detail::is_ordered_aux(in, nodes);
}

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_IS_ORDERED_HPP
