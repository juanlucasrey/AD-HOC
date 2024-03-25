#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_IS_ORDERED_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_IS_ORDERED_HPP

#include "../adhoc.hpp"
#include "../dependency.hpp"
#include "differential_operator.hpp"
#include "less_than.hpp"

namespace adhoc3 {

namespace detail {

template <class Nodes>
constexpr auto is_ordered_aux(Nodes /* nodes */, std::tuple<> /* in */) {
    return true;
}

template <class Nodes, class Id>
constexpr auto is_ordered_aux(Nodes /* nodes */, std::tuple<Id> /* in */) {
    return true;
}

template <class Nodes, class Id1, class Id2, class... Ids>
constexpr auto is_ordered_aux(Nodes nodes,
                              std::tuple<Id1, Id2, Ids...> /* in */) {
    if constexpr (less_than(nodes, Id2{}, Id1{})) {
        return is_ordered_aux(nodes, std::tuple<Id2, Ids...>{});
    } else {
        return false;
    }
}

} // namespace detail

template <class Nodes, class... Ids>
constexpr auto is_ordered(Nodes nodes, std::tuple<Ids...> in) {
    return detail::is_ordered_aux(nodes, in);
}

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_IS_ORDERED_HPP
