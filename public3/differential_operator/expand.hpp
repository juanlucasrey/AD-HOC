#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_EXPAND_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_EXPAND_HPP

#include "../adhoc.hpp"
#include "../dependency.hpp"
#include "differential_operator.hpp"
#include "select_root_derivatives.hpp"

namespace adhoc3 {

// namespace detail {

// template <class... Ids, std::size_t... Orders, std::size_t... Powers>
// // constexpr auto expand_aux(
// auto expand_aux(
//     std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...> /* ids1 */) {
//     return std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...>{};
// }

// template <class Id1, class Id2, std::size_t Order1, std::size_t Power1,
//           class... Ids, std::size_t... Orders, std::size_t... Powers>
// // constexpr auto expand_aux(
// auto expand_aux(
//     std::tuple<der2::p<Power1, der2::d<Order1, add_t<Id1, Id2>>>,
//                der2::p<Powers, der2::d<Orders, Ids>>...> /* ids1 */) {
//     return std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...>{};
// }

// } // namespace detail

template <class... Ids, std::size_t... Orders, std::size_t... Powers>
// constexpr auto expand(
auto expand(std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...> /* ids1 */) {
    if constexpr (is_root_class_v<
                      std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...>>) {
        return std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...>{};
    } else {
        auto intermediate =
            expand_aux(std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...>{});
        return expand(intermediate);
    }
}

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_EXPAND_HPP
