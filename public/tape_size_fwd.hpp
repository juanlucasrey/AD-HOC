#ifndef ADHOC_TAPE_SIZE_FWD_HPP
#define ADHOC_TAPE_SIZE_FWD_HPP

#include "calc_order.hpp"
#include "dependency.hpp"
#include "tape_nodes.hpp"
#include "utils.hpp"

#include <tuple>

namespace adhoc {

namespace detail {

template <class... Roots> struct tape_size_fwd_aux_t {
    template <class... Leafs, class... Nodes>
    constexpr static auto call(std::tuple<Leafs...> /* in1 */,
                               std::tuple<Nodes...> /* in2 */) {
        return std::tuple_cat(
            std::conditional_t<has_type<Nodes, Leafs..., Roots...>(),
                               std::tuple<>, std::tuple<Nodes>>{}...);
    }
};

} // namespace detail

template <class... Roots> constexpr auto tape_size_fwd_t(Roots... /* in */) {
    auto constexpr leafs = detail::leafs_t<Roots...>::template call();
    auto constexpr nodes =
        detail::calc_order_aux_t<true, Roots...>::template call();

    return detail::tape_size_fwd_aux_t<Roots...>::template call(leafs, nodes);
}

} // namespace adhoc

#endif // ADHOC_TAPE_SIZE_FWD_HPP