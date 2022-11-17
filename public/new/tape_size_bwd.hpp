#ifndef ADHOC_TAPE_SIZE_BWD_HPP
#define ADHOC_TAPE_SIZE_BWD_HPP

// #include "adhoc.hpp"
#include "calc_order.hpp"
#include "dependency.hpp"
// #include "utils.hpp"

#include <algorithm>
#include <tuple>

namespace adhoc2 {

namespace detail {

template <class... ActiveLeafsAndRoots, class... Nodes>
constexpr static auto active_nodes(std::tuple<Nodes...> /* in */) {
    return std::tuple_cat(
        std::conditional_t<(depends_many<Nodes, ActiveLeafsAndRoots...>()),
                           std::tuple<Nodes>, std::tuple<>>{}...);
}

} // namespace detail

template <class... ActiveLeafsAndRoots>
constexpr static auto tape_size_bwd() /* -> std::size_t */ {

    auto constexpr active_leaves_tuple = std::tuple_cat(
        std::conditional_t<
            (depends_many<ActiveLeafsAndRoots, ActiveLeafsAndRoots...>()),
            std::tuple<>, std::tuple<ActiveLeafsAndRoots>>{}...);

    auto constexpr nodes =
        detail::calc_order_aux_t<ActiveLeafsAndRoots...>::template call();

    auto constexpr active_nodes =
        detail::active_nodes<ActiveLeafsAndRoots...>(nodes);
    return active_nodes;
}

} // namespace adhoc2

#endif // ADHOC_TAPE_SIZE_BWD_HPP