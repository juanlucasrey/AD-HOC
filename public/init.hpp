#ifndef ADHOC_INIT_HPP
#define ADHOC_INIT_HPP

#include "adhoc.hpp"

#include <tuple>
#include <utility>

namespace adhoc {

namespace detail {

template <std::size_t... Idxs>
auto inline Init_aux(std::index_sequence<Idxs...> const & /* in */) {
    return std::tuple<double_t<Idxs>...>{};
}

} // namespace detail

template <std::size_t N> auto inline Init() {
    return detail::Init_aux(std::make_index_sequence<N>{});
}

} // namespace adhoc

#endif // ADHOC_INIT_HPP