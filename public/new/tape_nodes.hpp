#ifndef ADHOC_TAPE_NODES_HPP
#define ADHOC_TAPE_NODES_HPP

#include "base.hpp"
#include "leafs.hpp"
#include "utils.hpp"

#include <array>
#include <type_traits>
#include <utility>

namespace adhoc2 {

template <class... Nodes> class Tape {
    std::array<double, sizeof...(Nodes)> buff{};

  public:
    explicit Tape(Nodes... /* in */) {}
    template <class Derived> void inline set(Base<Derived> var, double val);
    template <class Derived> auto inline get(Base<Derived> var) const -> double;
    template <class Derived> auto inline get(Base<Derived> var) -> double &;
};

template <class... Nodes>
template <class Derived>
void Tape<Nodes...>::set(Base<Derived> /* in */, double val) {
    constexpr auto idx = idx_type<Derived, Nodes...>();
    this->buff[idx] = val;
}

template <class... Nodes>
template <class Derived>
auto Tape<Nodes...>::get(Base<Derived> /* in */) const -> double {
    constexpr auto idx = idx_type<Derived, Nodes...>();
    return this->buff[idx];
}

template <class... Nodes>
template <class Derived>
auto Tape<Nodes...>::get(Base<Derived> /* in */) -> double & {
    constexpr auto idx = idx_type<Derived, Nodes...>();
    return this->buff[idx];
}
namespace detail {
template <class... Roots, class... Leafs>
auto constexpr TapeRootsAndLeafs_aux(std::tuple<Leafs...> /* in */) {
    return Tape(Leafs{}..., Roots{}...);
}

} // namespace detail

template <class... Roots> auto inline TapeRootsAndLeafs(Roots &.../* out */) {
    auto constexpr leafs = detail::leafs_t<Roots...>::template call();
    return detail::TapeRootsAndLeafs_aux<Roots...>(leafs);
}

template <class... ActiveLeafsAndRoots>
auto inline TapeDerivatives(ActiveLeafsAndRoots &.../* in */) {
    return Tape(ActiveLeafsAndRoots{}...);
}

} // namespace adhoc2

#endif // ADHOC_TAPE_NODES_HPP