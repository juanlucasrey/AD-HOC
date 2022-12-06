#ifndef ADHOC_TAPE_HPP
#define ADHOC_TAPE_HPP

#include "calc_order.hpp"
#include "leafs.hpp"
#include "tape_size_bwd.hpp"
#include "tape_size_fwd.hpp"
#include "utils.hpp"

#include <array>
#include <iostream>
#include <tuple>

namespace adhoc {

namespace detail {

template <class ValuesTuple, class ActiveLeafsAndRootsDerivativesTuple,
          std::size_t MaxWidth>
class Tape2;

template <class... Values, class... ActiveLeafsAndRootsDerivatives,
          std::size_t MaxWidth>
class Tape2<std::tuple<Values...>,
            std::tuple<ActiveLeafsAndRootsDerivatives...>, MaxWidth> {
    std::array<double, sizeof...(Values)> m_values{};
    std::array<double, sizeof...(ActiveLeafsAndRootsDerivatives)>
        m_derivatives{};

  public:
    std::array<double, MaxWidth> m_buffer{};

    template <class Derived> auto inline val(Base<Derived> var) const -> double;
    template <class Derived> auto inline val(Base<Derived> var) -> double &;
    template <class Derived> auto inline der(Base<Derived> var) const -> double;
    template <class Derived> auto inline der(Base<Derived> var) -> double &;

    template <class Derived>
    auto inline val_check(Base<Derived> var) -> double &;
    template <class Derived>
    auto inline der_check(Base<Derived> var) -> double &;
};

template <class... Values, class... ActiveLeafsAndRootsDerivatives,
          std::size_t MaxWidth>
template <class Derived>
auto Tape2<std::tuple<Values...>, std::tuple<ActiveLeafsAndRootsDerivatives...>,
           MaxWidth>::val(Base<Derived> /* in */) const -> double {
    if constexpr (has_type<Derived, Values...>()) {
        constexpr auto idx = idx_type<Derived, Values...>();
        return this->m_values[idx];
    } else {
        return 0;
    }
}

template <class... Values, class... ActiveLeafsAndRootsDerivatives,
          std::size_t MaxWidth>
template <class Derived>
auto Tape2<std::tuple<Values...>, std::tuple<ActiveLeafsAndRootsDerivatives...>,
           MaxWidth>::val(Base<Derived> /* in */) -> double & {
    if constexpr (has_type<Derived, Values...>()) {
        constexpr auto idx = idx_type<Derived, Values...>();
        return this->m_values[idx];
    } else {
        return this->m_buffer[0];
    }
}

template <class... Values, class... ActiveLeafsAndRootsDerivatives,
          std::size_t MaxWidth>
template <class Derived>
auto Tape2<std::tuple<Values...>, std::tuple<ActiveLeafsAndRootsDerivatives...>,
           MaxWidth>::der(Base<Derived> /* in */) const -> double {
    if constexpr (has_type<Derived, ActiveLeafsAndRootsDerivatives...>()) {
        constexpr auto idx =
            idx_type<Derived, ActiveLeafsAndRootsDerivatives...>();
        return this->m_derivatives[idx];
    } else {
        return 0;
    }
}

template <class... Values, class... ActiveLeafsAndRootsDerivatives,
          std::size_t MaxWidth>
template <class Derived>
auto Tape2<std::tuple<Values...>, std::tuple<ActiveLeafsAndRootsDerivatives...>,
           MaxWidth>::der(Base<Derived> /* in */) -> double & {
    if constexpr (has_type<Derived, ActiveLeafsAndRootsDerivatives...>()) {
        constexpr auto idx =
            idx_type<Derived, ActiveLeafsAndRootsDerivatives...>();
        return this->m_derivatives[idx];
    } else {
        return this->m_buffer[0];
    }
}

template <class... Values, class... ActiveLeafsAndRootsDerivatives,
          std::size_t MaxWidth>
template <class Derived>
auto Tape2<std::tuple<Values...>, std::tuple<ActiveLeafsAndRootsDerivatives...>,
           MaxWidth>::val_check(Base<Derived> /* in */) -> double & {
    static_assert(has_type<Derived, Values...>());
    constexpr auto idx = idx_type<Derived, Values...>();
    return this->m_values[idx];
}

template <class... Values, class... ActiveLeafsAndRootsDerivatives,
          std::size_t MaxWidth>
template <class Derived>
auto Tape2<std::tuple<Values...>, std::tuple<ActiveLeafsAndRootsDerivatives...>,
           MaxWidth>::der_check(Base<Derived> /* in */) -> double & {
    static_assert(has_type<Derived, ActiveLeafsAndRootsDerivatives...>(),
                  "variable not active");
    constexpr auto idx = idx_type<Derived, ActiveLeafsAndRootsDerivatives...>();
    return this->m_derivatives[idx];
}

template <constants::ArgType D, class Tape>
inline auto val(Tape const & /* in */, constants::CD<D> /* node */) -> double {
    return constants::CD<D>::v();
}

template <class Node, class Tape>
inline auto val(Tape const &in, Node /* node */) -> double {
    return in.val(Node{});
}

} // namespace detail

template <class... ActiveLeafsAndRoots>
auto inline Tape(ActiveLeafsAndRoots... /* out */) {
    constexpr auto leafs =
        detail::leafs_t<ActiveLeafsAndRoots...>::template call();
    constexpr auto nodes =
        detail::calc_order_aux_t<true, ActiveLeafsAndRoots...>::template call();
    constexpr auto intermediate_deriv_tape_size =
        tape_size_bwd<ActiveLeafsAndRoots...>();
    return detail::Tape2<decltype(std::tuple_cat(leafs, nodes)),
                         std::tuple<ActiveLeafsAndRoots...>,
                         intermediate_deriv_tape_size>();
}

} // namespace adhoc

#endif // ADHOC_TAPE_HPP