#ifndef ADHOC3_DERS_HPP
#define ADHOC3_DERS_HPP

#include "adhoc.hpp"

namespace adhoc3 {

template <class... Types> struct mvar_t {};
template <class Id, std::size_t Order> struct IdAndOrder_t {};

namespace detail {

template <class... Args> auto duplicate_id(Args...) -> bool;

template <> auto constexpr duplicate_id() -> bool { return false; }

template <std::size_t Order1, class IdDerived, std::size_t... Orders,
          class... IdsDerived>
auto constexpr duplicate_id(mvar_t<IdAndOrder_t<Base<IdDerived>, Order1>>,
                            mvar_t<IdAndOrder_t<Base<IdsDerived>, Orders>>...)
    -> bool {
    bool thisisduplicate = (std::is_same<IdDerived, IdsDerived>::value || ...);
    bool duplicateinnext =
        duplicate_id(mvar_t<IdAndOrder_t<Base<IdsDerived>, Orders>>{}...);
    return thisisduplicate || duplicateinnext;
}

} // namespace detail

template <std::size_t Order = 1, class Derived> auto var(Base<Derived>) {
    return mvar_t<IdAndOrder_t<Base<Derived>, Order>>{};
}

template <std::size_t... Orders, class... IdsDerived>
auto var(mvar_t<IdAndOrder_t<Base<IdsDerived>, Orders>>...) {

    static_assert(!detail::duplicate_id(
                      mvar_t<IdAndOrder_t<Base<IdsDerived>, Orders>>{}...),
                  "duplicate ids on multivariate derivative declaration");

    return mvar_t<IdAndOrder_t<Base<IdsDerived>, Orders>...>{};
}

} // namespace adhoc3

#endif // ADHOC3_DERS_HPP