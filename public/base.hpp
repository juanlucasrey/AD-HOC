#ifndef BASE_HPP
#define BASE_HPP

#include <cstddef>

namespace adhoc2 {

template <class Input1, class Input2> class add;
template <class Input1, class Input2> class mul;
template <class Input1, class Input2> class div;

template <class Derived> struct Base {
    using is_adhoc_tag = void;

    template <class Derived2>
    auto operator+(Base<Derived2> const &rhs) const
        -> add<const Derived, const Derived2>;

    template <class Derived2>
    auto operator*(Base<Derived2> const &rhs) const
        -> mul<const Derived, const Derived2>;

    template <class Derived2>
    auto operator/(Base<Derived2> const &rhs) const
        -> div<const Derived, const Derived2>;

    template <class... Denom>
    constexpr static auto depends3() noexcept -> std::size_t;
};

template <class Derived>
template <class Derived2>
inline auto Base<Derived>::operator+(Base<Derived2> const &rhs) const
    -> add<const Derived, const Derived2> {
    return {*static_cast<Derived const *>(this),
            *static_cast<Derived2 const *>(&rhs)};
}

template <class Derived>
template <class Derived2>
inline auto Base<Derived>::operator*(Base<Derived2> const &rhs) const
    -> mul<const Derived, const Derived2> {
    return {*static_cast<Derived const *>(this),
            *static_cast<Derived2 const *>(&rhs)};
}

template <class Derived>
template <class Derived2>
inline auto Base<Derived>::operator/(Base<Derived2> const &rhs) const
    -> div<const Derived, const Derived2> {
    return {*static_cast<Derived const *>(this),
            *static_cast<Derived2 const *>(&rhs)};
}

template <class Derived>
template <class... Denom>
inline constexpr auto Base<Derived>::depends3() noexcept -> std::size_t {
    if constexpr (sizeof...(Denom) == 0) {
        return 0UL;
    } else {
        return (Derived::template depends<Denom>() + ...);
    }
}

} // namespace adhoc2

#endif // BASE_HPP