#ifndef BASE_HPP
#define BASE_HPP

#include <cstddef>

namespace adhoc2 {

template <class Input1, class Input2> class add;
template <class Input1, class Input2> class subs;
template <class Input1, class Input2> class mul;
template <class Input1, class Input2> class div;

// it would be nice to add operator to constexpr float types in this class
// however this proposal would have to pass
// https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2019/p1045r1.html
// https://github.com/cplusplus/papers/issues/603
template <class Derived> struct Base {
    using is_adhoc_tag = void;

    template <class Derived2>
    auto operator+(Base<Derived2> const &rhs) const
        -> add<const Derived, const Derived2>;

    template <class Derived2>
    auto operator-(Base<Derived2> const &rhs) const
        -> subs<const Derived, const Derived2>;

    template <class Derived2>
    auto operator*(Base<Derived2> const &rhs) const
        -> mul<const Derived, const Derived2>;

    template <class Derived2>
    auto operator/(Base<Derived2> const &rhs) const
        -> div<const Derived, const Derived2>;
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
inline auto Base<Derived>::operator-(Base<Derived2> const &rhs) const
    -> subs<const Derived, const Derived2> {
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

} // namespace adhoc2

#endif // BASE_HPP