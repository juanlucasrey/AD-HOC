#ifndef BASE3_HPP
#define BASE3_HPP

#include <cstddef>

namespace adhoc3 {

template <class Input1, class Input2> class add_t;
template <class Input1, class Input2> class sub_t;
template <class Input1, class Input2> class mul_t;
template <class Input1, class Input2> class div_t;

// it would be nice to add operator to constexpr float types in this class
// however this proposal would have to pass
// https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2019/p1045r1.html
// https://github.com/cplusplus/papers/issues/603
template <class Derived> struct Base {
    template <class Derived2>
    auto operator+(Base<Derived2> const &rhs) const
        -> add_t<const Derived, const Derived2>;

    template <class Derived2>
    auto operator-(Base<Derived2> const &rhs) const
        -> sub_t<const Derived, const Derived2>;

    template <class Derived2>
    auto operator*(Base<Derived2> const &rhs) const
        -> mul_t<const Derived, const Derived2>;

    template <class Derived2>
    auto operator/(Base<Derived2> const &rhs) const
        -> div_t<const Derived, const Derived2>;
};

template <class Derived>
template <class Derived2>
inline auto Base<Derived>::operator+(Base<Derived2> const & /* rhs */) const
    -> add_t<const Derived, const Derived2> {
    return add_t<const Derived, const Derived2>{};
}

template <class Derived>
template <class Derived2>
inline auto Base<Derived>::operator-(Base<Derived2> const & /* rhs */) const
    -> sub_t<const Derived, const Derived2> {
    return sub_t<const Derived, const Derived2>{};
}

template <class Derived>
template <class Derived2>
inline auto Base<Derived>::operator*(Base<Derived2> const & /* rhs */) const
    -> mul_t<const Derived, const Derived2> {
    return mul_t<const Derived, const Derived2>{};
}

template <class Derived>
template <class Derived2>
inline auto Base<Derived>::operator/(Base<Derived2> const & /* rhs */) const
    -> div_t<const Derived, const Derived2> {
    return div_t<const Derived, const Derived2>{};
}

} // namespace adhoc3

#endif // BASE3_HPP