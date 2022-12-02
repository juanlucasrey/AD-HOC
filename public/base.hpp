#ifndef ADHOC_BASE_HPP
#define ADHOC_BASE_HPP

#include <cstddef>

namespace adhoc {

template <class Input1, class Input2> struct add_t;
template <class Input1, class Input2> struct sub_t;
template <class Input1, class Input2> struct mul_t;
template <class Input1, class Input2> struct div_t;

// it would be nice to add operator to constexpr float types in this class
// however this proposal would have to pass
// https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2019/p1045r1.html
// https://github.com/cplusplus/papers/issues/603
template <class Derived> struct Base {
    template <class Derived2> auto operator+(Base<Derived2> /* rhs */) const {
        return add_t<Derived, Derived2>{};
    }

    template <class Derived2> auto operator-(Base<Derived2> /* rhs */) const {
        return sub_t<Derived, Derived2>{};
    }

    template <class Derived2> auto operator*(Base<Derived2> /* rhs */) const {
        return mul_t<Derived, Derived2>{};
    }

    template <class Derived2> auto operator/(Base<Derived2> /* rhs */) const {
        return div_t<Derived, Derived2>{};
    }
};

} // namespace adhoc

#endif // ADHOC_BASE_HPP