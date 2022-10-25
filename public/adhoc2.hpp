#ifndef ADHOC2_HPP
#define ADHOC2_HPP

#include <filter.hpp>
#include <utils.hpp>

#include <array>
#include <iostream>
#include <type_traits>

namespace adhoc2 {

template <class Input, class Derived> class Univariate {
  public:
    template <class Denom> constexpr static auto depends() noexcept -> bool;
    template <class Denom> constexpr static auto depends_in() noexcept -> bool;
};

template <class Input, class Derived>
template <class Denom>
inline constexpr auto Univariate<Input, Derived>::depends() noexcept -> bool {
    return std::is_same_v<Denom, Derived> ||
           std::is_same_v<Denom, Derived const> ||
           Input::template depends<Denom>();
}

template <class Input, class Derived>
template <class Denom>
inline constexpr auto Univariate<Input, Derived>::depends_in() noexcept
    -> bool {
    return Input::template depends<Denom>();
}

template <class Input1, class Input2, class Derived> class Bivariate {
  public:
    template <class Denom> constexpr static auto depends() noexcept -> bool;
    template <class Denom> constexpr static auto depends_in() noexcept -> bool;
};

template <class Input1, class Input2, class Derived>
template <class Denom>
inline constexpr auto Bivariate<Input1, Input2, Derived>::depends() noexcept
    -> bool {
    return std::is_same_v<Denom, Derived> ||
           std::is_same_v<Denom, Derived const> ||
           Input1::template depends<Denom>() ||
           Input2::template depends<Denom>();
}

template <class Input1, class Input2, class Derived>
template <class Denom>
inline constexpr auto Bivariate<Input1, Input2, Derived>::depends_in() noexcept
    -> bool {
    return Input1::template depends<Denom>() ||
           Input2::template depends<Denom>();
}

template <class Input1, class Input2> class mul;
template <class Input1, class Input2> class add;

template <class Derived> class Base {
  protected:
    double m_value{0};

  public:
    explicit Base(double value);
    [[nodiscard]] auto v() const noexcept -> double;

    template <class Derived2>
    auto operator+(Base<Derived2> const &rhs) const
        -> add<const Derived, const Derived2>;

    template <class Derived2>
    auto operator*(Base<Derived2> const &rhs) const
        -> mul<const Derived, const Derived2>;

    template <class... Denom>
    constexpr static auto depends3() noexcept -> std::size_t;
};

template <class Derived> Base<Derived>::Base(double value) : m_value(value) {}

template <class Derived>
inline auto Base<Derived>::v() const noexcept -> double {
    return this->m_value;
}

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
template <class... Denom>
inline constexpr auto Base<Derived>::depends3() noexcept -> std::size_t {
    if constexpr (sizeof...(Denom) == 0) {
        return 0UL;
    } else {
        return (Derived::template depends<Denom>() + ...);
    }
}

template <class Input>
class exp_t : public Base<exp_t<Input>>,
              public Univariate<Input, exp_t<Input>> {
  private:
    Input const &m_active;

  public:
    explicit exp_t(Input const &active);
    auto d() -> double;
};

template <class Input>
exp_t<Input>::exp_t(Input const &active)
    : Base<exp_t<Input>>(std::exp(active.v())), m_active(active) {}

template <class Input> auto exp_t<Input>::d() -> double { return this->v(); }

template <class Derived>
auto exp(Base<Derived> const &in) -> exp_t<const Derived> {
    return exp_t<const Derived>(*static_cast<const Derived *>(&in));
}

template <class Input>
class cos_t : public Base<cos_t<Input>>,
              public Univariate<Input, cos_t<Input>> {
  private:
    Input const &m_active;

  public:
    explicit cos_t(Input const &active);
    auto d() -> double;
};

template <class Input>
cos_t<Input>::cos_t(Input const &active)
    : Base<cos_t<Input>>(std::cos(active.v())), m_active(active) {}

template <class Input> auto cos_t<Input>::d() -> double {
    return -std::sin(this->m_active.v());
}

template <class Derived>
auto cos(Base<Derived> const &in) -> cos_t<const Derived> {
    return cos_t<const Derived>(*static_cast<const Derived *>(&in));
}

template <class Input1, class Input2>
class add : public Base<add<Input1, Input2>>,
            public Bivariate<Input1, Input2, add<Input1, Input2>> {
    Input1 const &m_active1;
    Input2 const &m_active2;

  public:
    add(Input1 const &active1, Input2 const &active2);
    auto d1() -> double;
    auto d2() -> double;
};

template <class Input1, class Input2>
add<Input1, Input2>::add(Input1 const &active1, Input2 const &active2)
    : Base<add<Input1, Input2>>(active1.v() + active2.v()), m_active1(active1),
      m_active2(active2) {}

template <class Input1, class Input2> auto add<Input1, Input2>::d1() -> double {
    return 1.0;
}

template <class Input1, class Input2> auto add<Input1, Input2>::d2() -> double {
    return 1.0;
}

template <class Input1, class Input2>
class mul : public Base<mul<Input1, Input2>>,
            public Bivariate<Input1, Input2, mul<Input1, Input2>> {
    Input1 const &m_active1;
    Input2 const &m_active2;

  public:
    mul(Input1 const &active1, Input2 const &active2);
    auto d1() -> double;
    auto d2() -> double;
};

template <class Input1, class Input2>
mul<Input1, Input2>::mul(Input1 const &active1, Input2 const &active2)
    : Base<mul<Input1, Input2>>(active1.v() * active2.v()), m_active1(active1),
      m_active2(active2) {}

template <class Input1, class Input2> auto mul<Input1, Input2>::d1() -> double {
    return this->m_active2.v();
}

template <class Input1, class Input2> auto mul<Input1, Input2>::d2() -> double {
    return this->m_active1.v();
}

namespace detail {
template <int N>
class adouble_aux : public Base<adouble_aux<N>>,
                    public Univariate<adouble_aux<N>, adouble_aux<N>> {
  public:
    explicit adouble_aux(double value);

    template <class Denom> constexpr static auto depends() noexcept -> bool;
};

template <int N>
inline adouble_aux<N>::adouble_aux(double value)
    : Base<adouble_aux<N>>(value) {}

template <int N>
template <class Denom>
inline constexpr auto adouble_aux<N>::depends() noexcept -> bool {
    return std::is_same_v<Denom, adouble_aux<N>>;
}

} // namespace detail

template <typename... Types> struct args {};

#define adouble detail::adouble_aux<__COUNTER__>

} // namespace adhoc2

#endif // ADHOC2_HPP