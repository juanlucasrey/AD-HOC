#ifndef ADHOC2_HPP
#define ADHOC2_HPP

#include <filter.hpp>
#include <functions/constants.hpp>
#include <utils.hpp>

#include <array>
#include <cmath>
#include <iostream>
#include <type_traits>

#ifndef NDEBUG
#include <iostream>
#endif

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

template <class Input1, class Input2> class add;
template <class Input1, class Input2> class mul;
template <class Input1, class Input2> class div;

template <class Derived> class Base {
  protected:
    double m_value{0};

  public:
    using is_adhoc_tag = void;
    explicit Base(double value);
    [[nodiscard]] auto v() const noexcept -> double;

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

template <class Input>
class exp_t : public Base<exp_t<Input>>,
              public Univariate<Input, exp_t<Input>> {
  private:
    Input const m_active;

  public:
    explicit exp_t(Input const &active);
    [[nodiscard]] auto d() const -> double;
    [[nodiscard]] auto input() const -> Input const &;
};

template <class Input>
exp_t<Input>::exp_t(Input const &active)
    : Base<exp_t<Input>>(std::exp(active.v())), m_active(active) {}

template <class Input> auto exp_t<Input>::d() const -> double {
#ifndef NDEBUG
    std::cout << this->v();
#endif
    return this->v();
}

template <class Input> auto exp_t<Input>::input() const -> Input const & {
    return this->m_active;
}

template <class Derived>
auto exp(Base<Derived> const &in) -> exp_t<const Derived> {
    return exp_t<const Derived>(*static_cast<const Derived *>(&in));
}

template <class Input>
class cos_t : public Base<cos_t<Input>>,
              public Univariate<Input, cos_t<Input>> {
  private:
    Input const m_active;

  public:
    explicit cos_t(Input const &active);
    [[nodiscard]] auto d() const -> double;
    [[nodiscard]] auto input() const -> Input const &;
};

template <class Input>
cos_t<Input>::cos_t(Input const &active)
    : Base<cos_t<Input>>(std::cos(active.v())), m_active(active) {}

template <class Input> auto cos_t<Input>::d() const -> double {
#ifndef NDEBUG
    std::cout << "(-std::sin(" << this->m_active.v() << "))";
#endif
    return -std::sin(this->m_active.v());
}

template <class Input> auto cos_t<Input>::input() const -> Input const & {
    return this->m_active;
}

template <class Derived>
auto cos(Base<Derived> const &in) -> cos_t<const Derived> {
    return cos_t<const Derived>(*static_cast<const Derived *>(&in));
}

template <class Input>
class sqrt_t : public Base<sqrt_t<Input>>,
               public Univariate<Input, sqrt_t<Input>> {
  private:
    Input const m_active;

  public:
    explicit sqrt_t(Input const &active);
    [[nodiscard]] auto d() const -> double;
    [[nodiscard]] auto input() const -> Input const &;
};

template <class Input>
sqrt_t<Input>::sqrt_t(Input const &active)
    : Base<sqrt_t<Input>>(std::sqrt(active.v())), m_active(active) {}

template <class Input> auto sqrt_t<Input>::d() const -> double {
#ifndef NDEBUG
    std::cout << "(0.5 * " << this->v() << " / " << this->m_active.v() << " )";
#endif
    return 0.5 * this->v() / this->m_active.v();
}

template <class Input> auto sqrt_t<Input>::input() const -> Input const & {
    return this->m_active;
}

template <class Derived>
auto sqrt(Base<Derived> const &in) -> sqrt_t<const Derived> {
    return sqrt_t<const Derived>(*static_cast<const Derived *>(&in));
}

template <class Input>
class log_t : public Base<log_t<Input>>,
              public Univariate<Input, log_t<Input>> {
  private:
    Input const m_active;

  public:
    explicit log_t(Input const &active);
    [[nodiscard]] auto d() const -> double;
    [[nodiscard]] auto input() const -> Input const &;
};

template <class Input>
log_t<Input>::log_t(Input const &active)
    : Base<log_t<Input>>(std::log(active.v())), m_active(active) {}

template <class Input> auto log_t<Input>::d() const -> double {
#ifndef NDEBUG
    std::cout << "(1.0 / " << this->m_active.v() << ")";
#endif
    return 1.0 / this->m_active.v();
}

template <class Input> auto log_t<Input>::input() const -> Input const & {
    return this->m_active;
}

template <class Derived>
auto log(Base<Derived> const &in) -> log_t<const Derived> {
    return log_t<const Derived>(*static_cast<const Derived *>(&in));
}

template <class Input>
class erfc_t : public Base<erfc_t<Input>>,
               public Univariate<Input, erfc_t<Input>> {
  private:
    Input const m_active;

  public:
    explicit erfc_t(Input const &active);
    [[nodiscard]] auto d() const -> double;
    [[nodiscard]] auto input() const -> Input const &;
};

template <class Input>
erfc_t<Input>::erfc_t(Input const &active)
    : Base<erfc_t<Input>>(std::erfc(active.v())), m_active(active) {}

template <class Input> auto erfc_t<Input>::d() const -> double {
    constexpr double two_over_root_pi =
        2.0 / constants::sqrt(constants::pi<double>());
#ifndef NDEBUG
    std::cout << "(-std::exp(" << -this->m_active.v() << " * "
              << this->m_active.v() << ") * " << two_over_root_pi << ")";
#endif
    return -std::exp(-this->m_active.v() * this->m_active.v()) *
           two_over_root_pi;
}

template <class Input> auto erfc_t<Input>::input() const -> Input const & {
    return this->m_active;
}

template <class Derived>
auto erfc(Base<Derived> const &in) -> erfc_t<const Derived> {
    return erfc_t<const Derived>(*static_cast<const Derived *>(&in));
}

template <class Input1, class Input2>
class add : public Base<add<Input1, Input2>>,
            public Bivariate<Input1, Input2, add<Input1, Input2>> {
    Input1 const m_active1;
    Input2 const m_active2;

  public:
    add(Input1 const &active1, Input2 const &active2);
    [[nodiscard]] auto d1() const -> double;
    [[nodiscard]] auto d2() const -> double;
    [[nodiscard]] auto input1() const -> Input1 const &;
    [[nodiscard]] auto input2() const -> Input2 const &;
};

template <class Input1, class Input2>
add<Input1, Input2>::add(Input1 const &active1, Input2 const &active2)
    : Base<add<Input1, Input2>>(active1.v() + active2.v()), m_active1(active1),
      m_active2(active2) {}

template <class Input1, class Input2>
auto add<Input1, Input2>::d1() const -> double {
#ifndef NDEBUG
    std::cout << "1.0";
#endif
    return 1.0;
}

template <class Input1, class Input2>
auto add<Input1, Input2>::d2() const -> double {
#ifndef NDEBUG
    std::cout << "1.0";
#endif
    return 1.0;
}

template <class Input1, class Input2>
auto add<Input1, Input2>::input1() const -> Input1 const & {
    return this->m_active1;
}

template <class Input1, class Input2>
auto add<Input1, Input2>::input2() const -> Input2 const & {
    return this->m_active2;
}

template <class Input1, class Input2>
class mul : public Base<mul<Input1, Input2>>,
            public Bivariate<Input1, Input2, mul<Input1, Input2>> {
    Input1 const m_active1;
    Input2 const m_active2;

  public:
    mul(Input1 const &active1, Input2 const &active2);
    [[nodiscard]] auto d1() const -> double;
    [[nodiscard]] auto d2() const -> double;
    [[nodiscard]] auto input1() const -> Input1 const &;
    [[nodiscard]] auto input2() const -> Input2 const &;
};

template <class Input1, class Input2>
mul<Input1, Input2>::mul(Input1 const &active1, Input2 const &active2)
    : Base<mul<Input1, Input2>>(active1.v() * active2.v()), m_active1(active1),
      m_active2(active2) {}

template <class Input1, class Input2>
auto mul<Input1, Input2>::d1() const -> double {
#ifndef NDEBUG
    std::cout << this->m_active2.v();
#endif
    return this->m_active2.v();
}

template <class Input1, class Input2>
auto mul<Input1, Input2>::d2() const -> double {
#ifndef NDEBUG
    std::cout << this->m_active1.v();
#endif
    return this->m_active1.v();
}

template <class Input1, class Input2>
auto mul<Input1, Input2>::input1() const -> Input1 const & {
    return this->m_active1;
}

template <class Input1, class Input2>
auto mul<Input1, Input2>::input2() const -> Input2 const & {
    return this->m_active2;
}

// NOTE: we could use div = mul<in1, inv<in2>> however this might lead to
// numerical errors on the valuation because in1 / in2 is not the same as in1 *
// (1.0/in2)
template <class Input1, class Input2>
class div : public Base<div<Input1, Input2>>,
            public Bivariate<Input1, Input2, div<Input1, Input2>> {
    Input1 const m_active1;
    Input2 const m_active2;

  public:
    div(Input1 const &active1, Input2 const &active2);
    [[nodiscard]] auto d1() const -> double;
    [[nodiscard]] auto d2() const -> double;
    [[nodiscard]] auto input1() const -> Input1 const &;
    [[nodiscard]] auto input2() const -> Input2 const &;
};

template <class Input1, class Input2>
div<Input1, Input2>::div(Input1 const &active1, Input2 const &active2)
    : Base<div<Input1, Input2>>(active1.v() / active2.v()), m_active1(active1),
      m_active2(active2) {}

template <class Input1, class Input2>
auto div<Input1, Input2>::d1() const -> double {
#ifndef NDEBUG
    std::cout << "(1.0 / " << this->m_active2.v() << ")";
#endif
    return 1.0 / this->m_active2.v();
}

template <class Input1, class Input2>
auto div<Input1, Input2>::d2() const -> double {
#ifndef NDEBUG
    std::cout << "(" << -this->v() << " / " << this->m_active2.v() << ")";
#endif
    // NOTE: this is when using mul<in1, inv<in2>> might be more efficient than
    // using div<in1, in2> but we choose this approach to keep evaluation values
    // the same
    return -this->v() / this->m_active2.v();
}

template <class Input1, class Input2>
auto div<Input1, Input2>::input1() const -> Input1 const & {
    return this->m_active1;
}

template <class Input1, class Input2>
auto div<Input1, Input2>::input2() const -> Input2 const & {
    return this->m_active2;
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
    return std::is_same_v<Denom, adouble_aux<N> const>;
}

} // namespace detail

template <typename... Types> struct args {};

#define adouble detail::adouble_aux<__COUNTER__>
template <int N> using adhoc = detail::adouble_aux<N>;
#define ID __COUNTER__

template <class T, class R = void> struct enable_if_type { using type = R; };

template <class T, class Enable = void> struct is_adhoc : std::false_type {};

// constant is a template useful for templetisation of functions
template <class T>
struct is_adhoc<T, typename enable_if_type<typename T::is_adhoc_tag>::type>
    : std::true_type {};

template <bool T> struct constant_type;

template <> struct constant_type<true> {
    template <int N> using type = adhoc<N>;
};

template <> struct constant_type<false> {
    template <int N> using type = double;
};

template <int N, class T>
using constant =
    typename constant_type<is_adhoc<T>::type::value>::template type<N>;

} // namespace adhoc2

#endif // ADHOC2_HPP