#ifndef ADHOC2_HPP
#define ADHOC2_HPP

#include <base.hpp>
#include <constants_constexpr.hpp>
#include <filter.hpp>
#include <utils.hpp>

#include <array>
#include <cmath>
#include <iostream>
#include <type_traits>

#ifdef CODELOGGER
#include <iostream>
#endif

namespace adhoc2 {

template <class Derived> class Val {
  protected:
    double m_value{0};

  public:
    explicit Val(double value);
    [[nodiscard]] auto v() const noexcept -> double;
};

template <class Derived> Val<Derived>::Val(double value) : m_value(value) {}

template <class Derived>
inline auto Val<Derived>::v() const noexcept -> double {
    return this->m_value;
}

template <class Input>
class exp_t : public Base<exp_t<Input>>, public Val<exp_t<Input>> {
  private:
    Input const m_active;

  public:
    explicit exp_t(Input const &active);
    [[nodiscard]] auto d() const -> double;
    [[nodiscard]] auto input() const -> Input const &;
};

template <class Input>
exp_t<Input>::exp_t(Input const &active)
    : Val<exp_t<Input>>(std::exp(active.v())), m_active(active) {}

template <class Input> auto exp_t<Input>::d() const -> double {
#ifdef CODELOGGER
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
class cos_t : public Base<cos_t<Input>>, public Val<cos_t<Input>> {
  private:
    Input const m_active;

  public:
    explicit cos_t(Input const &active);
    [[nodiscard]] auto d() const -> double;
    [[nodiscard]] auto input() const -> Input const &;
};

template <class Input>
cos_t<Input>::cos_t(Input const &active)
    : Val<cos_t<Input>>(std::cos(active.v())), m_active(active) {}

template <class Input> auto cos_t<Input>::d() const -> double {
#ifdef CODELOGGER
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
class sqrt_t : public Base<sqrt_t<Input>>, public Val<sqrt_t<Input>> {
  private:
    Input const m_active;

  public:
    explicit sqrt_t(Input const &active);
    [[nodiscard]] auto d() const -> double;
    [[nodiscard]] auto input() const -> Input const &;
};

template <class Input>
sqrt_t<Input>::sqrt_t(Input const &active)
    : Val<sqrt_t<Input>>(std::sqrt(active.v())), m_active(active) {}

template <class Input> auto sqrt_t<Input>::d() const -> double {
#ifdef CODELOGGER
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
class log_t : public Base<log_t<Input>>, public Val<log_t<Input>> {
  private:
    Input const m_active;

  public:
    explicit log_t(Input const &active);
    [[nodiscard]] auto d() const -> double;
    [[nodiscard]] auto input() const -> Input const &;
};

template <class Input>
log_t<Input>::log_t(Input const &active)
    : Val<log_t<Input>>(std::log(active.v())), m_active(active) {}

template <class Input> auto log_t<Input>::d() const -> double {
#ifdef CODELOGGER
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
class erfc_t : public Base<erfc_t<Input>>, public Val<erfc_t<Input>> {
  private:
    Input const m_active;

  public:
    explicit erfc_t(Input const &active);
    [[nodiscard]] auto d() const -> double;
    [[nodiscard]] auto input() const -> Input const &;
};

template <class Input>
erfc_t<Input>::erfc_t(Input const &active)
    : Val<erfc_t<Input>>(std::erfc(active.v())), m_active(active) {}

template <class Input> auto erfc_t<Input>::d() const -> double {
    constexpr double two_over_root_pi =
        2.0 / constexpression::sqrt(constexpression::pi<double>());
#ifdef CODELOGGER
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
class add_t : public Base<add_t<Input1, Input2>>,
              public Val<add_t<Input1, Input2>> {
    Input1 const m_active1;
    Input2 const m_active2;

  public:
    add_t(Input1 const &active1, Input2 const &active2);
    [[nodiscard]] auto d1() const -> double;
    [[nodiscard]] auto d2() const -> double;
    [[nodiscard]] auto input1() const -> Input1 const &;
    [[nodiscard]] auto input2() const -> Input2 const &;
};

template <class Input1, class Input2>
add_t<Input1, Input2>::add_t(Input1 const &active1, Input2 const &active2)
    : Val<add_t<Input1, Input2>>(active1.v() + active2.v()), m_active1(active1),
      m_active2(active2) {}

template <class Input1, class Input2>
auto add_t<Input1, Input2>::d1() const -> double {
#ifdef CODELOGGER
    std::cout << "1.0";
#endif
    return 1.0;
}

template <class Input1, class Input2>
auto add_t<Input1, Input2>::d2() const -> double {
#ifdef CODELOGGER
    std::cout << "1.0";
#endif
    return 1.0;
}

template <class Input1, class Input2>
auto add_t<Input1, Input2>::input1() const -> Input1 const & {
    return this->m_active1;
}

template <class Input1, class Input2>
auto add_t<Input1, Input2>::input2() const -> Input2 const & {
    return this->m_active2;
}

template <class Input1, class Input2>
class sub_t : public Base<sub_t<Input1, Input2>>,
              public Val<sub_t<Input1, Input2>> {
    Input1 const m_active1;
    Input2 const m_active2;

  public:
    sub_t(Input1 const &active1, Input2 const &active2);
    [[nodiscard]] auto d1() const -> double;
    [[nodiscard]] auto d2() const -> double;
    [[nodiscard]] auto input1() const -> Input1 const &;
    [[nodiscard]] auto input2() const -> Input2 const &;
};

template <class Input1, class Input2>
sub_t<Input1, Input2>::sub_t(Input1 const &active1, Input2 const &active2)
    : Val<sub_t<Input1, Input2>>(active1.v() - active2.v()), m_active1(active1),
      m_active2(active2) {}

template <class Input1, class Input2>
auto sub_t<Input1, Input2>::d1() const -> double {
#ifdef CODELOGGER
    std::cout << "1.0";
#endif
    return 1.0;
}

template <class Input1, class Input2>
auto sub_t<Input1, Input2>::d2() const -> double {
#ifdef CODELOGGER
    std::cout << "(-1.0)";
#endif
    return -1.0;
}

template <class Input1, class Input2>
auto sub_t<Input1, Input2>::input1() const -> Input1 const & {
    return this->m_active1;
}

template <class Input1, class Input2>
auto sub_t<Input1, Input2>::input2() const -> Input2 const & {
    return this->m_active2;
}

template <class Input1, class Input2>
class mul_t : public Base<mul_t<Input1, Input2>>,
              public Val<mul_t<Input1, Input2>> {
    Input1 const m_active1;
    Input2 const m_active2;

  public:
    mul_t(Input1 const &active1, Input2 const &active2);
    [[nodiscard]] auto d1() const -> double;
    [[nodiscard]] auto d2() const -> double;
    [[nodiscard]] auto input1() const -> Input1 const &;
    [[nodiscard]] auto input2() const -> Input2 const &;
};

template <class Input1, class Input2>
mul_t<Input1, Input2>::mul_t(Input1 const &active1, Input2 const &active2)
    : Val<mul_t<Input1, Input2>>(active1.v() * active2.v()), m_active1(active1),
      m_active2(active2) {}

template <class Input1, class Input2>
auto mul_t<Input1, Input2>::d1() const -> double {
#ifdef CODELOGGER
    std::cout << this->m_active2.v();
#endif
    return this->m_active2.v();
}

template <class Input1, class Input2>
auto mul_t<Input1, Input2>::d2() const -> double {
#ifdef CODELOGGER
    std::cout << this->m_active1.v();
#endif
    return this->m_active1.v();
}

template <class Input1, class Input2>
auto mul_t<Input1, Input2>::input1() const -> Input1 const & {
    return this->m_active1;
}

template <class Input1, class Input2>
auto mul_t<Input1, Input2>::input2() const -> Input2 const & {
    return this->m_active2;
}

// NOTE: we could use div_t = mul_t<in1, inv<in2>> however this might lead to
// numerical errors on the valuation because in1 / in2 is not the same as in1 *
// (1.0/in2)
template <class Input1, class Input2>
class div_t : public Base<div_t<Input1, Input2>>,
              public Val<div_t<Input1, Input2>> {
    Input1 const m_active1;
    Input2 const m_active2;

  public:
    div_t(Input1 const &active1, Input2 const &active2);
    [[nodiscard]] auto d1() const -> double;
    [[nodiscard]] auto d2() const -> double;
    [[nodiscard]] auto input1() const -> Input1 const &;
    [[nodiscard]] auto input2() const -> Input2 const &;
};

template <class Input1, class Input2>
div_t<Input1, Input2>::div_t(Input1 const &active1, Input2 const &active2)
    : Val<div_t<Input1, Input2>>(active1.v() / active2.v()), m_active1(active1),
      m_active2(active2) {}

template <class Input1, class Input2>
auto div_t<Input1, Input2>::d1() const -> double {
#ifdef CODELOGGER
    std::cout << "(1.0 / " << this->m_active2.v() << ")";
#endif
    return 1.0 / this->m_active2.v();
}

template <class Input1, class Input2>
auto div_t<Input1, Input2>::d2() const -> double {
#ifdef CODELOGGER
    std::cout << "(" << -this->v() << " / " << this->m_active2.v() << ")";
#endif
    // NOTE: this is when using mul_t<in1, inv<in2>> might be more efficient
    // than using div_t<in1, in2> but we choose this approach to keep evaluation
    // values the same
    return -this->v() / this->m_active2.v();
}

template <class Input1, class Input2>
auto div_t<Input1, Input2>::input1() const -> Input1 const & {
    return this->m_active1;
}

template <class Input1, class Input2>
auto div_t<Input1, Input2>::input2() const -> Input2 const & {
    return this->m_active2;
}

namespace detail {
template <int N>
class adouble_aux : public Base<adouble_aux<N>>, public Val<adouble_aux<N>> {
  public:
    explicit adouble_aux(double value = 0.);
};

template <int N>
inline adouble_aux<N>::adouble_aux(double value) : Val<adouble_aux<N>>(value) {}

} // namespace detail

template <typename... Types> struct args {};

#define adouble detail::adouble_aux<__COUNTER__>
template <int N> using adhoc = detail::adouble_aux<N>;
#define ID __COUNTER__

} // namespace adhoc2

#endif // ADHOC2_HPP