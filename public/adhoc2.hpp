#ifndef ADHOC2_HPP
#define ADHOC2_HPP

#include <array>
#include <type_traits>

namespace adhoc2 {

template <class Input> class Univariate {
  public:
    template <class Denom> constexpr static auto depends() noexcept -> bool;
    template <class... Denom>
    constexpr static auto depends2() noexcept -> std::size_t;
};

template <class Input>
template <class Denom>
inline constexpr auto Univariate<Input>::depends() noexcept -> bool {
    return Input::template depends<Denom>();
}

template <class Input>
template <class... Denom>
inline constexpr auto Univariate<Input>::depends2() noexcept -> std::size_t {
    return (Univariate<Input>::depends<Denom>() + ...);
}

template <class Input1, class Input2> class Bivariate {
  public:
    template <class Denom> constexpr static auto depends() noexcept -> bool;

    template <class... Denom>
    constexpr static auto depends2() noexcept -> std::size_t;
};

template <class Input1, class Input2>
template <class Denom>
inline constexpr auto Bivariate<Input1, Input2>::depends() noexcept -> bool {
    return Input1::template depends<Denom>() ||
           Input2::template depends<Denom>();
}

template <class Input1, class Input2>
template <class... Denom>
inline constexpr auto Bivariate<Input1, Input2>::depends2() noexcept
    -> std::size_t {
    return (Bivariate<Input1, Input2>::depends<Denom>() + ...);
}

template <class Input> class mul_scalar;
template <class Input> class add_scalar;
template <class Input1, class Input2> class mul_active;
template <class Input1, class Input2> class add_active;

template <class Derived> class Base {
  protected:
    double m_value{0};

  public:
    explicit Base(double value);
    [[nodiscard]] auto v() const noexcept -> double;
    auto operator+(double rhs) const -> add_scalar<const Derived>;
    auto operator*(double rhs) const -> mul_scalar<const Derived>;

    template <class Arg>
    auto operator+(Arg const &rhs) const
        -> add_active<const Derived, const Arg>;

    template <class Arg>
    auto operator*(Arg const &rhs) const
        -> mul_active<const Derived, const Arg>;
};

template <class Derived> Base<Derived>::Base(double value) : m_value(value) {}

template <class Derived>
inline auto Base<Derived>::v() const noexcept -> double {
    return this->m_value;
}

template <class Derived>
auto Base<Derived>::operator+(double rhs) const -> add_scalar<const Derived> {
    return {rhs, *static_cast<Derived const *>(this)};
}

template <class Derived>
auto Base<Derived>::operator*(double rhs) const -> mul_scalar<const Derived> {
    return {rhs, *static_cast<Derived const *>(this)};
}

template <class Derived>
template <class Arg>
inline auto Base<Derived>::operator+(Arg const &rhs) const
    -> add_active<const Derived, const Arg> {
    return {*static_cast<Derived const *>(this), rhs};
}

template <class Derived>
template <class Arg>
inline auto Base<Derived>::operator*(Arg const &rhs) const
    -> mul_active<const Derived, const Arg> {
    return {*static_cast<Derived const *>(this), rhs};
}

template <class Input>
class add_scalar : public Base<add_scalar<Input>>, public Univariate<Input> {
    Input const &m_active;

  public:
    add_scalar(double scalar, Input const &active);
    template <class Denom> auto d(Denom const &in) const noexcept -> double;
};

template <class Input>
add_scalar<Input>::add_scalar(double scalar, Input const &active)
    : Base<add_scalar<Input>>(scalar + active.v()), m_active(active) {}

template <class Input>
template <class Denom>
inline auto add_scalar<Input>::d(Denom const &in) const noexcept -> double {
    if constexpr (Univariate<Input>::template depends<Denom>()) {
        return this->m_active.d(in);
    } else {
        return 0.;
    }
}

template <class Input>
class mul_scalar : public Base<mul_scalar<Input>>, public Univariate<Input> {
    double m_scalar{0.};
    Input const &m_active;

  public:
    mul_scalar(double scalar, Input const &active);
    template <class Denom> auto d(Denom const &in) const noexcept -> double;
};

template <class Input>
mul_scalar<Input>::mul_scalar(double scalar, Input const &active)
    : Base<mul_scalar<Input>>(scalar * active.v()), m_scalar(scalar),
      m_active(active) {}

template <class Input>
template <class Denom>
inline auto mul_scalar<Input>::d(Denom const &in) const noexcept -> double {
    if constexpr (Univariate<Input>::template depends<Denom>()) {
        return this->m_scalar * this->m_active.d(in);
    } else {
        return 0.;
    }
}

template <class Input1, class Input2>
class add_active : public Base<add_active<Input1, Input2>>,
                   public Bivariate<Input1, Input2> {
    Input1 const &m_active1;
    Input2 const &m_active2;

  public:
    add_active(Input1 const &active1, Input2 const &active2);
    template <class Denom>
    [[nodiscard]] auto d(Denom const &in) const noexcept -> double;
    template <class Denom1, class Denom2>
    auto d2(Denom1 const &in1, Denom2 const &in2) const noexcept
        -> std::array<double, 2>;
};

template <class Input1, class Input2>
add_active<Input1, Input2>::add_active(Input1 const &active1,
                                       Input2 const &active2)
    : Base<add_active<Input1, Input2>>(active1.v() + active2.v()),
      m_active1(active1), m_active2(active2) {}

template <class Input1, class Input2>
template <class Denom>
inline auto add_active<Input1, Input2>::d(Denom const &in) const noexcept
    -> double {
    if constexpr (Input1::template depends<Denom>()) {
        if constexpr (Input2::template depends<Denom>()) {
            return this->m_active1.d(in) + this->m_active2.d(in);
        } else {
            return this->m_active1.d(in);
        }
    } else {
        if constexpr (Input2::template depends<Denom>()) {
            return this->m_active2.d(in);
        } else {
            return 0.;
        }
    }
}

template <class Input1, class Input2>
template <class Denom1, class Denom2>
inline auto add_active<Input1, Input2>::d2(Denom1 const &in1,
                                           Denom2 const &in2) const noexcept
    -> std::array<double, 2> {
    auto res1 = this->m_active1.d2(in1, in2);
    auto res2 = this->m_active2.d2(in1, in2);
    std::array<double, 2> res{};
    for (std::size_t i = 0; i < 2; ++i) {
        res[i] = res1[i] + res2[i];
    }
    return res;
}

template <class Input1, class Input2>
class mul_active : public Base<mul_active<Input1, Input2>>,
                   public Bivariate<Input1, Input2> {
    Input1 const &m_active1;
    Input2 const &m_active2;

  public:
    mul_active(Input1 const &active1, Input2 const &active2);
    template <class Denom> auto d(Denom const &in) const noexcept -> double;
    template <class Denom1, class Denom2>
    auto d2(Denom1 const &in1, Denom2 const &in2) const noexcept
        -> std::array<double, 2>;

    template <class... Denom>
    constexpr auto dmany(Denom const &...in) const noexcept -> std::array<
        double, mul_active<Input1, Input2>::template depends2<Denom...>()>;
};

template <class Input1, class Input2>
mul_active<Input1, Input2>::mul_active(Input1 const &active1,
                                       Input2 const &active2)
    : Base<mul_active<Input1, Input2>>(active1.v() * active2.v()),
      m_active1(active1), m_active2(active2) {}

template <class Input1, class Input2>
template <class Denom>
inline auto mul_active<Input1, Input2>::d(Denom const &in) const noexcept
    -> double {
    if constexpr (Input1::template depends<Denom>()) {
        if constexpr (Input2::template depends<Denom>()) {
            return this->m_active1.d(in) * this->m_active2.v() +
                   this->m_active2.d(in) * this->m_active1.v();
        } else {
            return this->m_active1.d(in) * this->m_active2.v();
        }
    } else {
        if constexpr (Input2::template depends<Denom>()) {
            return this->m_active2.d(in) * this->m_active1.v();
        } else {
            return 0.;
        }
    }
}

template <class Input1, class Input2>
template <class Denom1, class Denom2>
inline auto mul_active<Input1, Input2>::d2(Denom1 const &in1,
                                           Denom2 const &in2) const noexcept
    -> std::array<double, 2> {
    auto res1 = this->m_active1.d2(in1, in2);
    auto res2 = this->m_active2.d2(in1, in2);
    double v1 = this->m_active1.v();
    double v2 = this->m_active2.v();
    std::array<double, 2> res{};
    for (std::size_t i = 0; i < 2; ++i) {
        res[i] = res1[i] * v2 + res2[i] * v1;
    }
    return res;
}

template <class Input1, class Input2>
template <class... Denom>
inline constexpr auto
mul_active<Input1, Input2>::dmany(Denom const &.../* in */) const noexcept
    -> std::array<double,
                  mul_active<Input1, Input2>::template depends2<Denom...>()> {
    constexpr std::size_t M =
        mul_active<Input1, Input2>::template depends2<Denom...>();
    std::array<double, M> ret{};
    for (std::size_t i = 0; i < M; i++) {
        ret[i] = 1.0;
    }
    return ret;
}

namespace detail {
template <int N>
class adouble : public Base<adouble<N>>, public Univariate<adouble<N>> {
  public:
    explicit adouble(double value);

    template <class Denom> constexpr static auto depends() noexcept -> bool;

    template <class Denom>
    constexpr auto d(Denom const &in) const noexcept -> double;

    template <class Denom1, class Denom2>
    constexpr auto d2(Denom1 const &in1, Denom2 const &in2) const noexcept
        -> std::array<double, adouble<N>::template depends2<Denom1, Denom2>()>;

    template <class... Denom>
    constexpr auto dmany(Denom const &...in) const noexcept
        -> std::array<double, adouble<N>::template depends2<Denom...>()>;

    template <class... Denom>
    constexpr auto dinterface(Denom const &...in) const noexcept
        -> std::array<double, sizeof...(Denom)>;
};

template <int N>
inline adouble<N>::adouble(double value) : Base<adouble<N>>(value) {}

template <int N>
template <class Denom>
inline constexpr auto adouble<N>::depends() noexcept -> bool {
    return std::is_same_v<Denom, adouble<N>>;
}

template <int N>
template <class Denom>
inline constexpr auto adouble<N>::d(Denom const & /* in */) const noexcept
    -> double {
    return static_cast<double>(this->depends<Denom>());
}

template <int N>
template <class Denom1, class Denom2>
inline constexpr auto adouble<N>::d2(Denom1 const & /* in1 */,
                                     Denom2 const & /* in2 */) const noexcept
    -> std::array<double, adouble<N>::template depends2<Denom1, Denom2>()> {
    constexpr std::size_t M = adouble<N>::template depends2<Denom1, Denom2>();
    std::array<double, M> ret{};
    for (std::size_t i = 0; i < M; i++) {
        ret[i] = 1.0;
    }
    return ret;
}

template <int N>
template <class... Denom>
inline constexpr auto adouble<N>::dmany(Denom const &.../* in */) const noexcept
    -> std::array<double, adouble<N>::template depends2<Denom...>()> {
    constexpr std::size_t M = adouble<N>::template depends2<Denom...>();
    std::array<double, M> ret{};
    for (std::size_t i = 0; i < M; i++) {
        ret[i] = 1.0;
    }
    return ret;
}

template <int N>
template <class... Denom>
inline constexpr auto
adouble<N>::dinterface(Denom const &.../* in */) const noexcept
    -> std::array<double, sizeof...(Denom)> {
    return {static_cast<double>(adouble<N>::template depends<Denom>())...};
}

} // namespace detail

#define adouble detail::adouble<__COUNTER__>

} // namespace adhoc2

#endif // ADHOC2_HPP