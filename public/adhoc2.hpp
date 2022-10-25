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

template <class Input, int N> class mul_scalar;
template <class Input> class add_scalar;
template <class Input1, class Input2> class mul_active;
template <class Input1, class Input2> class add_active;

template <class Derived> class Base {
  protected:
    double m_value{0};

  public:
    explicit Base(double value);
    [[nodiscard]] auto v() const noexcept -> double;

    template <int N>
    auto mul_scalar2(double rhs) const -> mul_scalar<const Derived, N>;
    auto operator+(double rhs) const -> add_scalar<const Derived>;
    // auto operator*(double rhs) const -> mul_scalar<const Derived>;

    template <class Derived2>
    auto operator+(Base<Derived2> const &rhs) const
        -> add_active<const Derived, const Derived2>;

    template <class Derived2>
    auto operator*(Base<Derived2> const &rhs) const
        -> mul_active<const Derived, const Derived2>;

    template <class... Denom>
    constexpr static auto depends3() noexcept -> std::size_t;

    template <class... Denom>
    inline auto backward(Denom const &...in) const noexcept
        -> std::array<double, Base<Derived>::template depends3<Denom...>()>;
};

template <class Derived> Base<Derived>::Base(double value) : m_value(value) {}

template <class Derived>
inline auto Base<Derived>::v() const noexcept -> double {
    return this->m_value;
}

template <class Derived>
template <int N>
auto Base<Derived>::mul_scalar2(double rhs) const
    -> mul_scalar<const Derived, N> {
    return {rhs, *static_cast<Derived const *>(this)};
}

template <class Derived>
auto Base<Derived>::operator+(double rhs) const -> add_scalar<const Derived> {
    return {rhs, *static_cast<Derived const *>(this)};
}

// template <class Derived>
// auto Base<Derived>::operator*(double rhs) const -> mul_scalar<const Derived>
// {
//     return {rhs, *static_cast<Derived const *>(this)};
// }

template <class Derived>
template <class Derived2>
inline auto Base<Derived>::operator+(Base<Derived2> const &rhs) const
    -> add_active<const Derived, const Derived2> {
    return {*static_cast<Derived const *>(this),
            *static_cast<Derived2 const *>(&rhs)};
}

template <class Derived>
template <class Derived2>
inline auto Base<Derived>::operator*(Base<Derived2> const &rhs) const
    -> mul_active<const Derived, const Derived2> {
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

namespace detail {

template <typename T, typename... Types>
constexpr bool are_types_unique_v =
    (!std::is_same_v<T, Types> && ...) && are_types_unique_v<Types...>;
template <typename T> inline constexpr bool are_types_unique_v<T> = true;

} // namespace detail

template <class Derived>
template <class... Denom>
inline auto Base<Derived>::backward(Denom const &.../* in */) const noexcept
    -> std::array<double, Base<Derived>::template depends3<Denom...>()> {

    // let's avoid multiplying the same thing multiple times
    static_assert(detail::are_types_unique_v<Denom...>);
    std::array<double, Base<Derived>::template depends3<Denom...>()> res;
    res[0] = 1.0;

    return res;
}

template <class Input>
class add_scalar : public Base<add_scalar<Input>>,
                   public Univariate<Input, add_scalar<Input>> {
    Input const &m_active;

  public:
    add_scalar(double scalar, Input const &active);
    template <class Denom> auto d(Denom const &in) const noexcept -> double;

    template <class... Denom>
    auto dmany(Denom const &...in) const noexcept
        -> std::array<double, add_scalar<Input>::template depends3<Denom...>()>;
};

template <class Input>
add_scalar<Input>::add_scalar(double scalar, Input const &active)
    : Base<add_scalar<Input>>(scalar + active.v()), m_active(active) {}

template <class Input>
template <class Denom>
inline auto add_scalar<Input>::d(Denom const &in) const noexcept -> double {
    if constexpr (Univariate<Input,
                             add_scalar<Input>>::template depends<Denom>()) {
        return this->m_active.d(in);
    } else {
        return 0.;
    }
}

template <class Input>
template <class... Denom>
inline auto add_scalar<Input>::dmany(Denom const &...in) const noexcept
    -> std::array<double, add_scalar<Input>::template depends3<Denom...>()> {
    return this->m_active.dmany(in...);
}

template <class Input, int N>
class mul_scalar : public Base<mul_scalar<Input, N>>,
                   public Univariate<Input, mul_scalar<Input, N>> {
    double m_scalar{0.};
    Input const &m_active;

  public:
    mul_scalar(double scalar, Input const &active);
    template <class Denom> auto d(Denom const &in) const noexcept -> double;

    template <class... Denom>
    auto dmany(Denom const &...in) const noexcept
        -> std::array<double,
                      mul_scalar<Input, N>::template depends3<Denom...>()>;
};

template <class Input, int N>
mul_scalar<Input, N>::mul_scalar(double scalar, Input const &active)
    : Base<mul_scalar<Input, N>>(scalar * active.v()), m_scalar(scalar),
      m_active(active) {}

template <class Input, int N>
template <class Denom>
inline auto mul_scalar<Input, N>::d(Denom const &in) const noexcept -> double {
    if constexpr (Univariate<Input,
                             mul_scalar<Input, N>>::template depends<Denom>()) {
        return this->m_scalar * this->m_active.d(in);
    } else {
        return 0.;
    }
}

template <class Input, int N>
template <class... Denom>
inline auto mul_scalar<Input, N>::dmany(Denom const &...in) const noexcept
    -> std::array<double, mul_scalar<Input, N>::template depends3<Denom...>()> {
    constexpr std::size_t M =
        mul_scalar<Input, N>::template depends3<Denom...>();
    auto res = this->m_active.dmany(in...);

    // should be vectorised by compiler
    for (auto &r : res) {
        r *= this->m_scalar;
    }
    return res;
}

template <typename... Types> struct args {};

template <class Input1, class Input2>
class add_active
    : public Base<add_active<Input1, Input2>>,
      public Bivariate<Input1, Input2, add_active<Input1, Input2>> {
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
class mul_active
    : public Base<mul_active<Input1, Input2>>,
      public Bivariate<Input1, Input2, mul_active<Input1, Input2>> {
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
        double, mul_active<Input1, Input2>::template depends3<Denom...>()>;
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

namespace detail {
template <int N>
class adouble_aux : public Base<adouble_aux<N>>,
                    public Univariate<adouble_aux<N>, adouble_aux<N>> {
  public:
    explicit adouble_aux(double value);

    template <class Denom> constexpr static auto depends() noexcept -> bool;

    template <class Denom>
    constexpr auto d(Denom const &in) const noexcept -> double;

    template <class Denom1, class Denom2>
    constexpr auto d2(Denom1 const &in1, Denom2 const &in2) const noexcept
        -> std::array<double,
                      adouble_aux<N>::template depends3<Denom1, Denom2>()>;

    template <class... Denom>
    constexpr auto dmany(Denom const &...in) const noexcept
        -> std::array<double, adouble_aux<N>::template depends3<Denom...>()>;

    template <class... Denom>
    constexpr auto dinterface(Denom const &...in) const noexcept
        -> std::array<double, sizeof...(Denom)>;
};

template <int N>
inline adouble_aux<N>::adouble_aux(double value)
    : Base<adouble_aux<N>>(value) {}

template <int N>
template <class Denom>
inline constexpr auto adouble_aux<N>::depends() noexcept -> bool {
    return std::is_same_v<Denom, adouble_aux<N>>;
}

template <int N>
template <class Denom>
inline constexpr auto adouble_aux<N>::d(Denom const & /* in */) const noexcept
    -> double {
    return static_cast<double>(this->depends<Denom>());
}

template <int N>
template <class Denom1, class Denom2>
inline constexpr auto
adouble_aux<N>::d2(Denom1 const & /* in1 */,
                   Denom2 const & /* in2 */) const noexcept
    -> std::array<double, adouble_aux<N>::template depends3<Denom1, Denom2>()> {
    constexpr std::size_t M =
        adouble_aux<N>::template depends3<Denom1, Denom2>();
    std::array<double, M> ret{};
    for (std::size_t i = 0; i < M; i++) {
        ret[i] = 1.0;
    }
    return ret;
}

template <int N>
template <class... Denom>
inline constexpr auto
adouble_aux<N>::dmany(Denom const &.../* in */) const noexcept
    -> std::array<double, adouble_aux<N>::template depends3<Denom...>()> {
    constexpr std::size_t M = adouble_aux<N>::template depends3<Denom...>();
    std::array<double, M> ret{};
    for (std::size_t i = 0; i < M; i++) {
        ret[i] = 1.0;
    }
    return ret;
}

template <int N>
template <class... Denom>
inline constexpr auto
adouble_aux<N>::dinterface(Denom const &.../* in */) const noexcept
    -> std::array<double, sizeof...(Denom)> {
    return {static_cast<double>(adouble_aux<N>::template depends<Denom>())...};
}

} // namespace detail

#define adouble detail::adouble_aux<__COUNTER__>
#define mulscalar(in1, in2) in1.mul_scalar2<__COUNTER__>(in2)

} // namespace adhoc2

#endif // ADHOC2_HPP