#ifndef CONSTANTS_TYPE_HPP
#define CONSTANTS_TYPE_HPP

#include <base.hpp>
#include <constants_constexpr.hpp>

namespace adhoc2::constants {

template <class Derived> struct ConstBase {
    operator double() const { return Derived::v(); }
};

struct Pi : public Base<Pi>, public ConstBase<Pi> {
    constexpr auto static v() -> double { return detail::pi<double>(); }
};

template <std::size_t N>
struct Const : public Base<Const<N>>, public ConstBase<Const<N>> {
    constexpr auto static v() -> double { return static_cast<double>(N); }
};

template <class N, class D>
struct Frac : public Base<Frac<N, D>>, public ConstBase<Frac<N, D>> {
    constexpr auto static v() -> double { return N::v() / D::v(); }
};

template <class V>
struct Sqrt : public Base<Sqrt<V>>, public ConstBase<Sqrt<V>> {
    constexpr auto static v() -> double { return detail::sqrt(V::v()); }
};

template <class V>
struct Minus : public Base<Minus<V>>, public ConstBase<Minus<V>> {
    constexpr auto static v() -> double { return -V::v(); }
};

#if __cplusplus >= 202002L
namespace detail {
template <typename T> struct AsTemplateArg {
    T buffer{};
    constexpr AsTemplateArg(T t) : buffer(t) {}
    constexpr operator T() const { return this->buffer; }
};
} // namespace detail

template <typename T, detail::AsTemplateArg<T> F>
struct C : public Base<C<T, F>> {
    constexpr operator T() const { return F; }
    constexpr auto static v() -> T { return F; }
};

template <detail::AsTemplateArg<float> F> using CF = C<float, F>;
template <detail::AsTemplateArg<double> F> using CD = C<double, F>;
template <detail::AsTemplateArg<long double> F> using CLD = C<long double, F>;
#endif

} // namespace adhoc2::constants

#endif // CONSTANTS_TYPE_HPP
