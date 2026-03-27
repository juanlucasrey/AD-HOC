#ifndef BMC_TOOLS_ARRAY_OPERATORS
#define BMC_TOOLS_ARRAY_OPERATORS

#include <array>

namespace adhoc {

template<class Float, std::size_t N>
auto
operator*(Float lhs, std::array<Float, N> const& rhs) -> std::array<Float, N>
{
    std::array<Float, N> res;
    for (auto i = 0; i < N; ++i) {
        res[i] = lhs * rhs[i];
    }
    return res;
}

template<class Float, std::size_t N>
auto
operator*(std::array<Float, N> const& lhs, Float rhs) -> std::array<Float, N>
{
    std::array<Float, N> res;
    for (auto i = 0; i < N; ++i) {
        res[i] = rhs * lhs[i];
    }
    return res;
}

template<class Float, std::size_t N>
auto
operator*=(std::array<Float, N>& lhs, Float rhs) -> std::array<Float, N>&
{
    for (auto i = 0; i < N; ++i) {
        lhs[i] = rhs * lhs[i];
    }
    return lhs;
}

template<class Float, std::size_t N>
auto
operator+(std::array<Float, N> const& lhs, std::array<Float, N> const& rhs) -> std::array<Float, N>
{
    std::array<Float, N> res;
    for (auto i = 0; i < N; ++i) {
        res[i] = lhs[i] + rhs[i];
    }
    return res;
}

template<class Float, std::size_t NL, std::size_t NR>
auto
operator+=(std::array<Float, NL>& lhs, std::array<Float, NR> const& rhs) -> std::array<Float, NL>&
{
    std::size_t constexpr Nmin = std::min(NL, NR);
    for (auto i = 0; i < Nmin; ++i) {
        lhs[i] += rhs[i];
    }
    return lhs;
}

template<class Float, std::size_t NL, std::size_t NR>
auto
operator-=(std::array<Float, NL>& lhs, std::array<Float, NR> const& rhs) -> std::array<Float, NL>&
{
    std::size_t constexpr Nmin = std::min(NL, NR);
    for (auto i = 0; i < Nmin; ++i) {
        lhs[i] -= rhs[i];
    }
    return lhs;
}

template<class Float, std::size_t NL, std::size_t NR>
auto
operator*=(std::array<Float, NL>& lhs, std::array<Float, NR> const& rhs) -> std::array<Float, NL>&
{
    std::size_t constexpr Nmin = std::min(NL, NR);
    for (auto i = 0; i < Nmin; ++i) {
        lhs[i] *= rhs[i];
    }
    return lhs;
}

template<class Float, std::size_t NL, std::size_t NR>
auto
operator/=(std::array<Float, NL>& lhs, std::array<Float, NR> const& rhs) -> std::array<Float, NL>&
{
    std::size_t constexpr Nmin = std::min(NL, NR);
    for (auto i = 0; i < Nmin; ++i) {
        lhs[i] /= rhs[i];
    }
    return lhs;
}

namespace detail {

template<std::size_t I = 0, class Float, std::size_t N, class InputIt>
void
copy_aux(std::array<Float, N>& dest, InputIt& input)
{
    dest[I] = *input;
    ++input;
    if constexpr ((I + 1) < N) {
        copy_aux<I + 1>(dest, input);
    }
}

template<std::size_t I, class Float, std::size_t N, class InputIt>
void
reverse_copy_aux(std::array<Float, N>& dest, InputIt& input)
{
    dest[I] = *(--input);
    if constexpr (I > 0) {
        reverse_copy_aux<I - 1>(dest, input);
    }
}

template<std::size_t I = 0, class Float, std::size_t NL, std::size_t NR>
void
dot_aux(Float& res, std::array<Float, NL> const& lhs, std::array<Float, NR> const& rhs)
{
    res += lhs[I] * rhs[I];
    std::size_t constexpr Nmin = std::min(NL, NR);
    if constexpr ((I + 1) < Nmin) {
        dot_aux<I + 1>(res, lhs, rhs);
    }
}

template<std::size_t S, std::size_t I = 0, class Float, std::size_t N>
void
sum_aux(Float& res, std::array<Float, N> const& in)
{
    res += in[I];
    if constexpr ((I + 1) < S) {
        sum_aux<S, I + 1>(res, in);
    }
}

} // namespace detail

template<class Float, std::size_t N, class InputIt>
void
copy(std::array<Float, N>& dest, InputIt& input)
{
    detail::copy_aux(dest, input);
}

template<class Float, std::size_t N, class InputIt>
void
reverse_copy(std::array<Float, N>& dest, InputIt& input)
{
    detail::reverse_copy_aux<N - 1>(dest, input);
}

template<class Float, std::size_t NL, std::size_t NR>
auto
dot(std::array<Float, NL> const& lhs, std::array<Float, NR> const& rhs) -> Float
{
    Float res = 0;
    detail::dot_aux(res, lhs, rhs);
    return res;
}

template<std::size_t S, class Float, std::size_t N>
auto
sum(std::array<Float, N> const& in) -> Float
{
    Float res = 0;
    detail::sum_aux<S>(res, in);
    return res;
}

template<std::size_t M, typename T, std::size_t N>
auto
tail(std::array<T, N> const& arr) noexcept -> std::array<T, M>
{
    static_assert(M <= N, "Cannot extract more elements than array contains");

    return [&arr]<std::size_t... Idx>(std::index_sequence<Idx...>) {
        constexpr std::size_t offset = N - M;
        return std::array<T, M>{ { arr[offset + Idx]... } };
    }(std::make_index_sequence<M>{});
}

} // namespace adhoc

#endif // BMC_TOOLS_ARRAY_OPERATORS
