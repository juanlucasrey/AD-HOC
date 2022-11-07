#ifndef CONSTANTS_TAPE_INPUT_HPP
#define CONSTANTS_TAPE_INPUT_HPP

#include <adhoc2.hpp>

#include <array>
#include <type_traits>

namespace adhoc2 {

namespace detail {

template <int N, class... Vars> struct get_alias {
    constexpr static auto call() noexcept -> auto {
        return get_alias<N - 1, detail::adouble_aux<N - 1>, Vars...>::call();
    }
};

template <class... Vars> struct get_alias<0, Vars...> {
    constexpr static auto call() noexcept -> auto {
        return std::tuple<Vars...>{};
    }
};

} // namespace detail

template <int N> class TapeInput {
    std::array<double, N> buff{};

  public:
    constexpr static auto getalias() -> auto;

    template <int M> void set(detail::adouble_aux<M> var, double val);
};

template <int N> constexpr auto TapeInput<N>::getalias() -> auto {
    return detail::get_alias<N>::call();
}

template <int N>
template <int M>
void TapeInput<N>::set(detail::adouble_aux<M> var, double val) {
    static_assert(M < N);
    this->buff[M] = val;

    // to be removed
    var = detail::adouble_aux<M>(val);
}

} // namespace adhoc2

#endif // CONSTANTS_TAPE_INPUT_HPP