#ifndef CONSTANTS_TAPE_INPUT3_HPP
#define CONSTANTS_TAPE_INPUT3_HPP

#include <adhoc3.hpp>

#include <array>
#include <type_traits>

namespace adhoc3 {

namespace detail {

template <std::size_t N, class... Vars> struct get_alias {
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

template <std::size_t N> class TapeInput {
    std::array<double, N> buff{};

  public:
    constexpr static auto getalias() -> auto;

    template <std::size_t M> void set(detail::adouble_aux<M> var, double val);
};

template <std::size_t N> constexpr auto TapeInput<N>::getalias() -> auto {
    return detail::get_alias<N>::call();
}

template <std::size_t N>
template <std::size_t M>
void TapeInput<N>::set(detail::adouble_aux<M>, double val) {
    static_assert(M < N);
    this->buff[M] = val;
}

} // namespace adhoc3

#endif // CONSTANTS_TAPE_INPUT3_HPP