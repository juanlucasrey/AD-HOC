#ifndef TAPE_INTERMEDIATE_HPP
#define TAPE_INTERMEDIATE_HPP

#include <adhoc3.hpp>
#include <tape_size_fwd.hpp>
#include <utils.hpp>

#include <array>
#include <type_traits>

namespace adhoc3 {

namespace detail {

template <class... IntermediateTypes> class TapeIntermediate {
    std::array<double, sizeof...(IntermediateTypes)> buff{};

  public:
    constexpr static auto getalias() -> auto;
    template <class Derived> void inline set(Base<Derived> var, double val);
    template <class Derived> auto inline get(Base<Derived> var) const -> double;
};

template <class... IntermediateTypes>
constexpr auto TapeIntermediate<IntermediateTypes...>::getalias() -> auto {
    return std::tuple<IntermediateTypes...>{};
}

template <class... IntermediateTypes>
template <class Derived>
void TapeIntermediate<IntermediateTypes...>::set(Base<Derived>, double val) {
    constexpr auto idx = adhoc2::idx_type2<Derived, IntermediateTypes...>();
    this->buff[idx] = val;
}

template <class... IntermediateTypes>
template <class Derived>
auto TapeIntermediate<IntermediateTypes...>::get(Base<Derived>) const
    -> double {
    constexpr auto idx = adhoc2::idx_type2<Derived, IntermediateTypes...>();
    return this->buff[idx];
}

template <class... IntermediateTypes>
auto inline CreateTapeIntermediate_aux(std::tuple<IntermediateTypes...>) {
    return detail::TapeIntermediate<IntermediateTypes...>{};
}

template <std::size_t... Idxs>
auto inline CreateTapeInitial_aux(std::index_sequence<Idxs...> const &) {
    return detail::TapeIntermediate<detail::adouble_aux<Idxs>...>{};
}

} // namespace detail

template <class Output> auto inline CreateTapeIntermediate(Output & /* out */) {
    constexpr auto f = adhoc3::fwd_calc_order_t(Output{});
    return detail::CreateTapeIntermediate_aux(f);
}

template <std::size_t N> auto inline CreateTapeInitial() {
    return detail::CreateTapeInitial_aux(std::make_index_sequence<N>{});
}

} // namespace adhoc3

#endif // TAPE_INTERMEDIATE_HPP