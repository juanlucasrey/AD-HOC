#ifndef ADHOC3_PARTITION_MULTINOMIAL_COEFFICIENT_INDEX_SEQUENCE_HPP
#define ADHOC3_PARTITION_MULTINOMIAL_COEFFICIENT_INDEX_SEQUENCE_HPP

#include "multinomial_coefficient.hpp"

#include <array>
#include <cstddef>

namespace adhoc3 {

namespace detail {

template <std::size_t Bins, std::size_t Balls, std::size_t... I>
constexpr auto FirstMultinomialIS_aux(std::index_sequence<I...> /* i */) {
    constexpr auto first = FirstMultinomial<Bins, Balls>();
    return std::index_sequence<first[I]...>{};
}

} // namespace detail

template <std::size_t Bins, std::size_t Balls>
constexpr auto FirstMultinomialIS() {
    return detail::FirstMultinomialIS_aux<Bins, Balls>(
        std::make_index_sequence<Bins>{});
}

namespace detail {

template <std::size_t Bins, std::size_t Balls, std::size_t... I>
constexpr auto LastMultinomialIS_aux(std::index_sequence<I...> /* i */) {
    constexpr auto last = LastMultinomial<Bins, Balls>();
    return std::index_sequence<last[I]...>{};
}

} // namespace detail

template <std::size_t Bins, std::size_t Balls>
constexpr auto LastMultinomialIS() {
    return detail::LastMultinomialIS_aux<Bins, Balls>(
        std::make_index_sequence<Bins>{});
}

namespace detail {

template <std::size_t... I, std::size_t... I2>
constexpr auto NextMultinomialIS_aux(std::index_sequence<I...> /* i */,
                                     std::index_sequence<I2...> /* i2 */) {
    constexpr std::array<std::size_t, sizeof...(I)> temp{I...};
    constexpr auto next = NextMultinomial(temp);
    return std::index_sequence<next[I2]...>{};
}

} // namespace detail

template <std::size_t... I>
constexpr auto NextMultinomialIS(std::index_sequence<I...> i) {
    return detail::NextMultinomialIS_aux(
        i, std::make_index_sequence<sizeof...(I)>{});
}

} // namespace adhoc3

#endif // ADHOC3_PARTITION_MULTINOMIAL_COEFFICIENT_INDEX_SEQUENCE_HPP
