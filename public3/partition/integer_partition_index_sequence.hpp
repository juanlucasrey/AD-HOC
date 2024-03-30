#ifndef ADHOC3_PARTITION_INTEGER_PARTITION_INDEX_SEQUENCE_HPP
#define ADHOC3_PARTITION_INTEGER_PARTITION_INDEX_SEQUENCE_HPP

#include "integer_partition.hpp"

#include <array>
#include <cstddef>

namespace adhoc3 {

namespace detail {

template <std::size_t N, std::size_t... I>
constexpr auto FirstPartitionIS_aux(std::index_sequence<I...> /* i */) {
    constexpr auto first = FirstPartition<N>();
    return std::index_sequence<first[I]...>{};
}

} // namespace detail

template <std::size_t N> constexpr auto FirstPartitionIS() {
    return detail::FirstPartitionIS_aux<N>(std::make_index_sequence<N>{});
}

namespace detail {

template <std::size_t N, std::size_t... I>
constexpr auto LastPartitionIS_aux(std::index_sequence<I...> /* i */) {
    constexpr auto last = LastPartition<N>();
    return std::index_sequence<last[I]...>{};
}

} // namespace detail

template <std::size_t N> constexpr auto LastPartitionIS() {
    return detail::LastPartitionIS_aux<N>(std::make_index_sequence<N>{});
}

namespace detail {

template <std::size_t... I, std::size_t... I2>
constexpr auto NextPartitionIS_aux(std::index_sequence<I...> /* i */,
                                   std::index_sequence<I2...> /* i2 */) {
    constexpr std::array<std::size_t, sizeof...(I)> temp{I...};
    constexpr auto next = NextPartition(temp);
    return std::index_sequence<next[I2]...>{};
}

} // namespace detail

template <std::size_t... I>
constexpr auto NextPartitionIS(std::index_sequence<I...> i) {
    return detail::NextPartitionIS_aux(
        i, std::make_index_sequence<sizeof...(I)>{});
}

} // namespace adhoc3

#endif // ADHOC3_PARTITION_INTEGER_PARTITION_INDEX_SEQUENCE_HPP
