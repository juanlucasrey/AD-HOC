#ifndef ADHOC3_PARTITION_INTEGER_PARTITION_HPP
#define ADHOC3_PARTITION_INTEGER_PARTITION_HPP

#include <array>
#include <cstddef>

namespace adhoc3 {

template <std::size_t N> constexpr auto FirstPartition() {
    std::array<std::size_t, N> arr{0};
    arr.back() = 1;
    return arr;
}

template <std::size_t N> constexpr auto LastPartition() {
    std::array<std::size_t, N> arr{0};
    arr.front() = N;
    return arr;
}

template <std::size_t N>
constexpr auto NextPartition(std::array<std::size_t, N> const &prev) {
    if (prev == LastPartition<N>()) {
        return prev;
    }

    std::array<std::size_t, N> arr = prev;

    // Find the rightmost non-one value in arr. Also, update the
    // rem_val so that we know how much value can be accommodated
    std::size_t rem_val = arr.front();
    arr.front() = 0;
    std::size_t k = 1;
    while (arr[k] == 0) {
        k++;
    }

    arr[k]--;
    rem_val += (k + 1);

    // If rem_val is more, then the sorted order is violated. Divide
    // rem_val in different values of size p[k] and copy these values at
    // different positions after p[k]
    if (rem_val > k) {
        // std::div will be constexpr in C++23 so we use a manual version in
        // the meantime
        std::size_t const div = rem_val / k;
        arr[k - 1] += div;
        rem_val -= k * div;
    }

    if (rem_val) {
        arr[rem_val - 1]++;
    }

    return arr;
}

template <std::size_t N>
constexpr auto PrevPartition(std::array<std::size_t, N> const &prev) {
    if (prev == FirstPartition<N>()) {
        return prev;
    }

    std::array<std::size_t, N> arr = prev;

    // Find the rightmost non-one value in arr. Also, update the
    // rem_val so that we know how much value can be accommodated
    std::size_t rem_val = arr.front();
    arr.front() = 0;
    std::size_t k = 1;
    while (k >= rem_val) {
        rem_val += arr[k] * (k + 1);
        arr[k] = 0;
        k++;
    }

    arr[k]++;
    rem_val -= (k + 1);

    arr.front() = rem_val;
    return arr;
}

} // namespace adhoc3

#endif // ADHOC3_PARTITION_INTEGER_PARTITION_HPP
