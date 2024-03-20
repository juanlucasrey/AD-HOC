#ifndef ADHOC3_PARTITION_PARTITION_FUNCTION_HPP
#define ADHOC3_PARTITION_PARTITION_FUNCTION_HPP

// #include <cstddef>

namespace adhoc3 {

namespace detail {
constexpr int partition_function_aux(int n, int k) {
    if (k == 0) {
        return 0;
    }
    if (n == 0) {
        return 1;
    }
    if (n < 0) {
        return 0;
    }

    return partition_function_aux(n, k - 1) + partition_function_aux(n - k, k);
}
} // namespace detail

constexpr auto partition_function(int n) -> int {
    // convention
    if (n == 0) {
        return 1;
    }
    return detail::partition_function_aux(n, n);
}

} // namespace adhoc3

#endif // ADHOC3_PARTITION_PARTITION_FUNCTION_HPP
