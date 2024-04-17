#ifndef ADHOC3_PARTITION_MULTINOMIAL_COEFFICIENT_HPP
#define ADHOC3_PARTITION_MULTINOMIAL_COEFFICIENT_HPP

#include <array>
#include <cstddef>

namespace adhoc3 {

template <std::size_t Bins, std::size_t Balls>
constexpr auto FirstMultinomial() {
    std::array<std::size_t, Bins> arr{0};
    arr.front() = Balls;
    return arr;
}

template <std::size_t Bins, std::size_t Balls>
constexpr auto LastMultinomial() {
    std::array<std::size_t, Bins> arr{0};
    arr.back() = Balls;
    return arr;
}

template <std::size_t Bins>
constexpr auto NextMultinomial(std::array<std::size_t, Bins> const &prev) {
    // check the previous mutlinomial coefficient is not the last
    std::array<std::size_t, Bins> arrlast{0};
    arrlast.back() = prev.back();
    if (arrlast == prev) {
        return prev;
    }

    std::array<std::size_t, Bins> arr = prev;

    std::size_t rem_val = arr.back();
    arr.back() = 0;
    std::size_t k = static_cast<std::size_t>(Bins) - 2;
    while (arr[k] == 0 && k > 0) {
        k--;
    }
    arr[k]--;
    rem_val++;
    arr[k + 1] = rem_val;
    return arr;
}

template <std::size_t Order>
constexpr auto factorialarr() -> std::array<std::size_t, Order + 1> {
    std::array<std::size_t, Order + 1> res{};
    res[0] = 1;

    std::size_t k = 1;
    while (k <= Order) {
        res[k] = res[k - 1] * k;
        k++;
    }

    return res;
}

template <std::size_t Order, std::size_t Bins, std::size_t FactorialInputs>
constexpr auto
MultinomialCoeff(std::array<std::size_t, FactorialInputs> const &factorials,
                 std::array<std::size_t, Bins> const &coeffs) -> std::size_t {
    std::size_t res = factorials[Order];

    std::size_t k = 0;
    while (k < Bins) {
        res /= factorials[coeffs[k]];
        k++;
    }
    return res;
}

} // namespace adhoc3

#endif // ADHOC3_PARTITION_MULTINOMIAL_COEFFICIENT_HPP
