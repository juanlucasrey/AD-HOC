#ifndef CASE_STUDIES_BRNG_SUBTRACT_WITH_CARRY_ENGINE
#define CASE_STUDIES_BRNG_SUBTRACT_WITH_CARRY_ENGINE

#include "linear_congruential_engine.hpp"

#include <algorithm>
#include <array>
#include <cmath>

#include <iostream>

namespace adhoc {

template <class UIntType, std::size_t w, std::size_t s, std::size_t r>
class subtract_with_carry_engine final {
  public:
    using result_type = UIntType;

    static constexpr std::size_t word_size = w;
    static constexpr std::size_t short_lag = s;
    static constexpr std::size_t long_lag = r;
    static constexpr std::uint_least32_t default_seed = 19780503U;

    static constexpr result_type modulus = static_cast<result_type>(1U) << w;

    explicit subtract_with_carry_engine(result_type value = default_seed) {
        linear_congruential_engine<std::uint_least32_t, 40014U, 0U, 2147483563U>
            e(value == 0U ? default_seed : value);
        constexpr auto n = std::size_t(w / 32) + 1;
        std::array<std::uint_least32_t, n * r> temp;
        std::generate(temp.begin(), temp.end(), e);

        for (std::size_t j = 0; j < r; j++) {
            UIntType val = 0;
            for (std::size_t k = 0; k < n; ++k) {
                val += static_cast<UIntType>(temp[n * j + k]) << 32 * k;
            }
            x[j] = val % modulus;
        }

        carry = (x[long_lag - 1] == 0);
        k = 0;
    }

    result_type operator()() {
        std::size_t short_index =
            (k < short_lag) ? (k + long_lag - short_lag) : (k - short_lag);

        UIntType temp = x[k] + carry;
        if (x[short_index] >= temp) {
            x[k] = x[short_index] - temp;
            carry = 0;
        } else {
            x[k] = modulus - temp + x[short_index];
            carry = 1;
        }

        UIntType result = x[k];
        ++k;
        if (k >= long_lag) {
            k = 0;
        }
        return result;
    }

    void discard(unsigned long long z) {
        for (unsigned long long i = 0; i < z; ++i) {
            this->operator()();
        }
    }

    static constexpr auto min() -> UIntType {
        return static_cast<result_type>(0U);
    }
    static constexpr auto max() -> UIntType {
        return (static_cast<result_type>(1U) << w) - 1;
    }

  private:
    std::array<UIntType, long_lag> x;
    std::size_t k{};
    UIntType carry;
};

using ranlux24_base =
    subtract_with_carry_engine<std::uint_fast32_t, 24, 10, 24>;

using ranlux48_base = subtract_with_carry_engine<std::uint_fast64_t, 48, 5, 12>;

} // namespace adhoc

#endif // CASE_STUDIES_BRNG_MERSENNE_TWISTER_ENGINE
