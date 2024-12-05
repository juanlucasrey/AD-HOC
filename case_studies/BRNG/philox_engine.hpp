#ifndef CASE_STUDIES_BRNG_PHILOX_ENGINE
#define CASE_STUDIES_BRNG_PHILOX_ENGINE

#include <array>
#include <cassert>
#include <cstdint>
#include <iterator>
#include <ranges>

namespace adhoc {

namespace detail {
constexpr auto next_stdint(unsigned w) -> unsigned {
    if (w <= 8) {
        return 8;
    }
    if (w <= 16) {
        return 16;
    }
    if (w <= 32) {
        return 32;
    }
    if (w <= 64) {
        return 64;
    }
    if (w <= 128) {
        return 128; // DANGER - assuming __uint128_t
    }
    return 0;
}

template <unsigned W> struct ui {};
template <> struct ui<8> {
    using least = uint_least8_t;
    using fast = uint_fast8_t;
};
template <> struct ui<16> {
    using least = uint_least16_t;
    using fast = uint_fast16_t;
};
template <> struct ui<32> {
    using least = uint_least32_t;
    using fast = uint_fast32_t;
};
template <> struct ui<64> {
    using least = uint_least64_t;
    using fast = uint_fast64_t;
};
template <> struct ui<128> {
    using least = __uint128_t;
    using fast = __uint128_t;
}; // DANGER - __uint128_t

template <unsigned w> using uint_least = typename ui<next_stdint(w)>::least;
template <unsigned w> using uint_fast = typename ui<next_stdint(w)>::fast;

template <unsigned W, std::unsigned_integral U>
static auto mulhilo(U a, U b) -> std::pair<U, U> {
    using uwide = uint_fast<2 * W>;
    const size_t xwidth = std::numeric_limits<uwide>::digits;
    const uwide ab = uwide(a) * uwide(b);
    return {U(ab >> W), U((ab << (xwidth - W)) >> (xwidth - W))};
}

} // namespace detail

template <typename UIntType, std::size_t w, std::size_t n, std::size_t r,
          UIntType... consts>
class philox_engine {
    static_assert(sizeof...(consts) == n);
    static_assert(n == 2 || n == 4 || n == 8 || n == 16);
    static_assert(0 < r);
    static_assert(0 < w && w <= std::numeric_limits<UIntType>::digits);

  public:
    using result_type = UIntType;

    static constexpr std::size_t word_size = w;
    static constexpr std::size_t word_count = n;
    static constexpr std::size_t round_count = r;
    static constexpr result_type default_seed = 20111115U;

    static constexpr auto min() -> result_type {
        return static_cast<result_type>(0U);
    }
    static constexpr auto max() -> result_type {
        return (result_type(~(result_type(0))) >>
                (std::numeric_limits<result_type>::digits - w));

        // equivalent to the following but shift does not overflow
        // return (static_cast<result_type>(1U) << w) - 1U;
    }

    philox_engine() : philox_engine(default_seed) {}
    explicit philox_engine(result_type value) {
        constexpr result_type mask = max();
        this->K[0] = value & mask;
    }

    auto operator()() -> result_type {
        this->j++;
        if (this->j == n) {
            this->generate();

            constexpr auto in_mask = max();
            std::size_t i = 0;
            do {
                this->X[i] = (this->X[i] + 1) & in_mask;
                ++i;
            } while (i < n && !this->X[i - 1]);

            this->j = 0;
        }
        return Y[this->j];
    }
    void discard(unsigned long long z) {
        for (unsigned long long i = 0; i < z; ++i) {
            this->operator()();
        }
    }

  private:
    void generate() {
        constexpr auto in_mask = max();
        constexpr std::array<UIntType, n> consts_arr{consts...};
        if constexpr (n == 2) {
            result_type R0 = this->X[0];
            result_type L0 = this->X[1];
            result_type K0 = this->K[0];
            for (size_t i = 0; i < round_count; ++i) {
                auto [hi, lo] = detail::mulhilo<word_size>(R0, consts_arr[0]);
                R0 = hi ^ L0 ^ K0;
                L0 = lo;
                K0 = (K0 + consts_arr[1]) & in_mask;
            }
            Y[0] = R0;
            Y[1] = L0;
        } else if constexpr (n == 4) {
            result_type R0 = this->X[0];
            result_type L0 = this->X[1];
            result_type R1 = this->X[2];
            result_type L1 = this->X[3];
            result_type K0 = this->K[0];
            result_type K1 = this->K[1];
            for (size_t i = 0; i < r; ++i) {
                auto [hi0, lo0] = detail::mulhilo<word_size>(R0, consts_arr[0]);
                auto [hi1, lo1] = detail::mulhilo<word_size>(R1, consts_arr[2]);
                R0 = hi1 ^ L0 ^ K0;
                L0 = lo1;
                R1 = hi0 ^ L1 ^ K1;
                L1 = lo0;
                K0 = (K0 + consts_arr[1]) & in_mask;
                K1 = (K1 + consts_arr[3]) & in_mask;
            }
            Y[0] = R0;
            Y[1] = L0;
            Y[2] = R1;
            Y[3] = L1;
        }
    }

    std::array<UIntType, n> X{0};
    std::array<UIntType, n / 2> K{0};
    std::array<UIntType, n> Y{0};
    std::size_t j{n - 1};
};

using philox4x32 = philox_engine<std::uint_fast32_t, 32, 4, 10, 0xD2511F53,
                                 0x9E3779B9, 0xCD9E8D57, 0xBB67AE85>;

using philox4x64 =
    philox_engine<std::uint_fast64_t, 64, 4, 10, 0xD2E7470EE14C6C93,
                  0x9E3779B97F4A7C15, 0xCA5A826395121157, 0xBB67AE8584CAA73B>;

} // namespace adhoc

#endif // CASE_STUDIES_BRNG_PHILOX_ENGINE
