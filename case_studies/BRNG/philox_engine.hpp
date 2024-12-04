#ifndef CASE_STUDIES_BRNG_PHILOX_ENGINE
#define CASE_STUDIES_BRNG_PHILOX_ENGINE

#include <array>
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

template <typename UIntType, typename Tuple, size_t... Is>
constexpr auto get_even_array_from_tuple(Tuple t, std::index_sequence<Is...>) {
    return std::array<UIntType, std::index_sequence<Is...>::size()>{
        std::get<Is * 2>(t)...};
}

template <typename UIntType, typename Tuple, size_t... Is>
constexpr auto get_odd_array_from_tuple(Tuple t, std::index_sequence<Is...>) {
    return std::array<UIntType, std::index_sequence<Is...>::size()>{
        std::get<Is * 2 + 1>(t)...};
}

} // namespace detail

template <typename UIntType, std::size_t w, std::size_t n, std::size_t r,
          UIntType... consts>
class philox_engine {
    static_assert(sizeof...(consts) == n);
    static_assert(n == 2 || n == 4 || n == 8 || n == 16);
    static_assert(0 < r);
    static_assert(0 < w && w <= std::numeric_limits<UIntType>::digits);

    static constexpr std::size_t array_size = n / 2;
    // static constexpr std::size_t key_count = n / 2; // keys_count

  public:
    using result_type = UIntType;

    static constexpr std::size_t word_size = w;
    static constexpr std::size_t word_count = n;
    static constexpr std::size_t round_count = r;
    static constexpr std::array<result_type, array_size> multipliers =
        detail::get_even_array_from_tuple<UIntType>(
            std::make_tuple(consts...), std::make_index_sequence<array_size>{});
    static constexpr std::array<result_type, array_size> round_consts =
        detail::get_odd_array_from_tuple<UIntType>(
            std::make_tuple(consts...), std::make_index_sequence<array_size>{});
    static constexpr result_type default_seed = 20111115U;

    static constexpr result_type min() { return static_cast<result_type>(0U); }
    static constexpr result_type max() {
        return (result_type(~(result_type(0))) >>
                (std::numeric_limits<result_type>::digits - w));

        // equivalent to the following but shift does not overflow
        // return (static_cast<result_type>(1U) << w) - 1U;
    }

    // constructors and seeding functions
    philox_engine() : philox_engine(default_seed) {}
    explicit philox_engine(result_type value) {
        std::size_t i = 0;
        for (i = 0; i < word_count; i++) {
            state[i] = 0;
        }

        constexpr result_type mask = max();

        state[i] = value & mask;
        // auto start = value.begin();
        // auto end = value.end();
        // for (; i < state_size; i++) { // keys are set as seed
        //     state[i] = (start == end) ? 0 : (*start++) & mask;
        // }
        results[0] = 0;
    }
    // template <class Sseq> explicit philox_engine(Sseq &q);
    // void seed(result_type value = default_seed);
    // template <class Sseq> void seed(Sseq &q);
    // void set_counter(const std::array<result_type, n> &counter);
    // generating functions

    static constexpr std::size_t period_counter_count = 1;
    using counter_type = detail::uint_fast<period_counter_count *
                                           word_size>; // WARNING: doesn't scale
                                                       // for word_count > 4

  private:
    static constexpr size_t state_size = 3 * n / 2; // counters + keys
    std::array<result_type, state_size> state{
        0}; // state: [counter_0,..., counter_n, key_0, ..., key_n/2-1];

    std::array<result_type, word_count> results{0};

    counter_type get_counter_internal() const { // need to check
        std::uint64_t ret = 0;
        for (size_t i = 0; i < period_counter_count; ++i) {
            ret |= std::uint64_t(state[i]) << (word_size * i);
        }
        return ret;
    }

    void set_counter_internal(std::array<result_type, state_size> &s,
                              counter_type newctr) { // need to check
        constexpr auto in_mask = max();
        static_assert(word_size * period_counter_count <=
                      std::numeric_limits<counter_type>::digits);
        for (size_t i = 0; i < period_counter_count; ++i)
            s[i] = (newctr >> (word_size * i)) & in_mask;
    }

    void increase_counter_internal() { // need to check
        constexpr auto in_mask = max();
        state[0] = (state[0] + 1) & in_mask;
        for (size_t i = 1; i < period_counter_count; ++i) {
            if (state[i - 1]) {
                [[likely]] return;
            }
            state[i] = (state[i] + 1) & in_mask;
        }
    }

    template <std::ranges::input_range InRange,
              std::weakly_incrementable Output>
    requires std::ranges::sized_range<InRange> &&
             std::integral<
                 std::iter_value_t<std::ranges::range_value_t<InRange>>> &&
             std::integral<std::iter_value_t<Output>> &&
             std::indirectly_writable<Output, std::iter_value_t<Output>>
    Output generate(InRange &&inrange, Output output) {
        constexpr auto in_mask = max();
        for (auto initer : inrange) {
            if constexpr (n == 2) {
                result_type R0 = (*initer++) & in_mask;
                result_type L0 = (*initer++) & in_mask;
                result_type K0 = (*initer++) & in_mask;
                for (size_t i = 0; i < round_count; ++i) {
                    auto [hi, lo] =
                        detail::mulhilo<word_size>(R0, multipliers[0]);
                    R0 = hi ^ K0 ^ L0;
                    L0 = lo;
                    K0 = (K0 + round_consts[0]) & in_mask;
                }
                *output++ = R0;
                *output++ = L0;
            } else if constexpr (n == 4) {
                result_type R0 = (*initer++) & in_mask;
                result_type L0 = (*initer++) & in_mask;
                result_type R1 = (*initer++) & in_mask;
                result_type L1 = (*initer++) & in_mask;
                result_type K0 = (*initer++) & in_mask;
                result_type K1 = (*initer++) & in_mask;
                for (size_t i = 0; i < r; ++i) {
                    auto [hi0, lo0] =
                        detail::mulhilo<word_size>(R0, multipliers[0]);
                    auto [hi1, lo1] =
                        detail::mulhilo<word_size>(R1, multipliers[1]);
                    R0 = hi1 ^ L0 ^ K0;
                    L0 = lo1;
                    R1 = hi0 ^ L1 ^ K1;
                    L1 = lo0;
                    K0 = (K0 + round_consts[0]) & in_mask;
                    K1 = (K1 + round_consts[1]) & in_mask;
                }
                *output++ = R0;
                *output++ = L0;
                *output++ = R1;
                *output++ = L1;
            }
            // No more cases.  See the static_assert(n==2 || n==4) at the top of
            // the class
        }
        return output;
    }

  public:
    result_type operator()() {
        result_type ret;
        auto out = &ret;
        std::size_t len = 1;
        auto ri = results[0];
        if (ri && len) {
            while (ri < word_count && len) {
                *out++ = results[ri++];
                --len;
            }
            if (ri == word_count)
                ri = 0;
        }

        auto nprf = len / word_count;

        auto c0 = get_counter_internal();
        std::array<result_type, state_size> state_tmp;
        out = generate(std::ranges::views::iota(c0, c0 + nprf) |
                           std::ranges::views::transform([&](auto ctr) {
                               state_tmp = state;
                               set_counter_internal(state_tmp, ctr);
                               return std::ranges::begin(state_tmp);
                           }),
                       out);
        len -= nprf * word_count;
        set_counter_internal(state, c0 + nprf);

        // Restock the results array
        if (ri == 0 && len) {
            generate(std::ranges::single_view(state.data()), results.data());
            increase_counter_internal();
        }

        // Finish off any stragglers.
        while (len--)
            *out++ = results[ri++];
        results[0] = ri;
        // return *out;
        return ret;
    }
    void discard(unsigned long long z) {
        for (unsigned long long i = 0; i < z; ++i) {
            this->operator()();
        }
    }
};

using philox4x32 = philox_engine<std::uint_fast32_t, 32, 4, 10, 0xD2511F53,
                                 0x9E3779B9, 0xCD9E8D57, 0xBB67AE85>;

using philox4x64 =
    philox_engine<std::uint_fast64_t, 64, 4, 10, 0xD2E7470EE14C6C93,
                  0x9E3779B97F4A7C15, 0xCA5A826395121157, 0xBB67AE8584CAA73B>;

} // namespace adhoc

#endif // CASE_STUDIES_BRNG_PHILOX_ENGINE
