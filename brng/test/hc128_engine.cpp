#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "seed_seq_inserter.hpp"

#include <algorithm>
#include <hc128_engine.hpp>

#include "external/wuhjresearch/hc128_ref.h"
#include "read_csv.hpp"

#include <algorithm>
#include <chrono>
#include <ranges>

namespace {

template <class SeedSeq> auto init_finder(SeedSeq &seq) {
    std::size_t idx = 0;
    adhoc::circular_buffer<std::uint32_t, 512> P;
    adhoc::circular_buffer<std::uint32_t, 512> Q;

    seq.generate(P.data().begin(), P.data().end());
    seq.generate(Q.data().begin(), Q.data().end());

    std::vector<std::uint32_t> init_array(8);

    std::array<std::uint32_t, 1024 + 256> W{};

    for (std::size_t i = 0; i < 1024; i++) {
        constexpr auto mask9 = adhoc::mask<std::size_t>(9);
        constexpr auto mask10 = adhoc::mask<std::size_t>(10);

        --idx;
        idx &= mask10;

        bool const use_first = idx < 512U;
        auto &PQ = use_first ? P : Q;
        --PQ;
        auto const &QP_data = use_first ? Q.data() : P.data();

        auto const u = PQ.template at<512 - 12>();
        auto const a = static_cast<std::uint8_t>(u);
        auto const c = static_cast<std::uint8_t>(u >> 16U);

        PQ.at() ^= (QP_data[a] + QP_data[256 + c]);

        int const rot1 = use_first ? 10 : -10;
        int const rot2 = use_first ? 23 : -23;
        int const rot3 = use_first ? 8 : -8;

        PQ.at() -= (std::rotr((PQ.template at<512 - 3>()), rot1) ^
                    std::rotr((PQ.template at<1>()), rot2)) +
                   std::rotr((PQ.template at<512 - 10>()), rot3);
    }

    for (std::size_t i = 512; i-- > 0;) {
        --Q;
        W[i + 256 + 512] = Q.at();
    }

    for (std::size_t i = 512; i-- > 0;) {
        --P;
        W[i + 256] = P.at();
    }

    for (std::size_t i = 1024 + 256; i-- > 16;) {

        auto const tmp1 = std::rotr(W[i - 15], 7) ^ std::rotr(W[i - 15], 18) ^
                          (W[i - 15] >> 3U);

        auto const tmp2 = std::rotr(W[i - 2], 17) ^ std::rotr(W[i - 2], 19) ^
                          (W[i - 2] >> 10U);

        W[i - 16] = W[i] - (tmp2 + W[i - 7] + tmp1 + i);
        W[i] = 0;
    }

    for (std::size_t i = 0; i < 4; i++) {
        init_array[i + 4] = W[i + 8];
    }

    for (std::size_t i = 0; i < 4; i++) {
        init_array[i] = W[i];
    }

    return init_array;
}

// a class to "wrap" the hc128_engine iterator
template <class UIntType> class hc128 final {
  public:
    hc128(UIntType s) : arg(s) {}

    [[nodiscard]] auto begin() -> adhoc::hc128_engine<UIntType> {
        return adhoc::hc128_engine<UIntType>(arg);
    }

    [[nodiscard]] static auto end() -> std::default_sentinel_t { return {}; }

  private:
    UIntType arg;
};

} // namespace

auto main() -> int {

    // ideally should work like this
    // {
    //     auto taken =
    //         std::views::iota(1) | std::views::take(5); // yields 1 2 3 4 5
    //     for (auto const &v : taken | std::views::reverse) {
    //         std::cout << v << " ";
    //     }
    //     std::cout << std::endl;
    // }

    {
        hc128<std::uint32_t> range(1);

        // Wrap it with take_view: this view *does* have a sentinel
        auto taken = std::views::take(range, 5);

        for (auto const &v : taken) {
            std::cout << v << " ";
        }
        std::cout << std::endl;

        for (auto const &v : taken | std::views::reverse) {
            std::cout << v << " ";
        }

        std::cout << std::endl;
    }

    // {
    //     adhoc::hc128_engine<std::uint32_t> range;

    //     // Wrap it with take_view: this view *does* have a sentinel
    //     auto taken = std::views::take(range, 5);

    //     for (auto const &v : taken) {
    //         std::cout << v << " ";
    //     }

    //     std::cout << "\n";
    // }

    {
        adhoc::seed_seq_inserter seq(
            adhoc::readCSV<std::uint32_t>("./randomgen/hc128_state.txt"));
        adhoc::seed_seq_inserter seq2(init_finder(seq));
        adhoc::hc128_engine<std::uint32_t> rng(seq2);

        std::vector<std::uint32_t> values_from_python =
            adhoc::split_uint64_to_uint32(
                adhoc::readCSV<std::uint64_t>("./randomgen/hc128_vals.txt"));

        for (auto val1 : values_from_python) {
            auto val2 = rng();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        HC128_State state;
        unsigned char key[16], iv[16];
        for (std::size_t i = 0; i < 16; i++) {
            key[i] = 0;
            iv[i] = 0;
        }

        Initialization(&state, key, iv);

        adhoc::hc128_engine<std::uint32_t> rng(0);

        for (std::size_t i = 0; i < 1000000; ++i) {
            OneStep(&state);
            auto val1 = state.keystreamword;
            auto val2 = rng();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::hc128_engine<std::uint32_t> rng(1);
        check_fwd_and_back(rng, 1000000);
        adhoc::hc128_engine<std::uint32_t> rng2(1);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::hc128_engine<std::uint32_t> rng(1);
        check_back_and_fwd(rng, 1000000);
        adhoc::hc128_engine<std::uint32_t> rng2(1);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::hc128_engine<std::uint32_t> rng1(1);
        adhoc::hc128_engine<std::uint64_t> rng2(1);
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::hc128_engine<std::uint64_t> rng(1);
        check_fwd_and_back(rng, 1000000);
        adhoc::hc128_engine<std::uint64_t> rng2(1);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::hc128_engine<std::uint64_t> rng(1);
        check_back_and_fwd(rng, 1000000);
        adhoc::hc128_engine<std::uint64_t> rng2(1);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::hc128_engine<std::uint64_t> rng(1);
        check_discard(rng, 100);
    }

    static_assert(
        std::bidirectional_iterator<adhoc::hc128_engine<std::uint32_t>>);
    {
        adhoc::hc128_engine<std::uint32_t> rng1(1);
        std::vector<std::uint32_t> v(10);
        std::ranges::generate(v, rng1);

        adhoc::hc128_engine<std::uint32_t> rng2(1);
        std::vector<std::uint32_t> v2(10);

        for (auto &val : v2) {
            val = *(++rng2);
        }

        for (std::size_t i = 0; i < v.size(); ++i) {
            EXPECT_EQUAL(v[i], v2[i]);
        }
    }

    {
        adhoc::hc128_engine<std::uint32_t> rng(1);
        std::vector<std::uint32_t> v(10);
        for (auto &val : v) {
            val = *(++rng);
        }
        std::vector<std::uint32_t> v2(10);

        for (auto &val : v2 | std::views::reverse) {
            val = *(rng);
            --rng;
        }

        for (std::size_t i = 0; i < v.size(); ++i) {
            EXPECT_EQUAL(v[i], v2[i]);
        }
    }

    std::size_t sims = 0;
    if (auto env_p = std::getenv("TIMING_SIMS")) {
        sims = std::atoi(env_p);
    }

    if (sims) {
        std::uint64_t res1 = 0;
        adhoc::hc128_engine<std::uint32_t> rng;
        {
            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims; i++) {
                res1 += rng();
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();
            std::cout << "fwd" << std::endl;
            std::cout << time << std::endl;
        }

        std::uint64_t res1b = 0;
        {
            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims; i++) {
                res1b += rng.operator()<false>();
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();
            std::cout << "back" << std::endl;
            std::cout << time << std::endl;
        }
        EXPECT_EQUAL(res1, res1b);
    }

    TEST_END;
}
