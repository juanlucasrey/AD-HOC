#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "seed_seq_inserter.hpp"

#include "PractRand.h"
#include "PractRand/RNGs/all.h"

#include <hc256_engine.hpp>

#include "external/wuhjresearch/hc256_ref.h"

#include <algorithm>
#include <chrono>
#include <numeric>
#include <ranges>

auto main() -> int {
    {
        HC256_State state;
        unsigned char key[32], iv[32];
        for (std::size_t i = 0; i < 32; i++) {
            key[i] = 0;
            iv[i] = 0;
        }

        Initialization(&state, key, iv);

        adhoc::hc256_engine<std::uint32_t> rng(0);

        for (std::size_t i = 0; i < 1000000; ++i) {
            OneStep(&state);
            auto val1 = state.keystreamword;
            auto val2 = rng();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        HC256_State state;
        unsigned char key[32], iv[32];
        std::vector<std::uint8_t> init(64);
        std::iota(init.begin(), init.end(), 4);

        std::size_t idx = 0;
        for (std::size_t i = 0; i < 32; i++) {
            key[i] = static_cast<std::uint32_t>(init[idx++]);
        }

        for (std::size_t i = 0; i < 32; i++) {
            iv[i] = static_cast<std::uint32_t>(init[idx++]);
        }

        Initialization(&state, key, iv);

        adhoc::seed_seq_inserter seq(init);

        adhoc::hc256_engine<std::uint32_t> rng(seq);

        for (std::size_t i = 0; i < 1000000; ++i) {
            OneStep(&state);
            auto val1 = state.keystreamword;
            auto val2 = rng();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        PractRand::RNGs::LightWeight::hc256 rng1(1);
        adhoc::hc256_engine<std::uint32_t> rng2(1);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1.raw32();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::hc256_engine<std::uint32_t> rng1(1);
        adhoc::hc256_engine<std::uint64_t> rng2(1);

        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::hc256_engine<std::uint32_t> rng(1);
        for (std::size_t i = 0; i < 4096; ++i) {
            rng.template operator()<false>();
        }

        check_fwd_and_back(rng, 1000000);
        adhoc::hc256_engine<std::uint32_t> rng2(1);
        for (std::size_t i = 0; i < 4096; ++i) {
            rng2.template operator()<false>();
        }
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::hc256_engine<std::uint32_t> rng(1);
        check_fwd_and_back(rng, 1000000);
        adhoc::hc256_engine<std::uint32_t> rng2(1);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::hc256_engine<std::uint32_t> rng(1);
        check_back_and_fwd(rng, 1000000);
        adhoc::hc256_engine<std::uint32_t> rng2(1);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::hc256_engine<std::uint64_t> rng(1);
        check_fwd_and_back(rng, 1000000);
        adhoc::hc256_engine<std::uint64_t> rng2(1);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::hc256_engine<std::uint64_t> rng(1);
        check_back_and_fwd(rng, 1000000);
        adhoc::hc256_engine<std::uint64_t> rng2(1);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::hc256_engine<std::uint64_t> rng(1);
        check_discard(rng, 100);
    }

    static_assert(
        std::bidirectional_iterator<adhoc::hc256_engine<std::uint32_t>>);
    {
        adhoc::hc256_engine<std::uint32_t> rng1(1);
        std::vector<std::uint32_t> v(10);
        std::ranges::generate(v, rng1);

        adhoc::hc256_engine<std::uint32_t> rng2(1);
        std::vector<std::uint32_t> v2(10);

        for (auto &val : v2) {
            val = *rng2;
            ++rng2;
        }

        for (std::size_t i = 0; i < v.size(); ++i) {
            EXPECT_EQUAL(v[i], v2[i]);
        }
    }

    {
        adhoc::hc256_engine<std::uint32_t> rng(1);
        std::vector<std::uint32_t> v(10);
        for (auto &val : v) {
            val = *rng;
            ++rng;
        }
        std::vector<std::uint32_t> v2(10);

        for (auto &val : v2 | std::views::reverse) {
            --rng;
            val = *(rng);
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
        adhoc::hc256_engine<std::uint32_t> rng;
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
