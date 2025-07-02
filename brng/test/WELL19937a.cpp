#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

extern "C" {
#include "third_party/umontreal/WELL19937a.h"
}

#include "distribution/uniform_distribution.hpp"
#include "seed_seq_inserter.hpp"
#include <WELL19937a.hpp>

#include <chrono>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>

auto main() -> int {
    {
        std::vector<std::uint32_t> seqval;
        seqval.resize(624);
        std::iota(seqval.begin(), seqval.end(), 1000);
        adhoc::seed_seq_inserter seq(seqval);

        adhoc::canonical<adhoc::WELL19937a<std::uint_fast32_t>> rng(seq);

        std::array<unsigned int, 624> seed;
        std::iota(seed.begin(), seed.end(), 1000);
        InitWELLRNG19937a(seed.data());

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = WELLRNG19937a();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        std::seed_seq seq{1, 2, 3, 4, 5};
        adhoc::WELL19937a<std::uint32_t> rng(seq);
        check_fwd_and_back(rng, 1000000);
        adhoc::WELL19937a<std::uint32_t> rng2(seq);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        std::seed_seq seq{1, 2, 3, 4, 5};
        adhoc::WELL19937a<std::uint32_t> rng(seq);
        check_back_and_fwd(rng, 1000000);
        adhoc::WELL19937a<std::uint32_t> rng2(seq);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        std::seed_seq seq1{1, 2, 3, 4, 5};
        adhoc::WELL19937a<std::uint32_t> rng1(seq1);
        std::seed_seq seq2{1, 2, 3, 4, 5};
        adhoc::WELL19937a<std::uint64_t> rng2(seq2);
        compare_rng(rng1, rng2, 100000);
        compare_rng_limits(rng1, rng2);
    }

    {
        std::seed_seq seq{1, 2, 3, 4, 5};
        adhoc::WELL19937a<std::uint64_t> rng(seq);
        check_fwd_and_back(rng, 1000000);
        adhoc::WELL19937a<std::uint64_t> rng2(seq);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        std::seed_seq seq{1, 2, 3, 4, 5};
        adhoc::WELL19937a<std::uint64_t> rng(seq);
        check_back_and_fwd(rng, 1000000);
        adhoc::WELL19937a<std::uint64_t> rng2(seq);
        EXPECT_EQUAL(rng, rng2);
    }

    std::size_t sims = 0;
    if (auto env_p = std::getenv("TIMING_SIMS")) {
        sims = std::atoi(env_p);
    }

    if (sims) {
        double res1 = 0;
        std::vector<std::uint32_t> seqval;
        seqval.resize(624);
        std::iota(seqval.begin(), seqval.end(), 1000);
        adhoc::seed_seq_inserter seq(seqval);
        adhoc::canonical<adhoc::WELL19937a<std::uint_fast32_t>> rng(seq);
        {
            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims; i++) {
                res1 += rng();
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();
            std::cout << "new" << std::endl;
            std::cout << time << std::endl;
        }

        double res2 = 0;
        std::array<unsigned int, 624> seed;
        std::iota(seed.begin(), seed.end(), 1000);
        InitWELLRNG19937a(seed.data());
        {
            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims; i++) {
                res2 += WELLRNG19937a();
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();
            std::cout << "old" << std::endl;
            std::cout << time << std::endl;
        }
        EXPECT_EQUAL(res1, res2);

        double res1b = 0;
        {
            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims; i++) {
                res1b += rng.operator()<false>();
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();
            std::cout << "new back" << std::endl;
            std::cout << time << std::endl;
        }
        EXPECT_NEAR_REL(res1, res1b, 1e-5);
    }

    TEST_END;
}
