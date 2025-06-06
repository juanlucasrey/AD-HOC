#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

extern "C" {
#include "external/WELL1024a.h"
}

#include <WELL1024a.hpp>
#include <distribution/uniform_distribution.hpp>

#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>

auto main() -> int {
    {
        std::seed_seq seq{1, 2, 3, 4, 5};
        adhoc::canonical<adhoc::WELL1024a<std::uint_fast32_t>> rng(seq);

        std::array<unsigned int, 32> seed = {
            3414683984, 1693102273, 232375452,  2717977097, 3956103128,
            2219488797, 1777128801, 3230442083, 4044208335, 98442330,
            3372076691, 424038268,  3268355489, 3886808778, 2652920978,
            3124645489, 3117374989, 3282134275, 4187435356, 1688570390,
            1812527118, 4140562627, 749212602,  2313212922, 1745485075,
            962989454,  3962071522, 3754576413, 472654705,  1301407037,
            757983457,  897108390};
        InitWELLRNG1024a(seed.data());

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = WELLRNG1024a();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        std::seed_seq seq1{1, 2, 3, 4, 5};
        adhoc::WELL1024a<std::uint32_t> rng1(seq1);
        check_fwd_and_back(rng1, 1000000);
        std::seed_seq seq2{1, 2, 3, 4, 5};
        adhoc::WELL1024a<std::uint32_t> rng2(seq2);
        EXPECT_EQUAL(rng1, rng2);
    }

    {
        std::seed_seq seq1{1, 2, 3, 4, 5};
        adhoc::WELL1024a<std::uint32_t> rng1(seq1);
        check_back_and_fwd(rng1, 1000000);
        std::seed_seq seq2{1, 2, 3, 4, 5};
        adhoc::WELL1024a<std::uint32_t> rng2(seq2);
        EXPECT_EQUAL(rng1, rng2);
    }

    {
        std::seed_seq seq1{1, 2, 3, 4, 5};
        adhoc::WELL1024a<std::uint32_t> rng1(seq1);
        std::seed_seq seq2{1, 2, 3, 4, 5};
        adhoc::WELL1024a<std::uint64_t> rng2(seq2);
        compare_rng(rng1, rng2, 100000);
        compare_rng_limits(rng1, rng2);
    }

    {
        std::seed_seq seq1{1, 2, 3, 4, 5};
        adhoc::WELL1024a<std::uint64_t> rng1(seq1);
        check_fwd_and_back(rng1, 1000000);
        std::seed_seq seq2{1, 2, 3, 4, 5};
        adhoc::WELL1024a<std::uint64_t> rng2(seq2);
        EXPECT_EQUAL(rng1, rng2);
    }

    {
        std::seed_seq seq1{1, 2, 3, 4, 5};
        adhoc::WELL1024a<std::uint64_t> rng1(seq1);
        check_back_and_fwd(rng1, 1000000);
        std::seed_seq seq2{1, 2, 3, 4, 5};
        adhoc::WELL1024a<std::uint64_t> rng2(seq2);
        EXPECT_EQUAL(rng1, rng2);
    }

    std::size_t sims = 0;
    if (auto env_p = std::getenv("TIMING_SIMS")) {
        sims = std::atoi(env_p);
    }

    if (sims) {
        double res1 = 0;
        std::seed_seq seq{1, 2, 3, 4, 5};
        adhoc::canonical<adhoc::WELL1024a<std::uint32_t>> rng(seq);
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
        std::array<unsigned int, 32> seed = {
            3414683984, 1693102273, 232375452,  2717977097, 3956103128,
            2219488797, 1777128801, 3230442083, 4044208335, 98442330,
            3372076691, 424038268,  3268355489, 3886808778, 2652920978,
            3124645489, 3117374989, 3282134275, 4187435356, 1688570390,
            1812527118, 4140562627, 749212602,  2313212922, 1745485075,
            962989454,  3962071522, 3754576413, 472654705,  1301407037,
            757983457,  897108390};
        InitWELLRNG1024a(seed.data());
        {
            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims; i++) {
                res2 += WELLRNG1024a();
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();
            std::cout << "old" << std::endl;
            std::cout << time << std::endl;
        }
        EXPECT_EQUAL(res1, res2);

        double res3 = 0;
        {
            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims; i++) {
                res3 += rng.operator()<false>();
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();
            std::cout << "new back" << std::endl;
            std::cout << time << std::endl;
        }

        EXPECT_NEAR_REL(res1, res3, 1e-5);
    }

    TEST_END;
}
