#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include <philox_engine.hpp>

#include "read_csv.hpp"
#include "seed_seq_inserter.hpp"

#include <chrono>
#include <iostream>
#include <random>

int main() {
    // std check
    {
        adhoc::philox4x32 rng1; // overload (1)
        adhoc::philox4x64 rng2; // overload (1)
        rng1.discard(10000 - 1);
        rng2.discard(10000 - 1);
        EXPECT_EQUAL(rng1(), 1955073260U);
        EXPECT_EQUAL(rng2(), 3409172418970261260U);
    }

    // check full cycle equality
    {
        adhoc::philox4x64 rng1;
        adhoc::philox4x64 rng2;
        rng2.set_counter(std::array<std::uint_fast64_t, 4UL>{
            18446744073709551615U, 18446744073709551615U, 18446744073709551615U,
            18446744073709551615U});
        for (unsigned int i = 0; i < 4; ++i) {
            rng2();
        }
        EXPECT_EQUAL(rng1, rng2);
    }

    {
        adhoc::seed_seq_inserter seq(
            adhoc::readCSV<std::uint64_t>("./randomgen/philox_state.txt"));
        adhoc::philox4x64 rng(seq);

        auto const values_from_python =
            adhoc::readCSV<std::uint64_t>("./randomgen/philox_vals.txt");

        // c++ starts 4 sims before c++ impl so we skip those
        rng.discard(4);
        for (auto val1 : values_from_python) {
            auto val2 = rng();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        std::seed_seq s1{1, 2, 3};

        adhoc::philox4x32 rng1(s1);
        check_fwd_and_back(rng1, 1000000);

        std::seed_seq s2{1, 2, 3};
        adhoc::philox4x32 rng2(s2);
        EXPECT_EQUAL(rng1, rng2);
    }

    {
        adhoc::philox4x32 rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::philox4x32 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::philox4x32 rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::philox4x32 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::philox4x64 rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::philox4x64 rng2;
        EXPECT_EQUAL(rng, rng2);
    }
    {
        adhoc::philox4x64 rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::philox4x64 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    // check 32->64
    {
        adhoc::philox4x32 rng1; // overload (1)
        adhoc::philox_engine<std::uint64_t, 32, 4, 10, 0xCD9E8D57, 0x9E3779B9,
                             0xD2511F53, 0xBB67AE85>
            rng2; // overload (1)
        compare_rng(rng1, rng2, 100000);
        compare_rng_limits(rng1, rng2);
    }

    std::size_t sims = 0;
    if (auto env_p = std::getenv("TIMING_SIMS")) {
        sims = std::atoi(env_p);
    }

    if (sims) {
        std::uint64_t res1 = 0;
        adhoc::philox4x32 rng;
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
            std::cout << "new back" << std::endl;
            std::cout << time << std::endl;
        }
        EXPECT_EQUAL(res1, res1b);
    }

    TEST_END;
}
