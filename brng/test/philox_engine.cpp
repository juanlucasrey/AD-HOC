#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include <philox_engine.hpp>

#include "read_csv.hpp"
#include "seed_seq_inserter.hpp"

#include <chrono>
#include <iostream>
#include <random>

#include "third_party/Random123/philox.h"

int main() {

    {
        using namespace r123;
        Philox2x32 rng;
        Philox2x32::ctr_type c = {{0, 0}};
        Philox2x32::key_type k = {{5}};

        adhoc::philox2x32 rng1(5);
        for (unsigned int j = 0; j < 10; ++j) {
            auto r = rng(c, k);
            ++c[0];
            auto val1 = rng1();
            auto val2 = rng1();
            EXPECT_EQUAL(r[0], val1);
            EXPECT_EQUAL(r[1], val2);
        }
    }

    {
        using namespace r123;
        Philox4x32 rng;
        Philox4x32::ctr_type c = {{0, 0, 0, 0}};
        Philox4x32::key_type k = {{5, 1}};

        adhoc::philox4x32 rng1(5, 1);
        for (unsigned int j = 0; j < 10; ++j) {
            auto r = rng(c, k);
            ++c[0];
            auto val1 = rng1();
            auto val2 = rng1();
            auto val3 = rng1();
            auto val4 = rng1();
            EXPECT_EQUAL(r[0], val1);
            EXPECT_EQUAL(r[1], val2);
            EXPECT_EQUAL(r[2], val3);
            EXPECT_EQUAL(r[3], val4);
        }
    }

    {
        using namespace r123;
        Philox2x64 rng;
        Philox2x64::ctr_type c = {{0, 0}};
        Philox2x64::key_type k = {{5}};

        adhoc::philox2x64 rng1(5);
        for (unsigned int j = 0; j < 10; ++j) {
            auto r = rng(c, k);
            ++c[0];
            auto val1 = rng1();
            auto val2 = rng1();
            EXPECT_EQUAL(r[0], val1);
            EXPECT_EQUAL(r[1], val2);
        }
    }

    {
        using namespace r123;
        Philox4x64 rng;
        Philox4x64::ctr_type c = {{0, 0, 0, 0}};
        Philox4x64::key_type k = {{5, 1}};

        adhoc::philox4x64 rng1(5, 1);
        for (unsigned int j = 0; j < 10; ++j) {
            auto r = rng(c, k);
            ++c[0];
            auto val1 = rng1();
            auto val2 = rng1();
            auto val3 = rng1();
            auto val4 = rng1();
            EXPECT_EQUAL(r[0], val1);
            EXPECT_EQUAL(r[1], val2);
            EXPECT_EQUAL(r[2], val3);
            EXPECT_EQUAL(r[3], val4);
        }
    }

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

    static_assert(std::bidirectional_iterator<adhoc::philox2x32>);
    static_assert(std::bidirectional_iterator<adhoc::philox4x32>);
    static_assert(std::bidirectional_iterator<adhoc::philox2x64>);
    static_assert(std::bidirectional_iterator<adhoc::philox4x64>);

    std::size_t sims = 0;
    if (auto env_p = std::getenv("TIMING_SIMS")) {
        sims = std::atoi(env_p);
    }

    if (sims) {
        std::uint64_t res_ref = 0;
        using namespace r123;
        Philox4x64 rng_ref;
        Philox4x64::ctr_type c = {{0, 0, 0, 0}};
        Philox4x64::key_type k = {{5, 4}};
        std::size_t sims_ov_4 = sims / 4;

        {
            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims_ov_4; i++) {
                auto const r = rng_ref(c, k);
                ++c[0];
                for (std::size_t j = 0; j < 4; j++) {
                    res_ref += r[j];
                }
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();
            std::cout << "reference" << std::endl;
            std::cout << time << std::endl;
        }

        std::uint64_t res1 = 0;
        adhoc::philox4x64 rng(5, 4);
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
        EXPECT_EQUAL(res1, res_ref);
    }

    TEST_END;
}
