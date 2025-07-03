#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "read_csv.hpp"
#include "seed_seq_inserter.hpp"
#include <threefry_engine.hpp>
#include <tools/uint128.hpp>

#include <chrono>
#include <iostream>

#include "third_party/Random123/threefry.h"

auto main() -> int {

    {
        using namespace r123;
        Threefry2x32 rng;
        Threefry2x32::ctr_type c = {{0, 0}};
        Threefry2x32::key_type k = {{5, 4}};

        adhoc::threefry2x32 rng1(5, 4);

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
        Threefry4x32 rng;
        Threefry4x32::ctr_type c = {{0, 0, 0, 0}};
        Threefry4x32::key_type k = {{5, 4, 3, 2}};

        adhoc::threefry4x32 rng1(5, 4, 3, 2);

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
        Threefry2x64 rng;
        Threefry2x64::ctr_type c = {{0, 0}};
        Threefry2x64::key_type k = {{5, 4}};

        adhoc::threefry2x64 rng1(5, 4);

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
        Threefry4x64 rng;
        Threefry4x64::ctr_type c = {{0, 0, 0, 0}};
        Threefry4x64::key_type k = {{5, 4, 3, 2}};

        adhoc::threefry4x64 rng1(5, 4, 3, 2);

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
        adhoc::seed_seq_inserter seq(
            adhoc::readCSV<std::uint64_t>("./randomgen/threefry_state.txt"));
        adhoc::threefry4x64 rng(seq);

        auto const values_from_python =
            adhoc::readCSV<std::uint64_t>("./randomgen/threefry_vals.txt");

        // c++ starts 4 sims before c++ impl so we skip those
        rng.discard(4);
        for (auto val1 : values_from_python) {
            auto val2 = rng();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::threefry2x32 rng1(5, 4);
        check_fwd_and_back(rng1, 1000000);
        adhoc::threefry2x32 rng2(5, 4);
        EXPECT_EQUAL(rng1, rng2);
    }

    {
        adhoc::threefry2x32 rng1(5, 4);
        check_back_and_fwd(rng1, 1000000);
        adhoc::threefry2x32 rng2(5, 4);
        EXPECT_EQUAL(rng1, rng2);
    }

    {
        adhoc::threefry4x32 rng1(5, 4, 3, 2);
        check_fwd_and_back(rng1, 1000000);
        adhoc::threefry4x32 rng2(5, 4, 3, 2);
        EXPECT_EQUAL(rng1, rng2);
    }

    {
        adhoc::threefry4x32 rng1(5, 4, 3, 2);
        check_back_and_fwd(rng1, 1000000);
        adhoc::threefry4x32 rng2(5, 4, 3, 2);
        EXPECT_EQUAL(rng1, rng2);
    }

    // check 32->64
    {
        adhoc::threefry_engine<std::uint32_t, 32, 2, 20, 0x1BD11BDA, 13, 15, 26,
                               6, 17, 29, 16, 24>
            rng1(5, 4);
        adhoc::threefry_engine<std::uint64_t, 32, 2, 20, 0x1BD11BDA, 13, 15, 26,
                               6, 17, 29, 16, 24>
            rng2(5, 4);
        compare_rng(rng1, rng2, 100000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::threefry_engine<std::uint32_t, 32, 4, 20, 0x1BD11BDA, 10, 26, 11,
                               21, 13, 27, 23, 5, 6, 20, 17, 11, 25, 10, 18, 20>
            rng1(5, 4, 3, 2);
        adhoc::threefry_engine<std::uint64_t, 32, 4, 20, 0x1BD11BDA, 10, 26, 11,
                               21, 13, 27, 23, 5, 6, 20, 17, 11, 25, 10, 18, 20>
            rng2(5, 4, 3, 2);
        compare_rng(rng1, rng2, 100000);
        compare_rng_limits(rng1, rng2);
    }

    static_assert(std::bidirectional_iterator<adhoc::threefry2x32>);
    static_assert(std::bidirectional_iterator<adhoc::threefry4x32>);
    static_assert(std::bidirectional_iterator<adhoc::threefry2x64>);
    static_assert(std::bidirectional_iterator<adhoc::threefry4x64>);

    std::size_t sims = 0;
    if (auto env_p = std::getenv("TIMING_SIMS")) {
        sims = std::atoi(env_p);
    }

    if (sims) {
        std::uint64_t res_ref = 0;
        using namespace r123;
        Threefry4x64 rng_ref;
        Threefry4x64::ctr_type c = {{0, 0, 0, 0}};
        Threefry4x64::key_type k = {{5, 4, 3, 2}};
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
        adhoc::threefry4x64 rng(5, 4, 3, 2);
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
