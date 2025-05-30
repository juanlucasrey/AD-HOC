#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "linear_congruential_engine.hpp"

#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>

int main() {
    // std check
    {
        adhoc::minstd_rand0 gen1; // overload (1)
        adhoc::minstd_rand gen2;  // overload (1)
        gen1.discard(10000 - 1);
        gen2.discard(10000 - 1);
        EXPECT_EQUAL(gen1(), 1043618065);
        EXPECT_EQUAL(gen2(), 399268537);
    }

    // check against std minstd_rand0
    {
        if constexpr (std::numeric_limits<std::uint_fast32_t>::digits == 32) {
            std::minstd_rand0 rng1;
            adhoc::minstd_rand0 rng2;
            compare_rng(rng1, rng2, 100);
            compare_rng_limits(rng1, rng2);
        } else {
            std::minstd_rand0 rng1;
            adhoc::linear_congruential_engine<std::uint64_t, 64, 16807, 0,
                                              2147483647>
                rng2;
            compare_rng(rng1, rng2, 100);
            compare_rng_limits(rng1, rng2);
        }
    }

    {
        adhoc::minstd_rand0 rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::minstd_rand0 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::minstd_rand0 rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::minstd_rand0 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        if constexpr (std::numeric_limits<std::uint_fast32_t>::digits == 32) {
            std::minstd_rand rng1;
            adhoc::minstd_rand rng2;
            compare_rng(rng1, rng2, 100);
            compare_rng_limits(rng1, rng2);
        } else {
            std::minstd_rand rng1;
            adhoc::linear_congruential_engine<std::uint64_t, 64, 48271, 0,
                                              2147483647>
                rng2;
            compare_rng(rng1, rng2, 100);
            compare_rng_limits(rng1, rng2);
        }
    }

    {
        adhoc::minstd_rand rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::minstd_rand rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::minstd_rand rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::minstd_rand rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        std::linear_congruential_engine<std::uint16_t, 75, 74, 0> rng1;
        adhoc::ZX81 rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::ZX81 rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::ZX81 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::ZX81 rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::ZX81 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        std::linear_congruential_engine<std::uint32_t, 1664525, 1013904223, 0>
            rng1;
        adhoc::ranqd1 rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::ranqd1 rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::ranqd1 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::ranqd1 rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::ranqd1 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        std::linear_congruential_engine<std::uint32_t, 65539, 0, 2147483648>
            rng1;
        adhoc::RANDU rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::RANDU rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::RANDU rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::RANDU rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::RANDU rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        std::linear_congruential_engine<std::uint32_t, 22695477, 1, 2147483648>
            rng1;
        adhoc::Borland rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::Borland rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::Borland rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::Borland rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::Borland rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        std::linear_congruential_engine<std::uint64_t, 6364136223846793005, 1,
                                        9223372036854775808U>
            rng1;
        adhoc::Newlib rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::Newlib rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::Newlib rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::Newlib rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::Newlib rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::linear_congruential_engine<std::uint32_t, 32, 22695477, 1,
                                          2147483648>
            rng1;
        adhoc::linear_congruential_engine<std::uint64_t, 32, 22695477, 1,
                                          2147483648>
            rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    using whatever =
        adhoc::linear_congruential_engine<std::uint_fast32_t, 32, 1664525,
                                          1013904223, 2147483647>;

    {
        whatever rng;
        check_fwd_and_back(rng, 1000000);
        whatever rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        whatever rng;
        check_back_and_fwd(rng, 1000000);
        whatever rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    int sims = 0;
    if (auto env_p = std::getenv("TIMING_SIMS")) {
        sims = std::atoi(env_p);
    }

    if (sims) {
        std::uint64_t res1 = 0;
        adhoc::minstd_rand0 rng;
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

        std::uint64_t res2 = 0;
        std::minstd_rand0 rng2;
        {
            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims; i++) {
                res2 += rng2();
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();
            std::cout << "old" << std::endl;
            std::cout << time << std::endl;
        }
        EXPECT_EQUAL(res1, res2);

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
