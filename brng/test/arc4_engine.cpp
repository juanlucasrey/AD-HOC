#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "external/arc4.hpp"
#include <arc4_engine.hpp>

#include <chrono>
#include <iostream>
#include <random>

auto main() -> int {
    {
        arc4_rand rng;

        std::seed_seq init{0xd0beef41U};
        adhoc::arc4_rand rng2(init);

        for (std::size_t i = 0; i < 20; ++i) {
            auto val1 = rng();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        std::seed_seq init{0xd0beef41U};
        adhoc::arc4_rand rng(init);
        check_fwd_and_back(rng, 1000000);
        std::seed_seq init2{0xd0beef41U};
        adhoc::arc4_rand rng2(init2);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        std::seed_seq init{0xd0beef41U};
        adhoc::arc4_rand rng(init);
        check_back_and_fwd(rng, 1000000);
        std::seed_seq init2{0xd0beef41U};
        adhoc::arc4_rand rng2(init2);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        std::seed_seq init{0xd0beef41U};
        adhoc::arc4_rand64 rng(init);
        check_fwd_and_back(rng, 1000000);
        std::seed_seq init2{0xd0beef41U};
        adhoc::arc4_rand64 rng2(init2);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        std::seed_seq init{0xd0beef41U};
        adhoc::arc4_rand64 rng(init);
        check_back_and_fwd(rng, 1000000);
        std::seed_seq init2{0xd0beef41U};
        adhoc::arc4_rand64 rng2(init2);
        EXPECT_EQUAL(rng, rng2);
    }

    std::size_t sims = 0;
    if (auto env_p = std::getenv("TIMING_SIMS")) {
        sims = std::atoi(env_p);
    }

    if (sims) {
        std::uint64_t res1 = 0;
        std::seed_seq init{0xd0beef41U};
        adhoc::arc4_rand rng(init);
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

        std::uint64_t res2 = 0;
        {
            arc4_rand rng;

            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims; i++) {
                res2 += rng();
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();

            std::cout << "old" << std::endl;
            std::cout << time << std::endl;
        }
        // EXPECT_NEAR_REL(res1, res2, 1e-12);
        EXPECT_EQUAL(res1, res2);

        std::uint64_t res3 = 0;
        std::seed_seq init3{0xd0beef41U};
        adhoc::arc4_rand rng3(init3);
        {
            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims; i++) {
                res3 += rng3();
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();
            std::cout << "new" << std::endl;
            std::cout << time << std::endl;
        }
        EXPECT_EQUAL(res1, res3);
    }

    TEST_END;
}
