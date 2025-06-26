#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include <linear_congruential_engine.hpp>
#include <romu_engine.hpp>

#include "external/romu/romu.hpp"
#include "read_csv.hpp"

#include <algorithm>
#include <chrono>
#include <ranges>

auto main() -> int {

    {
        auto vals = adhoc::readCSV<std::uint64_t>("./randomgen/romu_state.txt");
        adhoc::romuquad rng(vals[0], vals[1], vals[2], vals[3]);
        auto values_from_python =
            adhoc::readCSV<std::uint64_t>("./randomgen/romu_vals.txt");

        for (auto val1 : values_from_python) {
            auto val2 = rng();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        auto vals =
            adhoc::readCSV<std::uint64_t>("./randomgen/romu_trio_state.txt");
        adhoc::romutrio rng(vals[0], vals[1], vals[2]);
        auto values_from_python =
            adhoc::readCSV<std::uint64_t>("./randomgen/romu_trio_vals.txt");

        for (auto val1 : values_from_python) {
            auto val2 = rng();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        romuQuad_init(123, 456, 789, 357);
        adhoc::romuquad rng(123, 456, 789, 357);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = romuQuad_random();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        romuQuad32_init(123, 456, 789, 357);
        adhoc::romuquad32 rng(123, 456, 789, 357);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = romuQuad32_random();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        romuTrio_init(123, 456, 789);
        adhoc::romutrio rng(123, 456, 789);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = romuTrio_random();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        romuDuo_init(123, 456);
        adhoc::romuduo rng(123, 456);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = romuDuo_random();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        romuDuoJr_init(123, 456);
        adhoc::romuduojr rng(123, 456);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = romuDuoJr_random();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        romuTrio32_init(123, 456, 789);
        adhoc::romutrio32 rng(123, 456, 789);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = romuTrio32_random();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        std::uint32_t seed = 1234;
        romuMono32_init(seed);
        adhoc::romumono32 rng(seed);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = romuMono32_random();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::romuquad rng(123, 456, 789, 357);
        check_fwd_and_back(rng, 1000000);
        adhoc::romuquad rng2(123, 456, 789, 357);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::romuquad rng(123, 456, 789, 357);
        check_back_and_fwd(rng, 1000000);
        adhoc::romuquad rng2(123, 456, 789, 357);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::romutrio rng(123, 456, 789);
        check_fwd_and_back(rng, 1000000);
        adhoc::romutrio rng2(123, 456, 789);
        EXPECT_EQUAL(rng, rng2);
        bool stop = true;
    }

    {
        adhoc::romutrio rng(123, 456, 789);
        check_back_and_fwd(rng, 1000000);
        adhoc::romutrio rng2(123, 456, 789);
        EXPECT_EQUAL(rng, rng2);
        bool stop = true;
    }

    {
        adhoc::romuduo rng(123, 456);
        check_fwd_and_back(rng, 1000000);
        adhoc::romuduo rng2(123, 456);
        EXPECT_EQUAL(rng, rng2);
        bool stop = true;
    }

    {
        adhoc::romuduo rng(123, 456);
        check_back_and_fwd(rng, 1000000);
        adhoc::romuduo rng2(123, 456);
        EXPECT_EQUAL(rng, rng2);
        bool stop = true;
    }

    {
        adhoc::romuduojr rng(123, 456);
        check_fwd_and_back(rng, 1000000);
        adhoc::romuduojr rng2(123, 456);
        EXPECT_EQUAL(rng, rng2);
        bool stop = true;
    }

    {
        adhoc::romuduojr rng(123, 456);
        check_back_and_fwd(rng, 1000000);
        adhoc::romuduojr rng2(123, 456);
        EXPECT_EQUAL(rng, rng2);
        bool stop = true;
    }

    {
        adhoc::romuquad32 rng(123, 456, 789, 357);
        check_fwd_and_back(rng, 1000000);
        adhoc::romuquad32 rng2(123, 456, 789, 357);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::romuquad32 rng(123, 456, 789, 357);
        check_back_and_fwd(rng, 1000000);
        adhoc::romuquad32 rng2(123, 456, 789, 357);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::romutrio32 rng(123, 456, 789);
        check_fwd_and_back(rng, 1000000);
        adhoc::romutrio32 rng2(123, 456, 789);
        EXPECT_EQUAL(rng, rng2);
        bool stop = true;
    }

    {
        adhoc::romutrio32 rng(123, 456, 789);
        check_back_and_fwd(rng, 1000000);
        adhoc::romutrio32 rng2(123, 456, 789);
        EXPECT_EQUAL(rng, rng2);
        bool stop = true;
    }

    {
        adhoc::romumono32 rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::romumono32 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::romumono32 rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::romumono32 rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    static_assert(std::bidirectional_iterator<adhoc::romuquad>);

    {
        adhoc::romuquad rng(123, 456, 789, 357);
        check_discard(rng, 100);
    }

    {
        adhoc::romuquad rng1(123, 456, 789, 357);
        std::vector<std::uint64_t> v(10);
        std::ranges::generate(v, rng1);

        adhoc::romuquad rng2(123, 456, 789, 357);
        std::vector<std::uint64_t> v2(10);

        for (auto &val : v2) {
            val = *rng2;
            ++rng2;
        }

        for (std::size_t i = 0; i < v.size(); ++i) {
            EXPECT_EQUAL(v[i], v2[i]);
        }
    }

    {
        adhoc::romuquad rng(123, 456, 789, 357);
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
        adhoc::romuquad rng(123, 456, 789, 357);
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
