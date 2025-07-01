#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "read_csv.hpp"

#include "seed_seq_inserter.hpp"
#include <dsfmt_engine.hpp>
#include <tools/uint128.hpp>

#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include <distribution/centered_canonical_distribution.hpp>
#include <tools/mask.hpp>

namespace {

template <class Float = double>
auto readSSV(const std::string &filename) -> std::vector<Float> {
    std::vector<Float> data;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return data;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;

        while (std::getline(ss, value, ' ')) {
            if (value.find('.') != std::string::npos) {
                Float number = std::stod(value);
                data.push_back(number);
            }
        }
    }

    file.close();
    return data;
}

template <class RNG> auto test_vals(std::string name, double tolerance) {
    {
        std::size_t i = 0;
        auto results = readSSV<double>(name);

        {
            adhoc::centered_canonical_default<RNG> rng(
                static_cast<std::uint32_t>(0));

            for (std::size_t j = 0; j < 1000; ++j) {
                auto val = rng();

                double val1 = val > 0 ? (val + 1.0) : (val + 2.0);
                double val2 = results[i++];
                auto [is_near, tol] = expect_near_rel(val1, val2, tolerance);
                if (!is_near) {
                    bool stop = true;
                }
                EXPECT_NEAR_REL(val1, val2, tolerance);
            }
        }

        std::uint32_t seed = 0;
        for (std::size_t rep = 0; rep < 5; ++rep) {
            {
                adhoc::centered_canonical_default<RNG> rng(
                    static_cast<std::uint32_t>(seed++));

                for (std::size_t j = 0; j < 12; ++j) {
                    auto val = rng();

                    double val1 = val > 0 ? val : (val + 1.0);
                    double val2 = results[i++];
                    EXPECT_NEAR_REL(val1, val2, tolerance);
                }
            }

            {
                adhoc::centered_canonical_default<RNG> rng(
                    static_cast<std::uint32_t>(seed++));

                for (std::size_t j = 0; j < 12; ++j) {
                    auto val = rng();

                    double val1 = 1.0 - (val > 0 ? val : (val + 1.0));
                    double val2 = results[i++];
                    EXPECT_NEAR_REL(val1, val2, tolerance);
                }
            }

            {
                adhoc::centered_canonical_open<RNG> rng(
                    static_cast<std::uint32_t>(seed++));

                for (std::size_t j = 0; j < 12; ++j) {
                    auto val = rng();

                    double val1 = val > 0 ? val : (val + 1.0);
                    double val2 = results[i++];
                    EXPECT_NEAR_REL(val1, val2, tolerance);
                }
            }

            {
                adhoc::centered_canonical_default<RNG> rng(
                    static_cast<std::uint32_t>(seed++));

                for (std::size_t j = 0; j < 12; ++j) {
                    auto val = rng();

                    double val1 = val > 0 ? val : (val + 1.0);
                    double val2 = results[i++];
                    EXPECT_NEAR_REL(val1, val2, tolerance);
                }
            }
        }

        {
            std::array<std::uint32_t, 4> ar{1, 2, 3, 4};
            adhoc::centered_canonical_default<RNG> rng(ar);

            for (std::size_t j = 0; j < 1000; ++j) {
                auto val = rng();

                double val1 = val > 0 ? (val + 1.0) : (val + 2.0);
                double val2 = results[i++];
                EXPECT_NEAR_REL(val1, val2, tolerance);
            }
        }
    }
}
} // namespace

auto main() -> int {
    test_vals<adhoc::dsfmt521>("./dSFMT-src-2.2.3/dSFMT.521.out.txt", 5e-14);
    test_vals<adhoc::dsfmt1279>("./dSFMT-src-2.2.3/dSFMT.1279.out.txt", 5e-14);
    test_vals<adhoc::dsfmt2203>("./dSFMT-src-2.2.3/dSFMT.2203.out.txt", 2e-13);
    test_vals<adhoc::dsfmt4253>("./dSFMT-src-2.2.3/dSFMT.4253.out.txt", 5e-14);
    test_vals<adhoc::dsfmt11213>("./dSFMT-src-2.2.3/dSFMT.11213.out.txt",
                                 6e-14);
    test_vals<adhoc::dsfmt19937>("./dSFMT-src-2.2.3/dSFMT.19937.out.txt",
                                 7e-14);
    test_vals<adhoc::dsfmt44497>("./dSFMT-src-2.2.3/dSFMT.44497.out.txt",
                                 5e-14);
    test_vals<adhoc::dsfmt86243>("./dSFMT-src-2.2.3/dSFMT.86243.out.txt",
                                 5e-14);
    test_vals<adhoc::dsfmt132049>("./dSFMT-src-2.2.3/dSFMT.132049.out.txt",
                                  5e-14);
    test_vals<adhoc::dsfmt216091>("./dSFMT-src-2.2.3/dSFMT.216091.out.txt",
                                  5e-14);

    {
        std::vector<std::uint64_t> values_uint =
            adhoc::readCSV<std::uint64_t>("./randomgen/DSFMT_vals.txt");

        adhoc::seed_seq_inserter seq(
            adhoc::readCSV<std::uint64_t>("./randomgen/DSFMT_state.txt"));
        adhoc::dsfmt19937 rng(seq);
        constexpr auto sign_and_exponent = adhoc::mask<std::uint64_t>(10, 52);
        std::size_t j = 0;
        for (std::size_t i = 0; i < values_uint.size(); ++i) {
            auto val1 = rng();

            auto val_union = val1 | sign_and_exponent;
            auto val2 = values_uint[j++];
            EXPECT_EQUAL(val_union, val2);
        }
    }

    {
        std::vector<double> values_double =
            readSSV<double>("./randomgen/DSFMT_double_vals.txt");

        adhoc::seed_seq_inserter seq(
            adhoc::readCSV<std::uint64_t>("./randomgen/DSFMT_state.txt"));
        adhoc::centered_canonical_default<adhoc::dsfmt19937> rng(seq);
        constexpr auto sign_and_exponent = adhoc::mask<std::uint64_t>(10, 52);
        std::size_t j = 0;
        for (std::size_t i = 0; i < values_double.size(); ++i) {
            auto val1 = rng();

            double val_shifted = val1 > 0 ? val1 : (val1 + 1.0);
            auto val2 = values_double[j++];
            EXPECT_EQUAL(val_shifted, val2);
        }
    }

    {
        adhoc::dsfmt521 rng(0);
        check_fwd_and_back(rng, 1000000);
        adhoc::dsfmt521 rng2(0);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::dsfmt521 rng(0);
        check_back_and_fwd(rng, 1000000);
        adhoc::dsfmt521 rng2(0);
        EXPECT_EQUAL(rng, rng2);
        bool tmp = rng != rng2;
        EXPECT_EQUAL(false, tmp);
    }

    {
        adhoc::dsfmt521 rng(0);
        check_discard(rng, 100);
    }

    {
        adhoc::dsfmt521 rng(0);
        auto rng2 = rng++;
        adhoc::dsfmt521 rng3(0);
        EXPECT_EQUAL(rng2, rng3);
        ++rng3;
        EXPECT_EQUAL(rng, rng3);
    }

    {
        adhoc::dsfmt521 rng(0);
        auto rng2 = rng--;
        adhoc::dsfmt521 rng3(0);
        EXPECT_EQUAL(rng2, rng3);
        --rng3;
        EXPECT_EQUAL(rng, rng3);
    }

    static_assert(std::bidirectional_iterator<adhoc::dsfmt521>);
    static_assert(std::bidirectional_iterator<adhoc::dsfmt1279>);
    static_assert(std::bidirectional_iterator<adhoc::dsfmt2203>);
    static_assert(std::bidirectional_iterator<adhoc::dsfmt4253>);
    static_assert(std::bidirectional_iterator<adhoc::dsfmt11213>);
    static_assert(std::bidirectional_iterator<adhoc::dsfmt19937>);
    static_assert(std::bidirectional_iterator<adhoc::dsfmt44497>);
    static_assert(std::bidirectional_iterator<adhoc::dsfmt86243>);
    static_assert(std::bidirectional_iterator<adhoc::dsfmt132049>);
    static_assert(std::bidirectional_iterator<adhoc::dsfmt216091>);

    {
        adhoc::dsfmt_engine<std::uint64_t, 521, 3, 25, 0x000fbfefff77efffUL,
                            0x000ffeebfbdfbfdfUL, 0xcfb393d661638469UL,
                            0xc166867883ae2adbUL, 0xccaa588000000000UL,
                            0x0000000000000001UL>
            rng1(12);

        adhoc::dsfmt_engine<adhoc::uint128, 521, 3, 25, 0x000fbfefff77efffUL,
                            0x000ffeebfbdfbfdfUL, 0xcfb393d661638469UL,
                            0xc166867883ae2adbUL, 0xccaa588000000000UL,
                            0x0000000000000001UL>
            rng2(12);

        compare_rng(rng1, rng2, 10000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::dsfmt_engine<adhoc::uint128, 521, 3, 25, 0x000fbfefff77efffUL,
                            0x000ffeebfbdfbfdfUL, 0xcfb393d661638469UL,
                            0xc166867883ae2adbUL, 0xccaa588000000000UL,
                            0x0000000000000001UL>
            rng(12);
        check_fwd_and_back(rng, 1000000);
        adhoc::dsfmt_engine<adhoc::uint128, 521, 3, 25, 0x000fbfefff77efffUL,
                            0x000ffeebfbdfbfdfUL, 0xcfb393d661638469UL,
                            0xc166867883ae2adbUL, 0xccaa588000000000UL,
                            0x0000000000000001UL>
            rng2(12);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::dsfmt_engine<adhoc::uint128, 521, 3, 25, 0x000fbfefff77efffUL,
                            0x000ffeebfbdfbfdfUL, 0xcfb393d661638469UL,
                            0xc166867883ae2adbUL, 0xccaa588000000000UL,
                            0x0000000000000001UL>
            rng(12);
        check_back_and_fwd(rng, 1000000);
        adhoc::dsfmt_engine<adhoc::uint128, 521, 3, 25, 0x000fbfefff77efffUL,
                            0x000ffeebfbdfbfdfUL, 0xcfb393d661638469UL,
                            0xc166867883ae2adbUL, 0xccaa588000000000UL,
                            0x0000000000000001UL>
            rng2(12);
        EXPECT_EQUAL(rng, rng2);
    }

    std::size_t sims = 0;
    if (auto env_p = std::getenv("TIMING_SIMS")) {
        sims = std::atoi(env_p);
    }

    if (sims) {
        double res1 = 0;
        adhoc::centered_canonical_default<adhoc::dsfmt521> rng(
            std::uint32_t(1234));
        {
            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims; i++) {
                auto const val = rng();
                res1 += val > 0 ? val : (val + 1.0);
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();
            std::cout << "fwd" << std::endl;
            std::cout << time << std::endl;
        }

        double res1b = 0;
        {
            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims; i++) {
                auto const val = rng.operator()<false>();
                res1b += val > 0 ? val : (val + 1.0);
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();
            std::cout << "back" << std::endl;
            std::cout << time << std::endl;
        }
        EXPECT_NEAR_REL(res1, res1b, 1e-12);
    }

    TEST_END;
}
