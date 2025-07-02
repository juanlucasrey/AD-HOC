#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "read_csv.hpp"

#include "seed_seq_inserter.hpp"
#include <sfmt_engine.hpp>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <iostream>
#include <ranges>
#include <vector>

namespace {

auto count_digits(const std::string &s) -> int {
    return std::count_if(
        s.begin(), s.end(),
        [](unsigned char c) { return std::isdigit(c); } // correct
    );
}

template <class UIntType>
auto readSSV(const std::string &filename) -> std::vector<UIntType> {
    std::vector<UIntType> data;
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
            if (count_digits(value) > 4) {
                UIntType number = std::stoull(value);
                data.push_back(number);
            }
        }
    }

    file.close();
    return data;
}

template <class RNG, class Arg> class range_rng final {
  public:
    range_rng(Arg &s) : arg(s) {}
    range_rng(Arg &&s) : arg(s) {}

    [[nodiscard]] auto begin() -> RNG { return RNG(arg); }
    [[nodiscard]] static auto end() -> std::default_sentinel_t { return {}; }

  private:
    Arg arg;
};

} // namespace

auto main() -> int {

    {
        adhoc::seed_seq_inserter seq(
            adhoc::readCSV<std::uint32_t>("./randomgen/SFMT_state.txt"));

        auto const values_from_python = adhoc::split_uint64_to_uint32(
            adhoc::readCSV<std::uint64_t>("./randomgen/SFMT_vals.txt"));

        adhoc::sfmt19937 rng(seq);
        std::size_t j = 0;
        for (std::size_t i = 0; i < values_from_python.size(); ++i) {
            auto val1 = rng();
            auto val2 = values_from_python[j++];
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::seed_seq_inserter seq(
            adhoc::readCSV<std::uint32_t>("./randomgen/SFMT_state.txt"));

        auto const values_from_python =
            adhoc::readCSV<std::uint64_t>("./randomgen/SFMT_vals.txt");

        adhoc::sfmt19937_64 rng(seq);
        std::size_t j = 0;
        for (std::size_t i = 0; i < values_from_python.size(); ++i) {
            auto val1 = rng();
            auto val2 = values_from_python[j++];
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        auto vals = readSSV<std::uint32_t>("./SFMT-src-1.5.1/SFMT.607.out.txt");
        // range_rng<adhoc::sfmt607, std::uint32_t> range(1234);
        // EXPECT_TRUE(std::ranges::equal(std::views::take(vals, 1000),
        //                                std::views::take(range, 1000)));

        adhoc::sfmt607 rng(1234);
        std::size_t j = 0;
        for (std::size_t i = 0; i < 1000; ++i) {
            auto val1 = rng();
            auto val2 = vals[j++];
            EXPECT_EQUAL(val1, val2);
        }

        std::array<uint32_t, 4> key{0x1234, 0x5678, 0x9abc, 0xdef0};
        adhoc::sfmt607 rng2(key);
        for (std::size_t i = 0; i < 1000; ++i) {
            auto val1 = rng2();
            auto val2 = vals[j++];
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        std::size_t j = 0;
        auto vals =
            readSSV<std::uint32_t>("./SFMT-src-1.5.1/SFMT.1279.out.txt");
        adhoc::sfmt1279 rng(1234);
        for (std::size_t i = 0; i < 1000; ++i) {
            auto val1 = rng();
            auto val2 = vals[j++];
            EXPECT_EQUAL(val1, val2);
        }

        std::array<uint32_t, 4> key{0x1234, 0x5678, 0x9abc, 0xdef0};
        adhoc::sfmt1279 rng2(key);
        for (std::size_t i = 0; i < 1000; ++i) {
            auto val1 = rng2();
            auto val2 = vals[j++];
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        std::size_t j = 0;
        auto vals =
            readSSV<std::uint32_t>("./SFMT-src-1.5.1/SFMT.2281.out.txt");
        adhoc::sfmt2281 rng(1234);
        for (std::size_t i = 0; i < 1000; ++i) {
            auto val1 = rng();
            auto val2 = vals[j++];
            EXPECT_EQUAL(val1, val2);
        }

        std::array<uint32_t, 4> key{0x1234, 0x5678, 0x9abc, 0xdef0};
        adhoc::sfmt2281 rng2(key);
        for (std::size_t i = 0; i < 1000; ++i) {
            auto val1 = rng2();
            auto val2 = vals[j++];
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        std::size_t j = 0;
        auto vals =
            readSSV<std::uint32_t>("./SFMT-src-1.5.1/SFMT.4253.out.txt");
        adhoc::sfmt4253 rng(1234);
        for (std::size_t i = 0; i < 1000; ++i) {
            auto val1 = rng();
            auto val2 = vals[j++];
            EXPECT_EQUAL(val1, val2);
        }

        std::array<uint32_t, 4> key{0x1234, 0x5678, 0x9abc, 0xdef0};
        adhoc::sfmt4253 rng2(key);
        for (std::size_t i = 0; i < 1000; ++i) {
            auto val1 = rng2();
            auto val2 = vals[j++];
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        std::size_t j = 0;
        auto vals =
            readSSV<std::uint32_t>("./SFMT-src-1.5.1/SFMT.11213.out.txt");
        adhoc::sfmt11213 rng(1234);
        for (std::size_t i = 0; i < 1000; ++i) {
            auto val1 = rng();
            auto val2 = vals[j++];
            EXPECT_EQUAL(val1, val2);
        }

        std::array<uint32_t, 4> key{0x1234, 0x5678, 0x9abc, 0xdef0};
        adhoc::sfmt11213 rng2(key);
        for (std::size_t i = 0; i < 1000; ++i) {
            auto val1 = rng2();
            auto val2 = vals[j++];
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        std::size_t j = 0;
        auto vals =
            readSSV<std::uint32_t>("./SFMT-src-1.5.1/SFMT.19937.out.txt");
        adhoc::sfmt19937 rng(1234);
        for (std::size_t i = 0; i < 1000; ++i) {
            auto val1 = rng();
            auto val2 = vals[j++];
            EXPECT_EQUAL(val1, val2);
        }

        std::array<uint32_t, 4> key{0x1234, 0x5678, 0x9abc, 0xdef0};
        adhoc::sfmt19937 rng2(key);
        for (std::size_t i = 0; i < 1000; ++i) {
            auto val1 = rng2();
            auto val2 = vals[j++];
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        std::size_t j = 0;
        auto vals =
            readSSV<std::uint32_t>("./SFMT-src-1.5.1/SFMT.44497.out.txt");
        adhoc::sfmt44497 rng(1234);
        for (std::size_t i = 0; i < 1000; ++i) {
            auto val1 = rng();
            auto val2 = vals[j++];
            EXPECT_EQUAL(val1, val2);
        }

        std::array<uint32_t, 4> key{0x1234, 0x5678, 0x9abc, 0xdef0};
        adhoc::sfmt44497 rng2(key);
        for (std::size_t i = 0; i < 1000; ++i) {
            auto val1 = rng2();
            auto val2 = vals[j++];
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        std::size_t j = 0;
        auto vals =
            readSSV<std::uint32_t>("./SFMT-src-1.5.1/SFMT.86243.out.txt");
        adhoc::sfmt86243 rng(1234);
        for (std::size_t i = 0; i < 1000; ++i) {
            auto val1 = rng();
            auto val2 = vals[j++];
            EXPECT_EQUAL(val1, val2);
        }

        std::array<uint32_t, 4> key{0x1234, 0x5678, 0x9abc, 0xdef0};
        adhoc::sfmt86243 rng2(key);
        for (std::size_t i = 0; i < 1000; ++i) {
            auto val1 = rng2();
            auto val2 = vals[j++];
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        std::size_t j = 0;
        auto vals =
            readSSV<std::uint32_t>("./SFMT-src-1.5.1/SFMT.132049.out.txt");
        adhoc::sfmt132049 rng(1234);
        for (std::size_t i = 0; i < 1000; ++i) {
            auto val1 = rng();
            auto val2 = vals[j++];
            EXPECT_EQUAL(val1, val2);
        }

        std::array<uint32_t, 4> key{0x1234, 0x5678, 0x9abc, 0xdef0};
        adhoc::sfmt132049 rng2(key);
        for (std::size_t i = 0; i < 1000; ++i) {
            auto val1 = rng2();
            auto val2 = vals[j++];
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        std::size_t j = 0;
        auto vals =
            readSSV<std::uint32_t>("./SFMT-src-1.5.1/SFMT.216091.out.txt");
        adhoc::sfmt216091 rng(1234);
        for (std::size_t i = 0; i < 1000; ++i) {
            auto val1 = rng();
            auto val2 = vals[j++];
            EXPECT_EQUAL(val1, val2);
        }

        std::array<uint32_t, 4> key{0x1234, 0x5678, 0x9abc, 0xdef0};
        adhoc::sfmt216091 rng2(key);
        for (std::size_t i = 0; i < 1000; ++i) {
            auto val1 = rng2();
            auto val2 = vals[j++];
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::sfmt607 rng(1234);
        check_fwd_and_back(rng, 1000000);
        adhoc::sfmt607 rng2(1234);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::sfmt607 rng(1234);
        check_back_and_fwd(rng, 1000000);
        adhoc::sfmt607 rng2(1234);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::sfmt607 rng(1234);
        check_discard(rng, 100);
    }

    {
        adhoc::sfmt_engine<std::uint64_t, 64, 607, 2, 15, 24U, 13, 24U,
                           0xfdff37ffU, 0xef7f3f7dU, 0xff777b7dU, 0x7ff7fb2fU,
                           0x00000001U, 0x00000000U, 0x00000000U, 0x5986f054U>
            rng1(12);

        adhoc::sfmt_engine<adhoc::uint128, 64, 607, 2, 15, 24U, 13, 24U,
                           0xfdff37ffU, 0xef7f3f7dU, 0xff777b7dU, 0x7ff7fb2fU,
                           0x00000001U, 0x00000000U, 0x00000000U, 0x5986f054U>
            rng2(12);

        compare_rng(rng1, rng2, 10000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::sfmt_engine<std::uint64_t, 64, 607, 2, 15, 24U, 13, 24U,
                           0xfdff37ffU, 0xef7f3f7dU, 0xff777b7dU, 0x7ff7fb2fU,
                           0x00000001U, 0x00000000U, 0x00000000U, 0x5986f054U>
            rng(12);
        check_fwd_and_back(rng, 1000000);
        adhoc::sfmt_engine<std::uint64_t, 64, 607, 2, 15, 24U, 13, 24U,
                           0xfdff37ffU, 0xef7f3f7dU, 0xff777b7dU, 0x7ff7fb2fU,
                           0x00000001U, 0x00000000U, 0x00000000U, 0x5986f054U>
            rng2(12);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::sfmt_engine<std::uint64_t, 64, 607, 2, 15, 24U, 13, 24U,
                           0xfdff37ffU, 0xef7f3f7dU, 0xff777b7dU, 0x7ff7fb2fU,
                           0x00000001U, 0x00000000U, 0x00000000U, 0x5986f054U>
            rng(12);
        check_back_and_fwd(rng, 1000000);
        adhoc::sfmt_engine<std::uint64_t, 64, 607, 2, 15, 24U, 13, 24U,
                           0xfdff37ffU, 0xef7f3f7dU, 0xff777b7dU, 0x7ff7fb2fU,
                           0x00000001U, 0x00000000U, 0x00000000U, 0x5986f054U>
            rng2(12);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::sfmt_engine<std::uint32_t, 32, 607, 2, 15, 24U, 13, 24U,
                           0xfdff37ffU, 0xef7f3f7dU, 0xff777b7dU, 0x7ff7fb2fU,
                           0x00000001U, 0x00000000U, 0x00000000U, 0x5986f054U>
            rng1(12);

        adhoc::sfmt_engine<std::uint64_t, 32, 607, 2, 15, 24U, 13, 24U,
                           0xfdff37ffU, 0xef7f3f7dU, 0xff777b7dU, 0x7ff7fb2fU,
                           0x00000001U, 0x00000000U, 0x00000000U, 0x5986f054U>
            rng2(12);

        compare_rng(rng1, rng2, 10000);
        compare_rng_limits(rng1, rng2);
    }

    {
        adhoc::sfmt_engine<std::uint64_t, 32, 607, 2, 15, 24U, 13, 24U,
                           0xfdff37ffU, 0xef7f3f7dU, 0xff777b7dU, 0x7ff7fb2fU,
                           0x00000001U, 0x00000000U, 0x00000000U, 0x5986f054U>
            rng(12);
        check_fwd_and_back(rng, 1000000);
        adhoc::sfmt_engine<std::uint64_t, 32, 607, 2, 15, 24U, 13, 24U,
                           0xfdff37ffU, 0xef7f3f7dU, 0xff777b7dU, 0x7ff7fb2fU,
                           0x00000001U, 0x00000000U, 0x00000000U, 0x5986f054U>
            rng2(12);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::sfmt_engine<std::uint64_t, 32, 607, 2, 15, 24U, 13, 24U,
                           0xfdff37ffU, 0xef7f3f7dU, 0xff777b7dU, 0x7ff7fb2fU,
                           0x00000001U, 0x00000000U, 0x00000000U, 0x5986f054U>
            rng(12);
        check_back_and_fwd(rng, 1000000);
        adhoc::sfmt_engine<std::uint64_t, 32, 607, 2, 15, 24U, 13, 24U,
                           0xfdff37ffU, 0xef7f3f7dU, 0xff777b7dU, 0x7ff7fb2fU,
                           0x00000001U, 0x00000000U, 0x00000000U, 0x5986f054U>
            rng2(12);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        std::array<uint32_t, 4> key{0x1234, 0x5678, 0x9abc, 0xdef0};
        adhoc::sfmt607_64 rng(key);
        adhoc::sfmt607 rng2(key);

        for (std::size_t i = 0; i < 100; ++i) {
            std::uint64_t val1a = rng2();
            std::uint64_t val1b = rng2();
            std::uint64_t val1 = val1a | (val1b << 32U);
            std::uint64_t val2 = rng();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::sfmt607_64 rng(12);
        adhoc::sfmt607 rng2(12);

        for (std::size_t i = 0; i < 100; ++i) {
            std::uint64_t val1a = rng2();
            std::uint64_t val1b = rng2();
            std::uint64_t val1 = val1a | (val1b << 32U);
            std::uint64_t val2 = rng();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::sfmt607_64 rng(12);
        check_fwd_and_back(rng, 1000000);
        adhoc::sfmt607_64 rng2(12);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::sfmt607_64 rng(12);
        check_back_and_fwd(rng, 1000000);
        adhoc::sfmt607_64 rng2(12);
        EXPECT_EQUAL(rng, rng2);
    }

    static_assert(std::bidirectional_iterator<adhoc::sfmt607>);
    static_assert(std::bidirectional_iterator<adhoc::sfmt1279>);
    static_assert(std::bidirectional_iterator<adhoc::sfmt2281>);
    static_assert(std::bidirectional_iterator<adhoc::sfmt4253>);
    static_assert(std::bidirectional_iterator<adhoc::sfmt11213>);
    static_assert(std::bidirectional_iterator<adhoc::sfmt19937>);
    static_assert(std::bidirectional_iterator<adhoc::sfmt44497>);
    static_assert(std::bidirectional_iterator<adhoc::sfmt86243>);
    static_assert(std::bidirectional_iterator<adhoc::sfmt132049>);
    static_assert(std::bidirectional_iterator<adhoc::sfmt216091>);

    {
        adhoc::sfmt607 rng1(1234);
        std::vector<std::uint32_t> v(10);
        std::ranges::generate(v, rng1);

        adhoc::sfmt607 rng2(1234);
        std::vector<std::uint32_t> v2(10);

        for (std::size_t i = 0; i < v2.size(); ++i, ++rng2) {
            v2[i] = *rng2;
        }

        for (std::size_t i = 0; i < v.size(); ++i) {
            EXPECT_EQUAL(v[i], v2[i]);
        }
    }

    {
        adhoc::sfmt607 rng(1234);
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
        {
            std::uint64_t res1 = 0;
            adhoc::sfmt607 rng(1234);
            {
                auto time1 = std::chrono::high_resolution_clock::now();
                for (std::size_t i = 0; i < sims; i++) {
                    res1 += rng();
                }
                auto time2 = std::chrono::high_resolution_clock::now();
                auto time =
                    std::chrono::duration_cast<std::chrono::milliseconds>(
                        time2 - time1)
                        .count();
                std::cout << "fwd 32" << std::endl;
                std::cout << time << std::endl;
            }

            std::uint64_t res1b = 0;
            {
                auto time1 = std::chrono::high_resolution_clock::now();
                for (std::size_t i = 0; i < sims; i++) {
                    res1b += rng.operator()<false>();
                }
                auto time2 = std::chrono::high_resolution_clock::now();
                auto time =
                    std::chrono::duration_cast<std::chrono::milliseconds>(
                        time2 - time1)
                        .count();
                std::cout << "back 32" << std::endl;
                std::cout << time << std::endl;
            }
            EXPECT_EQUAL(res1, res1b);
        }

        {
            std::uint64_t res1 = 0;
            adhoc::sfmt607_64 rng(1234);
            {
                auto time1 = std::chrono::high_resolution_clock::now();
                for (std::size_t i = 0; i < sims; i++) {
                    res1 += rng();
                }
                auto time2 = std::chrono::high_resolution_clock::now();
                auto time =
                    std::chrono::duration_cast<std::chrono::milliseconds>(
                        time2 - time1)
                        .count();
                std::cout << "fwd 64" << std::endl;
                std::cout << time << std::endl;
            }

            std::uint64_t res1b = 0;
            {
                auto time1 = std::chrono::high_resolution_clock::now();
                for (std::size_t i = 0; i < sims; i++) {
                    res1b += rng.operator()<false>();
                }
                auto time2 = std::chrono::high_resolution_clock::now();
                auto time =
                    std::chrono::duration_cast<std::chrono::milliseconds>(
                        time2 - time1)
                        .count();
                std::cout << "back 64" << std::endl;
                std::cout << time << std::endl;
            }
            EXPECT_EQUAL(res1, res1b);
        }
    }

    TEST_END;
}
