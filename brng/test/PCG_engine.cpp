#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "external/pcg_basic.h"

#include "external/pcg-cpp/pcg_random.hpp"

#include <PCG_engine.hpp>

#include <tools/mask.hpp>

namespace {

constexpr auto operator""_ULLL(const char *digits) -> __uint128_t {
    __uint128_t result{};

    while (*digits != 0) {
        result *= 10;
        result += *digits - '0';
        ++digits;
    }

    return result;
}

} // namespace

auto main() -> int {

    std::size_t const sims = 1000000;
    {
        adhoc::PCG_engine rng(1955588763U, 2235320806U, 2495175643U,
                              3661511115U);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng();
            auto val2 = pcg32_random();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32 rng1;
        adhoc::PCG_engine rng2;

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_oneseq rng1;
        adhoc::PCG_engine rng2;

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_fast rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rs>
            rng2(3512640999U, 3405705229U, 0U, 0U);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg64 rng1;
        adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::xsl_rr>
            rng2;

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg64_oneseq rng1;
        adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::xsl_rr>
            rng2;

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg64_fast rng1;
        adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::xsl_rr>
            rng2{14627392581883831783ULL, 0, 0, 0};

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_once_insecure rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::rxs_m_xs,
                          false>
            rng2{1186293367UL, 2891336453UL};

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg64_once_insecure rng1;
        adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::rxs_m_xs,
                          false>
            rng2{5573589319906701683ULL, 1442695040888963407ULL};

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg128_once_insecure rng1;
        adhoc::PCG_engine<__uint128_t, 128, adhoc::tempering_type::xsl_rr_rr,
                          false>
            rng2(245720598905631564143578724636268694099_ULLL,
                 117397592171526113268558934119004209487_ULLL);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_oneseq_once_insecure rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::rxs_m_xs,
                          false>
            rng2(1186293367UL, 2891336453UL);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_k2 rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rr,
                          true, 1, 16>
            rng2(5573589319906701683ULL, 1442695040888963407ULL);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_k2_fast rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rs,
                          true, 1, 32>
            rng2(5573589319906701683ULL, 1442695040888963407ULL);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_k64 rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rr,
                          true, 6, 16>
            rng2(5573589319906701683ULL, 1442695040888963407ULL);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_k64_oneseq rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rs,
                          true, 6, 32, true, true>
            rng2(14627392581883831783ULL, 0);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_k64_fast rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rs,
                          true, 6, 32>
            rng2(5573589319906701683ULL, 1442695040888963407ULL);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_c64 rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rr,
                          true, 6, 16, false, false, false>
            rng2(5573589319906701683ULL, 1442695040888963407ULL);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_c64_oneseq rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rs,
                          true, 6, 32, false, false, false>
            rng2(5573589319906701683ULL, 1442695040888963407ULL);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_c64_fast rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rs,
                          true, 6, 32, false, false, false>
            rng2(14627392581883831783ULL, 0);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg64_k32 rng1;
        adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::xsl_rr,
                          true, 5, 16>
            rng2;

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg64_k32_oneseq rng1;
        adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::xsl_rr,
                          true, 5, 128>
            rng2;

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg64_k32_fast rng1;
        adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::xsl_rr,
                          true, 5, 128, false, true>
            rng2{14627392581883831783ULL, 0, 0, 0};

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg64_c32 rng1;
        adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::xsl_rr,
                          true, 5, 16, false, true, false>
            rng2;

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg64_c32_oneseq rng1;
        adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::xsl_rr,
                          true, 5, 128, false, true, false>
            rng2;

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg64_c32_fast rng1;
        adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::xsl_rr,
                          true, 5, 128, false, true, false>
            rng2{14627392581883831783ULL, 0, 0, 0};

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_k1024 rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rr,
                          true, 10, 16>
            rng2(5573589319906701683ULL, 1442695040888963407ULL);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_k1024_fast rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rs,
                          true, 10, 32>
            rng2(5573589319906701683ULL, 1442695040888963407ULL);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_c1024 rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rr,
                          true, 10, 16, false, false, false>
            rng2(5573589319906701683ULL, 1442695040888963407ULL);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_c1024_fast rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rs,
                          true, 10, 32, false, false, false>
            rng2(5573589319906701683ULL, 1442695040888963407ULL);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg64_k1024 rng1;
        adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::xsl_rr,
                          true, 10, 16>
            rng2;

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg64_k1024_fast rng1;
        adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::xsl_rr,
                          true, 10, 128>
            rng2;

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg64_c1024 rng1;
        adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::xsl_rr,
                          true, 10, 16, false, true, false>
            rng2;

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg64_c1024_fast rng1;
        adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::xsl_rr,
                          true, 10, 128, false, true, false>
            rng2;

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_k16384 rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rr,
                          true, 14, 16>
            rng2(5573589319906701683ULL, 1442695040888963407ULL);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_k16384_fast rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rs,
                          true, 14, 32>
            rng2(5573589319906701683ULL, 1442695040888963407ULL);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }
}
