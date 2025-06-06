#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include "external/pcg_basic.h"

#include "external/pcg-cpp/pcg_random.hpp"

#include <PCG_engine.hpp>

#include <tools/mask.hpp>

auto main() -> int {

    std::size_t const sims = 1000000;
    {
        adhoc::PCG_engine rng(11906719003534950848UL, 15726070495360670683UL);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng();
            auto val2 = pcg32_random();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::PCG_engine rng(11906719003534950848UL, 15726070495360670683UL);
        check_fwd_and_back(rng, 1000000);
        adhoc::PCG_engine rng2(11906719003534950848UL, 15726070495360670683UL);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::PCG_engine rng(11906719003534950848UL, 15726070495360670683UL);
        check_back_and_fwd(rng, 1000000);
        adhoc::PCG_engine rng2(11906719003534950848UL, 15726070495360670683UL);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        pcg32 rng1;
        adhoc::PCG_engine rng2(16070087622772795188UL, 1442695040888963407UL);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_oneseq rng1;
        adhoc::PCG_engine rng2(16070087622772795188UL, 1442695040888963407UL);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_fast rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rs>
            rng2(7515307391605212643UL, 0UL);

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

    // investigate
    // {
    //     pcg32_once_insecure rng1;
    //     adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::rxs_m_xs,
    //                       false>
    //         rng2{2109010154UL, 2891336453UL};

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    {
        pcg64_once_insecure rng1;
        adhoc::PCG_engine<std::uint64_t, 64, adhoc::tempering_type::rxs_m_xs,
                          false>
            rng2{16070087622772795188ULL, 1442695040888963407ULL};

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg128_once_insecure rng1;
        adhoc::PCG_engine<adhoc::uint128, 128, adhoc::tempering_type::xsl_rr_rr,
                          false>
            rng2(245720598905631564143578724636268694099_ULLL,
                 117397592171526113268558934119004209487_ULLL);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            auto val1_low = static_cast<std::uint64_t>(val1);
            auto val2_low = static_cast<std::uint64_t>(val2);
            EXPECT_EQUAL(val1_low, val2_low);
            auto val1_high = static_cast<std::uint64_t>(val1 >> 64U);
            auto val2_high = static_cast<std::uint64_t>(val2 >> 64U);
            EXPECT_EQUAL(val1_high, val2_high);
        }
    }

    // investigate
    // {
    //     pcg32_oneseq_once_insecure rng1;
    //     adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::rxs_m_xs,
    //                       false>
    //         rng2(2109010154UL, 2891336453UL);

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    // {
    //     pcg32_k2 rng1;
    //     adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rr,
    //                       true, 1, 16, false, false, true, true>
    //         rng2(16070087622772795188ULL, 1442695040888963407ULL);

    //     for (std::size_t i = 0; i < sims; ++i) {
    //         auto val1 = rng1();
    //         auto val2 = rng2();
    //         EXPECT_EQUAL(val1, val2);
    //     }
    // }

    {
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rr,
                          true, 1, 16, false, false, true, true>
            rng(16070087622772795188ULL, 1442695040888963407ULL);
        check_fwd_and_back(rng, 1000000);
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rr,
                          true, 1, 16, false, false, true, true>
            rng2(16070087622772795188ULL, 1442695040888963407ULL);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rr,
                          true, 1, 16, false, false, true, true>
            rng(16070087622772795188ULL, 1442695040888963407ULL);
        check_back_and_fwd(rng, 1000000);
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rr,
                          true, 1, 16, false, false, true, true>
            rng2(16070087622772795188ULL, 1442695040888963407ULL);
        EXPECT_EQUAL(rng, rng2);
    }

    {
        pcg32_k2_fast rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rs,
                          true, 1, 32, false, false, true, true>
            rng2(16070087622772795188ULL, 1442695040888963407ULL);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_k64 rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rr,
                          true, 6, 16, false, false, true, true>
            rng2(16070087622772795188ULL, 1442695040888963407ULL);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_k64_oneseq rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rs,
                          true, 6, 32, true, true, true, true>
            rng2(7515307391605212643ULL, 0);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_k64_fast rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rs,
                          true, 6, 32, false, false, true, true>
            rng2(16070087622772795188ULL, 1442695040888963407ULL);

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_c64 rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rr,
                          true, 6, 16, false, false, false, true>
            rng2(16070087622772795188ULL, 1442695040888963407ULL);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_c64_oneseq rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rs,
                          true, 6, 32, false, false, false, true>
            rng2(16070087622772795188ULL, 1442695040888963407ULL);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_c64_fast rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rs,
                          true, 6, 32, false, false, false, true>
            rng2(7515307391605212643ULL, 0);

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
                          true, 10, 16, false, false, true, true>
            rng2(16070087622772795188ULL, 1442695040888963407ULL);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_k1024_fast rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rs,
                          true, 10, 32, false, false, true, true>
            rng2(16070087622772795188ULL, 1442695040888963407ULL);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_c1024 rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rr,
                          true, 10, 16, false, false, false, true>
            rng2(16070087622772795188ULL, 1442695040888963407ULL);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_c1024_fast rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rs,
                          true, 10, 32, false, false, false, true>
            rng2(16070087622772795188ULL, 1442695040888963407ULL);

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
                          true, 14, 16, false, false, true, true>
            rng2(16070087622772795188ULL, 1442695040888963407ULL);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        pcg32_k16384_fast rng1;
        adhoc::PCG_engine<std::uint32_t, 32, adhoc::tempering_type::xsh_rs,
                          true, 14, 32, false, false, true, true>
            rng2(16070087622772795188ULL, 1442695040888963407ULL);

        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = rng1();
            auto val2 = rng2();
            EXPECT_EQUAL(val1, val2);
        }
    }
}
