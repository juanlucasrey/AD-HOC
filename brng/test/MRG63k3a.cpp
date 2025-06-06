#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"
#include <iostream>

extern "C" {
#include "external/MRG63k3a.h"
}

#include "seed_seq.hpp"
#include <MRG63k3a.hpp>
#include <distribution/uniform_distribution.hpp>

#include <cstdint>

auto main() -> int {
    // this case is designed to have the max value and the uniform open
    // distribution does not touch 1.0
    {
        adhoc::seed_seq<std::uint_fast64_t> seq;
        seq.vals = {8612413762542751856ULL, 7482105212578211481ULL,
                    7509699444841106451ULL, 5670463998974203631ULL,
                    3539266752878466645ULL, 229387088646599462ULL};
        adhoc::open<adhoc::MRG63k3a<std::uint64_t, true, true>> rng(seq);
        for (std::size_t i = 0; i < 15; ++i) {
            auto val1 = rng();
            EXPECT_NOT_EQUAL(val1, 0.);
            EXPECT_NOT_EQUAL(val1, 1.);
        }
    }

    {
        // p2 == 0 case, using uint128
        adhoc::seed_seq<std::uint_fast64_t> seq;
        seq.vals = {8612413762542751856ULL, 7482105212578211481ULL,
                    7509699444841106451ULL, 2679205843562998531,
                    7615793543397192295,    2157354250027567738};
        adhoc::open<adhoc::MRG63k3a<std::uint64_t, true, true>> rng(seq);

        std::vector<long long> init = {
            8612413762542751856ULL, 7482105212578211481ULL,
            7509699444841106451ULL, 2679205843562998531,
            7615793543397192295,    2157354250027567738};
        InitMRG63k3a(init.data());

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = MRG63k3a();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        // p1 == 0 case, using uint128
        adhoc::seed_seq<std::uint_fast64_t> seq;
        seq.vals = {4206931007091270405,    6204676270937440418,
                    1069855422653178874,    5670463998974203631ULL,
                    3539266752878466645ULL, 229387088646599462ULL};
        adhoc::open<adhoc::MRG63k3a<std::uint64_t, true, true>> rng(seq);

        std::vector<long long> init = {
            4206931007091270405,    6204676270937440418,
            1069855422653178874,    5670463998974203631ULL,
            3539266752878466645ULL, 229387088646599462ULL};
        InitMRG63k3a(init.data());

        bool hasone = false;
        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = MRG63k3a();

            hasone |= (val2 == 1.0);
            if (val2 != 1.0) {
                EXPECT_EQUAL(val1, val2);
            }
        }

        // original MRG63k3a returns 1!!! should not happen
        EXPECT_EQUAL(hasone, true);
    }

    {
        // p1 == p2 case, using uint128
        adhoc::seed_seq<std::uint_fast64_t> seq;
        seq.vals = {8612413762542751856ULL, 7482105212578211481ULL,
                    7509699444841106451ULL, 5670463998974203631ULL,
                    3539266752878466645ULL, 229387088646599462ULL};
        adhoc::open<adhoc::MRG63k3a<std::uint64_t, true, true>> rng(seq);

        std::vector<long long> init = {
            8612413762542751856ULL, 7482105212578211481ULL,
            7509699444841106451ULL, 5670463998974203631ULL,
            3539266752878466645ULL, 229387088646599462ULL};
        InitMRG63k3a(init.data());

        bool hasone = false;
        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = MRG63k3a();

            hasone |= (val2 == 1.0);
            if (val2 != 1.0) {
                EXPECT_EQUAL(val1, val2);
            }
        }

        // original MRG63k3a returns 1!!! should not happen
        EXPECT_EQUAL(hasone, true);
    }

    {
        // p2 == 0 case, not using uint128
        adhoc::seed_seq<std::uint_fast64_t> seq;
        seq.vals = {8612413762542751856ULL, 7482105212578211481ULL,
                    7509699444841106451ULL, 2679205843562998531,
                    7615793543397192295,    2157354250027567738};
        adhoc::open<adhoc::MRG63k3a<std::uint64_t, false, true>> rng(seq);

        std::vector<long long> init = {
            8612413762542751856ULL, 7482105212578211481ULL,
            7509699444841106451ULL, 2679205843562998531,
            7615793543397192295,    2157354250027567738};
        InitMRG63k3a(init.data());

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = MRG63k3a();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        // p1 == 0 case, not using uint128
        adhoc::seed_seq<std::uint_fast64_t> seq;
        seq.vals = {4206931007091270405,    6204676270937440418,
                    1069855422653178874,    5670463998974203631ULL,
                    3539266752878466645ULL, 229387088646599462ULL};
        adhoc::open<adhoc::MRG63k3a<std::uint64_t, false, true>> rng(seq);

        std::vector<long long> init = {
            4206931007091270405,    6204676270937440418,
            1069855422653178874,    5670463998974203631ULL,
            3539266752878466645ULL, 229387088646599462ULL};
        InitMRG63k3a(init.data());

        bool hasone = false;
        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = MRG63k3a();
            hasone |= (val2 == 1.0);
            if (val2 != 1.0) {
                EXPECT_EQUAL(val1, val2);
            }
        }

        // original MRG63k3a returns 1!!! should not happen
        EXPECT_EQUAL(hasone, true);
    }

    {
        // p1 == p2 case, not using uint128
        adhoc::seed_seq<std::uint_fast64_t> seq;
        seq.vals = {8612413762542751856ULL, 7482105212578211481ULL,
                    7509699444841106451ULL, 5670463998974203631ULL,
                    3539266752878466645ULL, 229387088646599462ULL};
        adhoc::open<adhoc::MRG63k3a<std::uint64_t, false, true>> rng(seq);

        std::vector<long long> init = {
            8612413762542751856ULL, 7482105212578211481ULL,
            7509699444841106451ULL, 5670463998974203631ULL,
            3539266752878466645ULL, 229387088646599462ULL};
        InitMRG63k3a(init.data());

        bool hasone = false;
        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = MRG63k3a();
            hasone |= (val2 == 1.0);
            if (val2 != 1.0) {
                EXPECT_EQUAL(val1, val2);
            }
        }

        // original MRG63k3a returns 1!!! should not happen
        EXPECT_EQUAL(hasone, true);
    }

    // default not using 128
    {
        adhoc::open<adhoc::MRG63k3a<std::uint64_t, false, true>> rng;
        std::vector<long long> init = {123456789, 123456789, 123456789,
                                       123456789, 123456789, 123456789};
        InitMRG63k3a(init.data());

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = MRG63k3a();
            EXPECT_EQUAL(val1, val2);
        }
    }

    // default  using 128
    {
        adhoc::open<adhoc::MRG63k3a<std::uint64_t, true, true>> rng;
        std::vector<long long> init = {123456789, 123456789, 123456789,
                                       123456789, 123456789, 123456789};
        InitMRG63k3a(init.data());

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = MRG63k3a();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::MRG63k3a<std::uint64_t> rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::MRG63k3a<std::uint64_t> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::MRG63k3a<std::uint64_t> rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::MRG63k3a<std::uint64_t> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

#ifndef _MSC_VER
    {
        adhoc::MRG63k3a<std::uint64_t> rng1;
        adhoc::MRG63k3a<__uint128_t> rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }
#endif

    TEST_END;
}
