#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

extern "C" {
#include "external/MRG31k3p.h"
}

#include <MRG31k3p.hpp>
#include <distribution/uniform_distribution.hpp>

#include <cstdint>

auto main() -> int {
    {
        adhoc::canonical<adhoc::MRG31k3p<std::uint32_t>> rng;

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = MRG31k3p();
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        adhoc::MRG31k3p<std::uint32_t> rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::MRG31k3p<std::uint32_t> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::MRG31k3p<std::uint32_t> rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::MRG31k3p<std::uint32_t> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::MRG31k3p<std::uint64_t> rng;
        check_fwd_and_back(rng, 1000000);
        adhoc::MRG31k3p<std::uint64_t> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::MRG31k3p<std::uint64_t> rng;
        check_back_and_fwd(rng, 1000000);
        adhoc::MRG31k3p<std::uint64_t> rng2;
        EXPECT_EQUAL(rng, rng2);
    }

    {
        adhoc::MRG31k3p<std::uint32_t> rng1;
        adhoc::MRG31k3p<std::uint64_t> rng2;
        compare_rng(rng1, rng2, 1000000);
        compare_rng_limits(rng1, rng2);
    }

    TEST_END;
}
