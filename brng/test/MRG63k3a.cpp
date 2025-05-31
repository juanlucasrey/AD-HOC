#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

extern "C" {
#include "external/MRG63k3a.h"
}

#include <MRG63k3a.hpp>
#include <distribution/uniform_distribution.hpp>

#include <cstdint>

auto main() -> int {

    {
        adhoc::canonical<adhoc::MRG63k3a<std::uint64_t>> rng;

        for (std::size_t i = 0; i < 1000000; ++i) {
            auto val1 = rng();
            auto val2 = MRG63k3a();
            EXPECT_EQUAL(val1, val2);
            if (val1 != val2) {
                std::cout << val1 << std::endl;
                std::cout << val2 << std::endl;
                std::cout << i << std::endl;
            }
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
