#include "../../../test_simple/test_simple_include.hpp"

#include <tools/unshift.hpp>

#include <random>

int main() {

    {
        std::mt19937 gen32;
        std::size_t sims = 100;
        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = gen32();
            auto val_s = val1 ^ (val1 << 16);
            auto val2 = adhoc::unshift_left_xor<32, 16>(val_s);
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        std::mt19937 gen32;
        std::size_t sims = 100;
        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = gen32();
            auto val_s = val1 ^ (val1 << 2);
            auto val2 = adhoc::unshift_left_xor<32, 2>(val_s);
            EXPECT_EQUAL(val1, val2);
        }
    }

    {
        std::mt19937 gen32;
        std::size_t sims = 100;
        for (std::size_t i = 0; i < sims; ++i) {
            auto val1 = gen32();
            auto val_s = val1 ^ (val1 << 5);
            auto val2 = adhoc::unshift_left_xor<32, 5>(val_s);
            EXPECT_EQUAL(val1, val2);
        }
    }

    TEST_END;
}
