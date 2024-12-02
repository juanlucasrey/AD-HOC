#include "../../test/test_tools.hpp"
#include "test_tools_rng.hpp"

#include "split_mix.hpp"

#include <vector>

int main() {
    TEST_START;

    // splitmix32 fwd and back
    {
        splitmix32 rng;
        TEST_FUNC(check_fwd_and_back(rng, 100));
    }

    // splitmix32 back and fwd
    {
        splitmix32 rng;
        TEST_FUNC(check_back_and_fwd(rng, 100));
    }

    // splitmix64 fwd and back
    {
        splitmix64 rng;
        TEST_FUNC(check_fwd_and_back(rng, 100));
    }

    // splitmix64 back and fwd
    {
        splitmix64 rng;
        TEST_FUNC(check_back_and_fwd(rng, 100));
    }

    TEST_END;
}
