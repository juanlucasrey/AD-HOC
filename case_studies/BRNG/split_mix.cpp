#include "split_mix.hpp"

#include <vector>

int main() {

    // minstd_rand0 fwd and back
    {
        splitmix32 rng;

        std::size_t n = 100;
        std::vector<std::uint32_t> vals(100);

        for (std::size_t i = 0; i < n; ++i) {
            vals[i] = rng();
        }

        std::vector<std::uint32_t> vals2(100);
        for (std::size_t i = n; i--;) {
            vals2[i] = rng.operator()<false>();
        }

        for (std::size_t i = 0; i < n; ++i) {
            if (vals[i] != vals2[i]) {
                return 1;
            }
        }
    }

    // minstd_rand0 back and fwd
    {
        splitmix32 rng;

        std::size_t n = 100;
        std::vector<std::uint32_t> vals(100);

        for (std::size_t i = 0; i < n; ++i) {
            vals[i] = rng.operator()<false>();
        }

        std::vector<std::uint32_t> vals2(100);
        for (std::size_t i = n; i--;) {
            vals2[i] = rng();
        }

        for (std::size_t i = 0; i < n; ++i) {
            if (vals[i] != vals2[i]) {
                return 1;
            }
        }
    }

    // minstd_rand fwd and back
    {
        splitmix64 rng;

        std::size_t n = 100;
        std::vector<std::uint64_t> vals(100);

        for (std::size_t i = 0; i < n; ++i) {
            vals[i] = rng();
        }

        std::vector<std::uint64_t> vals2(100);
        for (std::size_t i = n; i--;) {
            vals2[i] = rng.operator()<false>();
        }

        for (std::size_t i = 0; i < n; ++i) {
            if (vals[i] != vals2[i]) {
                return 1;
            }
        }
    }

    // minstd_rand back and fwd
    {
        splitmix64 rng;

        std::size_t n = 100;
        std::vector<std::uint64_t> vals(100);

        for (std::size_t i = 0; i < n; ++i) {
            vals[i] = rng.operator()<false>();
        }

        std::vector<std::uint64_t> vals2(100);
        for (std::size_t i = n; i--;) {
            vals2[i] = rng();
        }

        for (std::size_t i = 0; i < n; ++i) {
            if (vals[i] != vals2[i]) {
                return 1;
            }
        }
    }

    return 0;
}
