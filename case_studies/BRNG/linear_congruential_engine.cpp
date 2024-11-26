#include "linear_congruential_engine.hpp"

#include <random>
#include <vector>

int main() {

    // check against std minstd_rand0
    {
        std::minstd_rand0 mt1;
        adhoc::minstd_rand0 mt2;

        for (std::size_t i = 0; i < 100; ++i) {
            if (mt1() != mt2()) {
                return 1;
            }
        }

        std::uniform_real_distribution<> dis(1.0, 2.0);
        for (std::size_t i = 0; i < 100; ++i) {
            if (dis(mt1) != dis(mt2)) {
                return 1;
            }
        }

        std::uniform_int_distribution<> distrib(1, 6);
        for (std::size_t i = 0; i < 100; ++i) {
            if (distrib(mt1) != distrib(mt2)) {
                return 1;
            }
        }
    }

    // check against std minstd_rand
    {
        std::minstd_rand mt;
        adhoc::minstd_rand mt2;

        for (std::size_t i = 0; i < 100; ++i) {
            if (mt() != mt2()) {
                return 1;
            }
        }

        std::uniform_real_distribution<> dis(1.0, 2.0);
        for (std::size_t i = 0; i < 100; ++i) {
            if (dis(mt) != dis(mt2)) {
                return 1;
            }
        }

        std::uniform_int_distribution<> distrib(1, 6);
        for (std::size_t i = 0; i < 100; ++i) {
            if (distrib(mt) != distrib(mt2)) {
                return 1;
            }
        }
    }

    // minstd_rand0 fwd and back
    {
        adhoc::minstd_rand0 mt;

        std::size_t n = 100;
        std::vector<std::uint_fast32_t> vals(100);

        for (std::size_t i = 0; i < n; ++i) {
            vals[i] = mt();
        }

        std::vector<std::uint_fast32_t> vals2(100);
        for (std::size_t i = n; i--;) {
            vals2[i] = mt.operator()<false>();
        }

        for (std::size_t i = 0; i < n; ++i) {
            if (vals[i] != vals2[i]) {
                return 1;
            }
        }
    }

    // minstd_rand0 back and fwd
    {
        adhoc::minstd_rand0 mt;

        std::size_t n = 100;
        std::vector<std::uint_fast32_t> vals(100);

        for (std::size_t i = 0; i < n; ++i) {
            vals[i] = mt.operator()<false>();
        }

        std::vector<std::uint_fast32_t> vals2(100);
        for (std::size_t i = n; i--;) {
            vals2[i] = mt();
        }

        for (std::size_t i = 0; i < n; ++i) {
            if (vals[i] != vals2[i]) {
                return 1;
            }
        }
    }

    // minstd_rand fwd and back
    {
        adhoc::minstd_rand mt64;

        std::size_t n = 100;
        std::vector<std::uint_fast64_t> vals(100);

        for (std::size_t i = 0; i < n; ++i) {
            vals[i] = mt64();
        }

        std::vector<std::uint_fast64_t> vals2(100);
        for (std::size_t i = n; i--;) {
            vals2[i] = mt64.operator()<false>();
        }

        for (std::size_t i = 0; i < n; ++i) {
            if (vals[i] != vals2[i]) {
                return 1;
            }
        }
    }

    // minstd_rand back and fwd
    {
        adhoc::minstd_rand mt64;

        std::size_t n = 100;
        std::vector<std::uint_fast64_t> vals(100);

        for (std::size_t i = 0; i < n; ++i) {
            vals[i] = mt64.operator()<false>();
        }

        std::vector<std::uint_fast64_t> vals2(100);
        for (std::size_t i = n; i--;) {
            vals2[i] = mt64();
        }

        for (std::size_t i = 0; i < n; ++i) {
            if (vals[i] != vals2[i]) {
                return 1;
            }
        }
    }

    return 0;
}
