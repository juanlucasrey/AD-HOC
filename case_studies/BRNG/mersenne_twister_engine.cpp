#include "mersenne_twister_engine.hpp"

#include <random>
#include <vector>

int main() {
    // check against std
    {
        std::mt19937 mt;
        adhoc::mt19937 mt2;

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

    // check against std 64
    {
        std::mt19937_64 mt;
        adhoc::mt19937_64 mt2;

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

    // 32 fwd and back
    {
        adhoc::mt19937 mt;

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

    // 32 back and fwd
    {
        adhoc::mt19937 mt;

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

    // 64 fwd and back
    {
        adhoc::mt19937_64 mt64;

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

    // 64 back and fwd
    {
        adhoc::mt19937_64 mt64;

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
