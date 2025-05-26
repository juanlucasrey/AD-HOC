#include "../../include/tools/modular_arithmetic.hpp"

#include <cstdint>

int main() {
    constexpr std::uint64_t a13n = 810728U;
    constexpr std::uint64_t a23n = 1370589U;
    constexpr std::uint64_t m1 = 4294967087U;
    constexpr std::uint64_t m2 = 4294944443U;

    constexpr auto a13n_inv = adhoc::modular_multiplicative_inverse(m1, a13n);
    static_assert(a13n_inv == 2349796154U);
    constexpr auto a23n_inv = adhoc::modular_multiplicative_inverse(m2, a23n);
    static_assert(a23n_inv == 69372715U);
}
