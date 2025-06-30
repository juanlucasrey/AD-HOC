#include <tools/mask.hpp>

#include <cstdint>

int main() {
    constexpr auto m1 = adhoc::mask<std::uint32_t>(5, 4);
    static_assert(m1 == 496U);

    constexpr auto m2 = adhoc::mask<std::uint32_t>(32);
    static_assert(m2 == 4294967295U);

    constexpr auto m3 = adhoc::mask<std::uint64_t>(64);
    static_assert(m3 == 18446744073709551615U);

    constexpr auto m4 = adhoc::mask<std::uint64_t>(62, 2);
    static_assert(m4 == 18446744073709551612U);
}
