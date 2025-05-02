#ifndef CASE_STUDIES_BRNG_BIT_FIND
#define CASE_STUDIES_BRNG_BIT_FIND

#include <array>

namespace adhoc {

namespace detail {
inline auto keepHighestBit(unsigned int n) -> unsigned int {
    n |= (n >> 1);
    n |= (n >> 2);
    n |= (n >> 4);
    n |= (n >> 8);
    n |= (n >> 16);
    return n - (n >> 1);
}
} // namespace detail

inline auto highestBitIndex(unsigned int b) -> std::size_t {
    constexpr unsigned int deBruijnMagic = 0x06EB14F9;
    constexpr std::array<std::size_t, 32> deBruijnTable = {
        0,  1,  16, 2,  29, 17, 3,  22, 30, 20, 18, 11, 13, 4, 7,  23,
        31, 15, 28, 21, 19, 10, 12, 6,  14, 27, 9,  5,  26, 8, 25, 24,
    };
    return deBruijnTable[(detail::keepHighestBit(b) * deBruijnMagic) >> 27];
}

namespace detail {
template <class UIntType>
inline auto findLeastSignificantZeroNaive(UIntType n) -> std::size_t {
    std::size_t c = 0;
    while ((n & 1)) {
        n >>= 1;
        ++c;
    }
    return c;
}

// Returns the position of the least significant 0 bit (0-based index from
// right) Returns 32 if no 0 bits exist (i.e., input is all 1s)
inline auto findLeastSignificantZero32(uint32_t n) -> std::size_t {
    // Invert bits: 1s become 0s and 0s become 1s
    n = ~n;

    constexpr std::array<std::size_t, 32> deBruijnTable = {
        0,  1,  28, 2,  29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4,  8,
        31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6,  11, 5,  10, 9};

    uint32_t isolated = n & -n;
    // Prevent multiplication by 0 and ensure valid table index
    uint32_t adjusted = isolated | static_cast<uint64_t>(isolated == 0);
    uint32_t index = ((adjusted * 0x077CB531U) >> 27U);
    // Use bitwise ops to select between table value and 32
    auto is_zero = static_cast<std::size_t>(isolated == 0);
    return (deBruijnTable[index] & ~(-is_zero)) | (32U & (-is_zero));
}

// Returns the position of the least significant 0 bit (0-based index from
// right) Returns 64 if no 0 bits exist (i.e., input is all 1s)
inline auto findLeastSignificantZero64(uint64_t n) -> std::size_t {
    // Invert bits: 1s become 0s and 0s become 1s
    n = ~n;

    constexpr std::array<std::size_t, 64> deBruijnTable = {
        0,  1,  2,  7,  3,  13, 8,  19, 4,  25, 14, 28, 9,  34, 20, 40,
        5,  17, 26, 38, 15, 46, 29, 48, 10, 31, 35, 54, 21, 50, 41, 57,
        63, 6,  12, 18, 24, 27, 33, 39, 16, 37, 45, 47, 30, 53, 49, 56,
        62, 11, 23, 32, 36, 44, 52, 55, 61, 22, 43, 51, 60, 42, 59, 58};

    uint64_t isolated = n & -n;
    // Prevent multiplication by 0 and ensure valid table index
    uint64_t adjusted = isolated | static_cast<uint64_t>(isolated == 0);
    uint64_t index = ((adjusted * 0x218A392CD3D5DBFULL) >> 58U);
    // Use bitwise ops to select between table value and 64
    auto is_zero = static_cast<std::size_t>(isolated == 0);
    return (deBruijnTable[index] & ~(-is_zero)) | (64U & (-is_zero));
}
} // namespace detail

template <bool Fast = true, class UIntType>
inline auto findLeastSignificantZero(UIntType n) -> std::size_t {
    if constexpr (Fast) {
        if constexpr (std::numeric_limits<UIntType>::digits == 32) {
            return detail::findLeastSignificantZero32(n);
        } else {
            static_assert(std::numeric_limits<UIntType>::digits == 64);
            return detail::findLeastSignificantZero64(n);
        }
    } else {
        return detail::findLeastSignificantZeroNaive(n);
    }
}

} // namespace adhoc

#endif // CASE_STUDIES_BRNG_BIT_FIND
