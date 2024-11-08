#include <cmath>

namespace {

constexpr auto half_square_size(std::size_t max) -> std::size_t {
    return (max * (max + 1)) / 2;
}

constexpr auto half_square_index(std::size_t y, std::size_t x, std::size_t max)
    -> std::size_t {
    return x + max * y - (y * (y - 1)) / 2;
}

} // namespace

int main() {
    // EXPERIMENTAL, for operator * values
    static_assert(half_square_size(5) == 15);
    static_assert(half_square_index(0, 0, 5) == 0);
    static_assert(half_square_index(0, 1, 5) == 1);
    static_assert(half_square_index(0, 2, 5) == 2);
    static_assert(half_square_index(0, 3, 5) == 3);
    static_assert(half_square_index(0, 4, 5) == 4);
    static_assert(half_square_index(1, 0, 5) == 5);
    static_assert(half_square_index(1, 1, 5) == 6);
    static_assert(half_square_index(1, 2, 5) == 7);
    static_assert(half_square_index(1, 3, 5) == 8);
    static_assert(half_square_index(2, 0, 5) == 9);
    static_assert(half_square_index(2, 1, 5) == 10);
    static_assert(half_square_index(2, 2, 5) == 11);
    static_assert(half_square_index(3, 0, 5) == 12);
    static_assert(half_square_index(3, 1, 5) == 13);
    static_assert(half_square_index(4, 0, 5) == 14);
    return 0;
}
