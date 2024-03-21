#include <partition/integer_partition.hpp>
#include <partition/partition_function.hpp>

#include <gtest/gtest.h>

namespace adhoc3 {

TEST(Partition, IntegerPartition) {

    {
        constexpr std::array<std::array<std::size_t, 1>, partition_function(1)>
            parts{{{1}}};

        constexpr auto first = FirstPartition<1>();
        static_assert(first.arr == parts.front());

        constexpr auto last = LastPartition<1>();
        static_assert(last.arr == parts.back());

        constexpr auto next = NextPartition(first.arr, last.arr);
        static_assert(next.arr == std::array<std::size_t, 1>{});

        constexpr auto prev = PrevPartition(last.arr, first.arr);
        static_assert(prev.arr == std::array<std::size_t, 1>{});
    }

    {
        constexpr std::array<std::array<std::size_t, 2>, partition_function(2)>
            parts{{{0, 1}, {2, 0}}};

        constexpr auto first = FirstPartition<2>();
        static_assert(first.arr == parts.front());

        constexpr auto last = LastPartition<2>();
        static_assert(last.arr == parts.back());

        constexpr auto next = NextPartition(first.arr, last.arr);
        static_assert(next.arr == parts[1]);

        constexpr auto zero = NextPartition(next.arr, last.arr);
        static_assert(zero.arr == std::array<std::size_t, 2>{});

        constexpr auto prev = PrevPartition(last.arr, first.arr);
        static_assert(prev.arr == parts[0]);

        constexpr auto zero2 = PrevPartition(prev.arr, first.arr);
        static_assert(zero2.arr == std::array<std::size_t, 2>{});
    }

    {
        constexpr std::array<std::array<std::size_t, 3>, partition_function(3)>
            parts{{{0, 0, 1}, {1, 1, 0}, {3, 0, 0}}};

        constexpr auto first = FirstPartition<3>();
        static_assert(first.arr == parts.front());

        constexpr auto last = LastPartition<3>();
        static_assert(last.arr == parts.back());

        constexpr auto next1 = NextPartition(first.arr, last.arr);
        static_assert(next1.arr == parts[1]);

        constexpr auto next2 = NextPartition(next1.arr, last.arr);
        static_assert(next2.arr == parts[2]);

        constexpr auto zero = NextPartition(next2.arr, last.arr);
        static_assert(zero.arr == std::array<std::size_t, 3>{});

        constexpr auto prev1 = PrevPartition(last.arr, first.arr);
        static_assert(prev1.arr == parts[1]);

        constexpr auto prev2 = PrevPartition(prev1.arr, first.arr);
        static_assert(prev2.arr == parts[0]);

        constexpr auto zero2 = PrevPartition(prev2.arr, first.arr);
        static_assert(zero2.arr == std::array<std::size_t, 3>{});
    }

    {
        constexpr std::array<std::array<std::size_t, 4>, partition_function(4)>
            parts{{{0, 0, 0, 1},
                   {1, 0, 1, 0},
                   {0, 2, 0, 0},
                   {2, 1, 0, 0},
                   {4, 0, 0, 0}}};

        constexpr auto first = FirstPartition<4>();
        static_assert(first.arr == parts.front());

        constexpr auto last = LastPartition<4>();
        static_assert(last.arr == parts.back());

        constexpr auto next1 = NextPartition(first.arr, last.arr);
        static_assert(next1.arr == parts[1]);

        constexpr auto next2 = NextPartition(next1.arr, last.arr);
        static_assert(next2.arr == parts[2]);

        constexpr auto next3 = NextPartition(next2.arr, last.arr);
        static_assert(next3.arr == parts[3]);

        constexpr auto next4 = NextPartition(next3.arr, last.arr);
        static_assert(next4.arr == parts[4]);

        constexpr auto zero = NextPartition(next4.arr, last.arr);
        static_assert(zero.arr == std::array<std::size_t, 4>{});

        constexpr auto prev1 = PrevPartition(last.arr, first.arr);
        static_assert(prev1.arr == parts[3]);

        constexpr auto prev2 = PrevPartition(prev1.arr, first.arr);
        static_assert(prev2.arr == parts[2]);

        constexpr auto prev3 = PrevPartition(prev2.arr, first.arr);
        static_assert(prev3.arr == parts[1]);

        constexpr auto prev4 = PrevPartition(prev3.arr, first.arr);
        static_assert(prev4.arr == parts[0]);

        constexpr auto zero2 = PrevPartition(prev4.arr, first.arr);
        static_assert(zero2.arr == std::array<std::size_t, 4>{});
    }

    {
        constexpr std::array<std::array<std::size_t, 5>, partition_function(5)>
            parts{{{0, 0, 0, 0, 1},
                   {1, 0, 0, 1, 0},
                   {0, 1, 1, 0, 0},
                   {2, 0, 1, 0, 0},
                   {1, 2, 0, 0, 0},
                   {3, 1, 0, 0, 0},
                   {5, 0, 0, 0, 0}}};

        constexpr auto first = FirstPartition<5>();
        static_assert(first.arr == parts.front());

        constexpr auto last = LastPartition<5>();
        static_assert(last.arr == parts.back());

        constexpr auto next1 = NextPartition(first.arr, last.arr);
        static_assert(next1.arr == parts[1]);

        constexpr auto next2 = NextPartition(next1.arr, last.arr);
        static_assert(next2.arr == parts[2]);

        constexpr auto next3 = NextPartition(next2.arr, last.arr);
        static_assert(next3.arr == parts[3]);

        constexpr auto next4 = NextPartition(next3.arr, last.arr);
        static_assert(next4.arr == parts[4]);

        constexpr auto next5 = NextPartition(next4.arr, last.arr);
        static_assert(next5.arr == parts[5]);

        constexpr auto next6 = NextPartition(next5.arr, last.arr);
        static_assert(next6.arr == parts[6]);

        constexpr auto zero = NextPartition(next6.arr, last.arr);
        static_assert(zero.arr == std::array<std::size_t, 5>{});

        constexpr auto prev1 = PrevPartition(last.arr, first.arr);
        static_assert(prev1.arr == parts[5]);

        constexpr auto prev2 = PrevPartition(prev1.arr, first.arr);
        static_assert(prev2.arr == parts[4]);

        constexpr auto prev3 = PrevPartition(prev2.arr, first.arr);
        static_assert(prev3.arr == parts[3]);

        constexpr auto prev4 = PrevPartition(prev3.arr, first.arr);
        static_assert(prev4.arr == parts[2]);

        constexpr auto prev5 = PrevPartition(prev4.arr, first.arr);
        static_assert(prev5.arr == parts[1]);

        constexpr auto prev6 = PrevPartition(prev5.arr, first.arr);
        static_assert(prev6.arr == parts[0]);

        constexpr auto zero2 = PrevPartition(prev6.arr, first.arr);
        static_assert(zero2.arr == std::array<std::size_t, 5>{});
    }

    {
        constexpr std::array<std::array<std::size_t, 6>, partition_function(6)>
            parts{{{0, 0, 0, 0, 0, 1},
                   {1, 0, 0, 0, 1, 0},
                   {0, 1, 0, 1, 0, 0},
                   {2, 0, 0, 1, 0, 0},
                   {0, 0, 2, 0, 0, 0},
                   {1, 1, 1, 0, 0, 0},
                   {3, 0, 1, 0, 0, 0},
                   {0, 3, 0, 0, 0, 0},
                   {2, 2, 0, 0, 0, 0},
                   {4, 1, 0, 0, 0, 0},
                   {6, 0, 0, 0, 0, 0}}};

        constexpr auto first = FirstPartition<6>();
        static_assert(first.arr == parts.front());

        constexpr auto last = LastPartition<6>();
        static_assert(last.arr == parts.back());

        constexpr auto next1 = NextPartition(first.arr, last.arr);
        static_assert(next1.arr == parts[1]);

        constexpr auto next2 = NextPartition(next1.arr, last.arr);
        static_assert(next2.arr == parts[2]);

        constexpr auto next3 = NextPartition(next2.arr, last.arr);
        static_assert(next3.arr == parts[3]);

        constexpr auto next4 = NextPartition(next3.arr, last.arr);
        static_assert(next4.arr == parts[4]);

        constexpr auto next5 = NextPartition(next4.arr, last.arr);
        static_assert(next5.arr == parts[5]);

        constexpr auto next6 = NextPartition(next5.arr, last.arr);
        static_assert(next6.arr == parts[6]);

        constexpr auto next7 = NextPartition(next6.arr, last.arr);
        static_assert(next7.arr == parts[7]);

        constexpr auto next8 = NextPartition(next7.arr, last.arr);
        static_assert(next8.arr == parts[8]);

        constexpr auto next9 = NextPartition(next8.arr, last.arr);
        static_assert(next9.arr == parts[9]);

        constexpr auto next10 = NextPartition(next9.arr, last.arr);
        static_assert(next10.arr == parts[10]);

        constexpr auto zero = NextPartition(next10.arr, last.arr);
        static_assert(zero.arr == std::array<std::size_t, 6>{});

        constexpr auto prev1 = PrevPartition(last.arr, first.arr);
        static_assert(prev1.arr == parts[9]);

        constexpr auto prev2 = PrevPartition(prev1.arr, first.arr);
        static_assert(prev2.arr == parts[8]);

        constexpr auto prev3 = PrevPartition(prev2.arr, first.arr);
        static_assert(prev3.arr == parts[7]);

        constexpr auto prev4 = PrevPartition(prev3.arr, first.arr);
        static_assert(prev4.arr == parts[6]);

        constexpr auto prev5 = PrevPartition(prev4.arr, first.arr);
        static_assert(prev5.arr == parts[5]);

        constexpr auto prev6 = PrevPartition(prev5.arr, first.arr);
        static_assert(prev6.arr == parts[4]);

        constexpr auto prev7 = PrevPartition(prev6.arr, first.arr);
        static_assert(prev7.arr == parts[3]);

        constexpr auto prev8 = PrevPartition(prev7.arr, first.arr);
        static_assert(prev8.arr == parts[2]);

        constexpr auto prev9 = PrevPartition(prev8.arr, first.arr);
        static_assert(prev9.arr == parts[1]);

        constexpr auto prev10 = PrevPartition(prev9.arr, first.arr);
        static_assert(prev10.arr == parts[0]);

        constexpr auto zero2 = PrevPartition(prev10.arr, first.arr);
        static_assert(zero2.arr == std::array<std::size_t, 6>{});
    }

    {
        constexpr std::array<std::array<std::size_t, 7>, partition_function(7)>
            parts{{{0, 0, 0, 0, 0, 0, 1},
                   {1, 0, 0, 0, 0, 1, 0},
                   {0, 1, 0, 0, 1, 0, 0},
                   {2, 0, 0, 0, 1, 0, 0},
                   {0, 0, 1, 1, 0, 0, 0},
                   {1, 1, 0, 1, 0, 0, 0},
                   {3, 0, 0, 1, 0, 0, 0},
                   {1, 0, 2, 0, 0, 0, 0},
                   {0, 2, 1, 0, 0, 0, 0},
                   {2, 1, 1, 0, 0, 0, 0},
                   {4, 0, 1, 0, 0, 0, 0},
                   {1, 3, 0, 0, 0, 0, 0},
                   {3, 2, 0, 0, 0, 0, 0},
                   {5, 1, 0, 0, 0, 0, 0},
                   {7, 0, 0, 0, 0, 0, 0}}};

        constexpr auto first = FirstPartition<7>();
        static_assert(first.arr == parts.front());

        constexpr auto last = LastPartition<7>();
        static_assert(last.arr == parts.back());

        constexpr auto next1 = NextPartition(first.arr, last.arr);
        static_assert(next1.arr == parts[1]);

        constexpr auto next2 = NextPartition(next1.arr, last.arr);
        static_assert(next2.arr == parts[2]);

        constexpr auto next3 = NextPartition(next2.arr, last.arr);
        static_assert(next3.arr == parts[3]);

        constexpr auto next4 = NextPartition(next3.arr, last.arr);
        static_assert(next4.arr == parts[4]);

        constexpr auto next5 = NextPartition(next4.arr, last.arr);
        static_assert(next5.arr == parts[5]);

        constexpr auto next6 = NextPartition(next5.arr, last.arr);
        static_assert(next6.arr == parts[6]);

        constexpr auto next7 = NextPartition(next6.arr, last.arr);
        static_assert(next7.arr == parts[7]);

        constexpr auto next8 = NextPartition(next7.arr, last.arr);
        static_assert(next8.arr == parts[8]);

        constexpr auto next9 = NextPartition(next8.arr, last.arr);
        static_assert(next9.arr == parts[9]);

        constexpr auto next10 = NextPartition(next9.arr, last.arr);
        static_assert(next10.arr == parts[10]);

        constexpr auto next11 = NextPartition(next10.arr, last.arr);
        static_assert(next11.arr == parts[11]);

        constexpr auto next12 = NextPartition(next11.arr, last.arr);
        static_assert(next12.arr == parts[12]);

        constexpr auto next13 = NextPartition(next12.arr, last.arr);
        static_assert(next13.arr == parts[13]);

        constexpr auto next14 = NextPartition(next13.arr, last.arr);
        static_assert(next14.arr == parts[14]);

        constexpr auto zero = NextPartition(next14.arr, last.arr);
        static_assert(zero.arr == std::array<std::size_t, 7>{});

        constexpr auto prev1 = PrevPartition(last.arr, first.arr);
        static_assert(prev1.arr == parts[13]);

        constexpr auto prev2 = PrevPartition(prev1.arr, first.arr);
        static_assert(prev2.arr == parts[12]);

        constexpr auto prev3 = PrevPartition(prev2.arr, first.arr);
        static_assert(prev3.arr == parts[11]);

        constexpr auto prev4 = PrevPartition(prev3.arr, first.arr);
        static_assert(prev4.arr == parts[10]);

        constexpr auto prev5 = PrevPartition(prev4.arr, first.arr);
        static_assert(prev5.arr == parts[9]);

        constexpr auto prev6 = PrevPartition(prev5.arr, first.arr);
        static_assert(prev6.arr == parts[8]);

        constexpr auto prev7 = PrevPartition(prev6.arr, first.arr);
        static_assert(prev7.arr == parts[7]);

        constexpr auto prev8 = PrevPartition(prev7.arr, first.arr);
        static_assert(prev8.arr == parts[6]);

        constexpr auto prev9 = PrevPartition(prev8.arr, first.arr);
        static_assert(prev9.arr == parts[5]);

        constexpr auto prev10 = PrevPartition(prev9.arr, first.arr);
        static_assert(prev10.arr == parts[4]);

        constexpr auto prev11 = PrevPartition(prev10.arr, first.arr);
        static_assert(prev11.arr == parts[3]);

        constexpr auto prev12 = PrevPartition(prev11.arr, first.arr);
        static_assert(prev12.arr == parts[2]);

        constexpr auto prev13 = PrevPartition(prev12.arr, first.arr);
        static_assert(prev13.arr == parts[1]);

        constexpr auto prev14 = PrevPartition(prev13.arr, first.arr);
        static_assert(prev14.arr == parts[0]);

        constexpr auto zero2 = PrevPartition(prev14.arr, first.arr);
        static_assert(zero2.arr == std::array<std::size_t, 7>{});
    }
}

} // namespace adhoc3
