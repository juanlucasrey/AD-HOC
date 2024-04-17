#include <partition/integer_partition.hpp>
#include <partition/partition_function.hpp>

#include <gtest/gtest.h>

namespace adhoc3 {

TEST(Partition, IntegerPartition) {

    {
        constexpr std::array<std::array<std::size_t, 1>, partition_function(1)>
            parts{{{1}}};

        constexpr auto first = FirstPartition<1>();
        static_assert(first == parts.front());

        constexpr auto last = LastPartition<1>();
        static_assert(last == parts.back());

        constexpr auto next = NextPartition(first);
        static_assert(next == last);

        constexpr auto prev = PrevPartition(last);
        static_assert(prev == last);
    }

    {
        constexpr std::array<std::array<std::size_t, 2>, partition_function(2)>
            parts{{{0, 1}, {2, 0}}};

        constexpr auto first = FirstPartition<2>();
        static_assert(first == parts.front());

        constexpr auto last = LastPartition<2>();
        static_assert(last == parts.back());

        constexpr auto next = NextPartition(first);
        static_assert(next == parts[1]);

        constexpr auto zero = NextPartition(next);
        static_assert(zero == last);

        constexpr auto prev = PrevPartition(last);
        static_assert(prev == parts[0]);

        constexpr auto zero2 = PrevPartition(prev);
        static_assert(zero2 == first);
    }

    {
        constexpr std::array<std::array<std::size_t, 3>, partition_function(3)>
            parts{{{0, 0, 1}, {1, 1, 0}, {3, 0, 0}}};

        constexpr auto first = FirstPartition<3>();
        static_assert(first == parts.front());

        constexpr auto last = LastPartition<3>();
        static_assert(last == parts.back());

        constexpr auto next1 = NextPartition(first);
        static_assert(next1 == parts[1]);

        constexpr auto next2 = NextPartition(next1);
        static_assert(next2 == parts[2]);

        constexpr auto zero = NextPartition(next2);
        static_assert(zero == last);

        constexpr auto prev1 = PrevPartition(last);
        static_assert(prev1 == parts[1]);

        constexpr auto prev2 = PrevPartition(prev1);
        static_assert(prev2 == parts[0]);

        constexpr auto zero2 = PrevPartition(prev2);
        static_assert(zero2 == first);
    }

    {
        constexpr std::array<std::array<std::size_t, 4>, partition_function(4)>
            parts{{{0, 0, 0, 1},
                   {1, 0, 1, 0},
                   {0, 2, 0, 0},
                   {2, 1, 0, 0},
                   {4, 0, 0, 0}}};

        constexpr auto first = FirstPartition<4>();
        static_assert(first == parts.front());

        constexpr auto last = LastPartition<4>();
        static_assert(last == parts.back());

        constexpr auto next1 = NextPartition(first);
        static_assert(next1 == parts[1]);

        constexpr auto next2 = NextPartition(next1);
        static_assert(next2 == parts[2]);

        constexpr auto next3 = NextPartition(next2);
        static_assert(next3 == parts[3]);

        constexpr auto next4 = NextPartition(next3);
        static_assert(next4 == parts[4]);

        constexpr auto zero = NextPartition(next4);
        static_assert(zero == last);

        constexpr auto prev1 = PrevPartition(last);
        static_assert(prev1 == parts[3]);

        constexpr auto prev2 = PrevPartition(prev1);
        static_assert(prev2 == parts[2]);

        constexpr auto prev3 = PrevPartition(prev2);
        static_assert(prev3 == parts[1]);

        constexpr auto prev4 = PrevPartition(prev3);
        static_assert(prev4 == parts[0]);

        constexpr auto zero2 = PrevPartition(prev4);
        static_assert(zero2 == first);
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
        static_assert(first == parts.front());

        constexpr auto last = LastPartition<5>();
        static_assert(last == parts.back());

        constexpr auto next1 = NextPartition(first);
        static_assert(next1 == parts[1]);

        constexpr auto next2 = NextPartition(next1);
        static_assert(next2 == parts[2]);

        constexpr auto next3 = NextPartition(next2);
        static_assert(next3 == parts[3]);

        constexpr auto next4 = NextPartition(next3);
        static_assert(next4 == parts[4]);

        constexpr auto next5 = NextPartition(next4);
        static_assert(next5 == parts[5]);

        constexpr auto next6 = NextPartition(next5);
        static_assert(next6 == parts[6]);

        constexpr auto zero = NextPartition(next6);
        static_assert(zero == last);

        constexpr auto prev1 = PrevPartition(last);
        static_assert(prev1 == parts[5]);

        constexpr auto prev2 = PrevPartition(prev1);
        static_assert(prev2 == parts[4]);

        constexpr auto prev3 = PrevPartition(prev2);
        static_assert(prev3 == parts[3]);

        constexpr auto prev4 = PrevPartition(prev3);
        static_assert(prev4 == parts[2]);

        constexpr auto prev5 = PrevPartition(prev4);
        static_assert(prev5 == parts[1]);

        constexpr auto prev6 = PrevPartition(prev5);
        static_assert(prev6 == parts[0]);

        constexpr auto zero2 = PrevPartition(prev6);
        static_assert(zero2 == first);
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
        static_assert(first == parts.front());

        constexpr auto last = LastPartition<6>();
        static_assert(last == parts.back());

        constexpr auto next1 = NextPartition(first);
        static_assert(next1 == parts[1]);

        constexpr auto next2 = NextPartition(next1);
        static_assert(next2 == parts[2]);

        constexpr auto next3 = NextPartition(next2);
        static_assert(next3 == parts[3]);

        constexpr auto next4 = NextPartition(next3);
        static_assert(next4 == parts[4]);

        constexpr auto next5 = NextPartition(next4);
        static_assert(next5 == parts[5]);

        constexpr auto next6 = NextPartition(next5);
        static_assert(next6 == parts[6]);

        constexpr auto next7 = NextPartition(next6);
        static_assert(next7 == parts[7]);

        constexpr auto next8 = NextPartition(next7);
        static_assert(next8 == parts[8]);

        constexpr auto next9 = NextPartition(next8);
        static_assert(next9 == parts[9]);

        constexpr auto next10 = NextPartition(next9);
        static_assert(next10 == parts[10]);

        constexpr auto zero = NextPartition(next10);
        static_assert(zero == last);

        constexpr auto prev1 = PrevPartition(last);
        static_assert(prev1 == parts[9]);

        constexpr auto prev2 = PrevPartition(prev1);
        static_assert(prev2 == parts[8]);

        constexpr auto prev3 = PrevPartition(prev2);
        static_assert(prev3 == parts[7]);

        constexpr auto prev4 = PrevPartition(prev3);
        static_assert(prev4 == parts[6]);

        constexpr auto prev5 = PrevPartition(prev4);
        static_assert(prev5 == parts[5]);

        constexpr auto prev6 = PrevPartition(prev5);
        static_assert(prev6 == parts[4]);

        constexpr auto prev7 = PrevPartition(prev6);
        static_assert(prev7 == parts[3]);

        constexpr auto prev8 = PrevPartition(prev7);
        static_assert(prev8 == parts[2]);

        constexpr auto prev9 = PrevPartition(prev8);
        static_assert(prev9 == parts[1]);

        constexpr auto prev10 = PrevPartition(prev9);
        static_assert(prev10 == parts[0]);

        constexpr auto zero2 = PrevPartition(prev10);
        static_assert(zero2 == first);
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
        static_assert(first == parts.front());

        constexpr auto last = LastPartition<7>();
        static_assert(last == parts.back());

        constexpr auto next1 = NextPartition(first);
        static_assert(next1 == parts[1]);

        constexpr auto next2 = NextPartition(next1);
        static_assert(next2 == parts[2]);

        constexpr auto next3 = NextPartition(next2);
        static_assert(next3 == parts[3]);

        constexpr auto next4 = NextPartition(next3);
        static_assert(next4 == parts[4]);

        constexpr auto next5 = NextPartition(next4);
        static_assert(next5 == parts[5]);

        constexpr auto next6 = NextPartition(next5);
        static_assert(next6 == parts[6]);

        constexpr auto next7 = NextPartition(next6);
        static_assert(next7 == parts[7]);

        constexpr auto next8 = NextPartition(next7);
        static_assert(next8 == parts[8]);

        constexpr auto next9 = NextPartition(next8);
        static_assert(next9 == parts[9]);

        constexpr auto next10 = NextPartition(next9);
        static_assert(next10 == parts[10]);

        constexpr auto next11 = NextPartition(next10);
        static_assert(next11 == parts[11]);

        constexpr auto next12 = NextPartition(next11);
        static_assert(next12 == parts[12]);

        constexpr auto next13 = NextPartition(next12);
        static_assert(next13 == parts[13]);

        constexpr auto next14 = NextPartition(next13);
        static_assert(next14 == parts[14]);

        constexpr auto zero = NextPartition(next14);
        static_assert(zero == last);

        constexpr auto prev1 = PrevPartition(last);
        static_assert(prev1 == parts[13]);

        constexpr auto prev2 = PrevPartition(prev1);
        static_assert(prev2 == parts[12]);

        constexpr auto prev3 = PrevPartition(prev2);
        static_assert(prev3 == parts[11]);

        constexpr auto prev4 = PrevPartition(prev3);
        static_assert(prev4 == parts[10]);

        constexpr auto prev5 = PrevPartition(prev4);
        static_assert(prev5 == parts[9]);

        constexpr auto prev6 = PrevPartition(prev5);
        static_assert(prev6 == parts[8]);

        constexpr auto prev7 = PrevPartition(prev6);
        static_assert(prev7 == parts[7]);

        constexpr auto prev8 = PrevPartition(prev7);
        static_assert(prev8 == parts[6]);

        constexpr auto prev9 = PrevPartition(prev8);
        static_assert(prev9 == parts[5]);

        constexpr auto prev10 = PrevPartition(prev9);
        static_assert(prev10 == parts[4]);

        constexpr auto prev11 = PrevPartition(prev10);
        static_assert(prev11 == parts[3]);

        constexpr auto prev12 = PrevPartition(prev11);
        static_assert(prev12 == parts[2]);

        constexpr auto prev13 = PrevPartition(prev12);
        static_assert(prev13 == parts[1]);

        constexpr auto prev14 = PrevPartition(prev13);
        static_assert(prev14 == parts[0]);

        constexpr auto zero2 = PrevPartition(prev14);
        static_assert(zero2 == first);
    }
}

TEST(Partition, BellCoeffs) {
    // https://en.wikipedia.org/wiki/Bell_polynomials coefficients
    // order 1
    {
        static_assert(partition_function(1) == 1);
        constexpr std::array<std::array<std::size_t, 1>, 1> parts{{{1}}};

        static_assert(BellCoeff(parts[0]) == 1);
    }

    // order 2
    {
        static_assert(partition_function(2) == 2);
        constexpr std::array<std::array<std::size_t, 2>, 2> parts{
            {{0, 1}, {2, 0}}};

        static_assert(BellCoeff(parts[0]) == 1);
        static_assert(BellCoeff(parts[1]) == 1);
    }

    // order 3
    {
        static_assert(partition_function(3) == 3);
        constexpr std::array<std::array<std::size_t, 3>, 3> parts{
            {{0, 0, 1}, {1, 1, 0}, {3, 0, 0}}};

        static_assert(BellCoeff(parts[0]) == 1);
        static_assert(BellCoeff(parts[1]) == 3);
        static_assert(BellCoeff(parts[2]) == 1);
    }

    // order 4
    {
        static_assert(partition_function(4) == 5);
        constexpr std::array<std::array<std::size_t, 4>, 5> parts{
            {{0, 0, 0, 1},
             {1, 0, 1, 0},
             {0, 2, 0, 0},
             {2, 1, 0, 0},
             {4, 0, 0, 0}}};

        static_assert(BellCoeff(parts[0]) == 1);
        static_assert(BellCoeff(parts[1]) == 4);
        static_assert(BellCoeff(parts[2]) == 3);
        static_assert(BellCoeff(parts[3]) == 6);
        static_assert(BellCoeff(parts[4]) == 1);
    }

    // order 5
    {
        static_assert(partition_function(5) == 7);
        constexpr std::array<std::array<std::size_t, 5>, 7> parts{
            {{0, 0, 0, 0, 1},
             {1, 0, 0, 1, 0},
             {0, 1, 1, 0, 0},
             {2, 0, 1, 0, 0},
             {1, 2, 0, 0, 0},
             {3, 1, 0, 0, 0},
             {5, 0, 0, 0, 0}}};

        static_assert(BellCoeff(parts[0]) == 1);
        static_assert(BellCoeff(parts[1]) == 5);
        static_assert(BellCoeff(parts[2]) == 10);
        static_assert(BellCoeff(parts[3]) == 10);
        static_assert(BellCoeff(parts[4]) == 15);
        static_assert(BellCoeff(parts[5]) == 10);
        static_assert(BellCoeff(parts[6]) == 1);
    }

    // order 6
    {
        static_assert(partition_function(6) == 11);
        constexpr std::array<std::array<std::size_t, 6>, 11> parts{
            {{0, 0, 0, 0, 0, 1},
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

        static_assert(BellCoeff(parts[0]) == 1);
        static_assert(BellCoeff(parts[1]) == 6);
        static_assert(BellCoeff(parts[2]) == 15);
        static_assert(BellCoeff(parts[3]) == 15);
        static_assert(BellCoeff(parts[4]) == 10);
        static_assert(BellCoeff(parts[5]) == 60);
        static_assert(BellCoeff(parts[6]) == 20);
        static_assert(BellCoeff(parts[7]) == 15);
        static_assert(BellCoeff(parts[8]) == 45);
        static_assert(BellCoeff(parts[9]) == 15);
        static_assert(BellCoeff(parts[10]) == 1);
    }

    // order 7
    {
        static_assert(partition_function(7) == 15);
        constexpr std::array<std::array<std::size_t, 7>, 15> parts{
            {{0, 0, 0, 0, 0, 0, 1},
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

        static_assert(BellCoeff(parts[0]) == 1);
        static_assert(BellCoeff(parts[1]) == 7);
        static_assert(BellCoeff(parts[2]) == 21);
        static_assert(BellCoeff(parts[3]) == 21);
        static_assert(BellCoeff(parts[4]) == 35);
        static_assert(BellCoeff(parts[5]) == 105);
        static_assert(BellCoeff(parts[6]) == 35);
        static_assert(BellCoeff(parts[7]) == 70);
        static_assert(BellCoeff(parts[8]) == 105);
        static_assert(BellCoeff(parts[9]) == 210);
        static_assert(BellCoeff(parts[10]) == 35);
        static_assert(BellCoeff(parts[11]) == 105);
        static_assert(BellCoeff(parts[12]) == 105);
        static_assert(BellCoeff(parts[13]) == 21);
        static_assert(BellCoeff(parts[14]) == 1);
    }
}

} // namespace adhoc3
