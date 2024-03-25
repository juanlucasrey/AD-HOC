#include <partition/combinations.hpp>
#include <partition/multinomial_coefficient.hpp>

#include <gtest/gtest.h>

namespace adhoc3 {

TEST(Partition, Multinomial) {

    {
        constexpr std::array<std::array<std::size_t, 1>, combinations(1, 1)>
            parts{{{1}}};

        constexpr auto first = FirstMultinomial<1, 1>();
        static_assert(first == parts.front());

        constexpr auto last = LastMultinomial<1, 1>();
        static_assert(last == parts.back());

        constexpr auto last2 = NextMultinomial(last);
        static_assert(last == last2);
    }

    {
        constexpr std::array<std::array<std::size_t, 1>, combinations(1, 2)>
            parts{{{2}}};

        constexpr auto first = FirstMultinomial<1, 2>();
        static_assert(first == parts.front());

        constexpr auto last = LastMultinomial<1, 2>();
        static_assert(last == parts.back());

        constexpr auto last2 = NextMultinomial(last);
        static_assert(last == last2);
    }

    {
        constexpr std::array<std::array<std::size_t, 1>, combinations(1, 3)>
            parts{{{3}}};

        constexpr auto first = FirstMultinomial<1, 3>();
        static_assert(first == parts.front());

        constexpr auto last = LastMultinomial<1, 3>();
        static_assert(last == parts.back());

        constexpr auto last2 = NextMultinomial(last);
        static_assert(last == last2);
    }

    {
        constexpr std::array<std::array<std::size_t, 2>, combinations(2, 1)>
            parts{{{1, 0}, {0, 1}}};

        constexpr auto first = FirstMultinomial<2, 1>();
        static_assert(first == parts.front());

        constexpr auto next1 = NextMultinomial(first);
        static_assert(next1 == parts[1]);

        constexpr auto last = LastMultinomial<2, 1>();
        static_assert(last == parts.back());

        constexpr auto last2 = NextMultinomial(last);
        static_assert(last == last2);
    }

    {
        constexpr std::array<std::array<std::size_t, 2>, combinations(2, 2)>
            parts{{{2, 0}, {1, 1}, {0, 2}}};

        constexpr auto first = FirstMultinomial<2, 2>();
        static_assert(first == parts.front());

        constexpr auto next1 = NextMultinomial(first);
        static_assert(next1 == parts[1]);

        constexpr auto next2 = NextMultinomial(next1);
        static_assert(next2 == parts[2]);

        constexpr auto last = LastMultinomial<2, 2>();
        static_assert(last == parts.back());
    }

    {
        constexpr std::array<std::array<std::size_t, 2>, combinations(2, 3)>
            parts{{{3, 0}, {2, 1}, {1, 2}, {0, 3}}};

        constexpr auto first = FirstMultinomial<2, 3>();
        static_assert(first == parts.front());

        constexpr auto next1 = NextMultinomial(first);
        static_assert(next1 == parts[1]);

        constexpr auto next2 = NextMultinomial(next1);
        static_assert(next2 == parts[2]);

        constexpr auto next3 = NextMultinomial(next2);
        static_assert(next3 == parts[3]);

        constexpr auto last = LastMultinomial<2, 3>();
        static_assert(last == parts.back());
    }

    {
        constexpr std::array<std::array<std::size_t, 4>, combinations(4, 3)>
            parts{{{3, 0, 0, 0}, {2, 1, 0, 0}, {2, 0, 1, 0}, {2, 0, 0, 1},
                   {1, 2, 0, 0}, {1, 1, 1, 0}, {1, 1, 0, 1}, {1, 0, 2, 0},
                   {1, 0, 1, 1}, {1, 0, 0, 2}, {0, 3, 0, 0}, {0, 2, 1, 0},
                   {0, 2, 0, 1}, {0, 1, 2, 0}, {0, 1, 1, 1}, {0, 1, 0, 2},
                   {0, 0, 3, 0}, {0, 0, 2, 1}, {0, 0, 1, 2}, {0, 0, 0, 3}}};

        constexpr auto first = FirstMultinomial<4, 3>();
        static_assert(first == parts.front());

        constexpr auto next1 = NextMultinomial(first);
        static_assert(next1 == parts[1]);

        constexpr auto next2 = NextMultinomial(next1);
        static_assert(next2 == parts[2]);

        constexpr auto next3 = NextMultinomial(next2);
        static_assert(next3 == parts[3]);

        constexpr auto next4 = NextMultinomial(next3);
        static_assert(next4 == parts[4]);

        constexpr auto next5 = NextMultinomial(next4);
        static_assert(next5 == parts[5]);

        constexpr auto next6 = NextMultinomial(next5);
        static_assert(next6 == parts[6]);

        constexpr auto next7 = NextMultinomial(next6);
        static_assert(next7 == parts[7]);

        constexpr auto next8 = NextMultinomial(next7);
        static_assert(next8 == parts[8]);

        constexpr auto next9 = NextMultinomial(next8);
        static_assert(next9 == parts[9]);

        constexpr auto next10 = NextMultinomial(next9);
        static_assert(next10 == parts[10]);

        constexpr auto next11 = NextMultinomial(next10);
        static_assert(next11 == parts[11]);

        constexpr auto next12 = NextMultinomial(next11);
        static_assert(next12 == parts[12]);

        constexpr auto next13 = NextMultinomial(next12);
        static_assert(next13 == parts[13]);

        constexpr auto next14 = NextMultinomial(next13);
        static_assert(next14 == parts[14]);

        constexpr auto next15 = NextMultinomial(next14);
        static_assert(next15 == parts[15]);

        constexpr auto next16 = NextMultinomial(next15);
        static_assert(next16 == parts[16]);

        constexpr auto next17 = NextMultinomial(next16);
        static_assert(next17 == parts[17]);

        constexpr auto next18 = NextMultinomial(next17);
        static_assert(next18 == parts[18]);

        constexpr auto next19 = NextMultinomial(next18);
        static_assert(next19 == parts[19]);

        constexpr auto last = LastMultinomial<4, 3>();
        static_assert(last == parts.back());
        static_assert(last == next19);

        constexpr auto last2 = NextMultinomial(last);
        static_assert(last == last2);
    }
}

} // namespace adhoc3
