#include <combinatorics/combinations.hpp>
#include <combinatorics/multinomial_coefficient.hpp>
#include <combinatorics/multinomial_coefficient_index_sequence.hpp>
#include <combinatorics/trinomial.hpp>
#include <dependency.hpp>
#include <utils/tuple.hpp>

#include <tuple>
#include <type_traits>
#include <utility>

int main() {
    using namespace adhoc4;

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

    {
        constexpr std::array<std::array<std::size_t, 1>, combinations(1, 1)>
            parts{{{1}}};

        constexpr auto factorials = factorialarr<1>();
        static_assert(MultinomialCoeff<1, 1>(factorials, parts[0]) == 1);
    }

    {
        constexpr std::array<std::array<std::size_t, 1>, combinations(1, 2)>
            parts{{{2}}};

        constexpr auto factorials = factorialarr<2>();
        static_assert(MultinomialCoeff<2, 1>(factorials, parts[0]) == 1);
    }

    {
        constexpr std::array<std::array<std::size_t, 1>, combinations(1, 3)>
            parts{{{3}}};

        constexpr auto factorials = factorialarr<3>();
        static_assert(MultinomialCoeff<3, 1>(factorials, parts[0]) == 1);
    }

    {
        constexpr std::array<std::array<std::size_t, 2>, combinations(2, 1)>
            parts{{{1, 0}, {0, 1}}};

        constexpr auto factorials = factorialarr<1>();
        static_assert(MultinomialCoeff<1, 2>(factorials, parts[0]) == 1);
        static_assert(MultinomialCoeff<1, 2>(factorials, parts[1]) == 1);
    }

    {
        constexpr std::array<std::array<std::size_t, 2>, combinations(2, 2)>
            parts{{{2, 0}, {1, 1}, {0, 2}}};

        constexpr auto factorials = factorialarr<2>();
        static_assert(MultinomialCoeff<2, 2>(factorials, parts[0]) == 1);
        static_assert(MultinomialCoeff<2, 2>(factorials, parts[1]) == 2);
        static_assert(MultinomialCoeff<2, 2>(factorials, parts[2]) == 1);
    }

    {
        constexpr std::array<std::array<std::size_t, 2>, combinations(2, 3)>
            parts{{{3, 0}, {2, 1}, {1, 2}, {0, 3}}};

        constexpr auto factorials = factorialarr<3>();
        static_assert(MultinomialCoeff<3, 2>(factorials, parts[0]) == 1);
        static_assert(MultinomialCoeff<3, 2>(factorials, parts[1]) == 3);
        static_assert(MultinomialCoeff<3, 2>(factorials, parts[2]) == 3);
        static_assert(MultinomialCoeff<3, 2>(factorials, parts[3]) == 1);
    }

    {
        constexpr std::array<std::array<std::size_t, 4>, combinations(4, 3)>
            parts{{{3, 0, 0, 0}, {2, 1, 0, 0}, {2, 0, 1, 0}, {2, 0, 0, 1},
                   {1, 2, 0, 0}, {1, 1, 1, 0}, {1, 1, 0, 1}, {1, 0, 2, 0},
                   {1, 0, 1, 1}, {1, 0, 0, 2}, {0, 3, 0, 0}, {0, 2, 1, 0},
                   {0, 2, 0, 1}, {0, 1, 2, 0}, {0, 1, 1, 1}, {0, 1, 0, 2},
                   {0, 0, 3, 0}, {0, 0, 2, 1}, {0, 0, 1, 2}, {0, 0, 0, 3}}};

        constexpr auto factorials = factorialarr<3>();
        static_assert(MultinomialCoeff<3, 4>(factorials, parts[0]) == 1);
        static_assert(MultinomialCoeff<3, 4>(factorials, parts[1]) == 3);
        static_assert(MultinomialCoeff<3, 4>(factorials, parts[2]) == 3);
        static_assert(MultinomialCoeff<3, 4>(factorials, parts[3]) == 3);
        static_assert(MultinomialCoeff<3, 4>(factorials, parts[4]) == 3);
        static_assert(MultinomialCoeff<3, 4>(factorials, parts[5]) == 6);
        static_assert(MultinomialCoeff<3, 4>(factorials, parts[6]) == 6);
        static_assert(MultinomialCoeff<3, 4>(factorials, parts[7]) == 3);
        static_assert(MultinomialCoeff<3, 4>(factorials, parts[8]) == 6);
        static_assert(MultinomialCoeff<3, 4>(factorials, parts[9]) == 3);
        static_assert(MultinomialCoeff<3, 4>(factorials, parts[10]) == 1);
        static_assert(MultinomialCoeff<3, 4>(factorials, parts[11]) == 3);
        static_assert(MultinomialCoeff<3, 4>(factorials, parts[12]) == 3);
        static_assert(MultinomialCoeff<3, 4>(factorials, parts[13]) == 3);
        static_assert(MultinomialCoeff<3, 4>(factorials, parts[14]) == 6);
        static_assert(MultinomialCoeff<3, 4>(factorials, parts[15]) == 3);
        static_assert(MultinomialCoeff<3, 4>(factorials, parts[16]) == 1);
        static_assert(MultinomialCoeff<3, 4>(factorials, parts[17]) == 3);
        static_assert(MultinomialCoeff<3, 4>(factorials, parts[18]) == 3);
        static_assert(MultinomialCoeff<3, 4>(factorials, parts[19]) == 1);
    }

    {
        constexpr auto parts = std::tuple<std::index_sequence<1>>{};
        static_assert(MultinomialCoeff2(std::get<0>(parts)) == 1);
    }

    {
        constexpr auto parts = std::tuple<std::index_sequence<2>>{};
        static_assert(MultinomialCoeff2(std::get<0>(parts)) == 1);
    }

    {
        constexpr auto parts = std::tuple<std::index_sequence<3>>{};
        static_assert(MultinomialCoeff2(std::get<0>(parts)) == 1);
    }

    {
        constexpr auto parts =
            std::tuple<std::index_sequence<1, 0>, std::index_sequence<0, 1>>{};
        static_assert(MultinomialCoeff2(std::get<0>(parts)) == 1);
        static_assert(MultinomialCoeff2(std::get<1>(parts)) == 1);
    }

    {
        constexpr auto parts =
            std::tuple<std::index_sequence<2, 0>, std::index_sequence<1, 1>,
                       std::index_sequence<0, 2>>{};
        static_assert(MultinomialCoeff2(std::get<0>(parts)) == 1);
        static_assert(MultinomialCoeff2(std::get<1>(parts)) == 2);
        static_assert(MultinomialCoeff2(std::get<2>(parts)) == 1);
    }

    {
        constexpr auto parts =
            std::tuple<std::index_sequence<3, 0>, std::index_sequence<2, 1>,
                       std::index_sequence<1, 2>, std::index_sequence<0, 3>>{};
        static_assert(MultinomialCoeff2(std::get<0>(parts)) == 1);
        static_assert(MultinomialCoeff2(std::get<1>(parts)) == 3);
        static_assert(MultinomialCoeff2(std::get<2>(parts)) == 3);
        static_assert(MultinomialCoeff2(std::get<3>(parts)) == 1);
    }

    {
        constexpr auto parts = std::tuple<
            std::index_sequence<3, 0, 0, 0>, std::index_sequence<2, 1, 0, 0>,
            std::index_sequence<2, 0, 1, 0>, std::index_sequence<2, 0, 0, 1>,
            std::index_sequence<1, 2, 0, 0>, std::index_sequence<1, 1, 1, 0>,
            std::index_sequence<1, 1, 0, 1>, std::index_sequence<1, 0, 2, 0>,
            std::index_sequence<1, 0, 1, 1>, std::index_sequence<1, 0, 0, 2>,
            std::index_sequence<0, 3, 0, 0>, std::index_sequence<0, 2, 1, 0>,
            std::index_sequence<0, 2, 0, 1>, std::index_sequence<0, 1, 2, 0>,
            std::index_sequence<0, 1, 1, 1>, std::index_sequence<0, 1, 0, 2>,
            std::index_sequence<0, 0, 3, 0>, std::index_sequence<0, 0, 2, 1>,
            std::index_sequence<0, 0, 1, 2>, std::index_sequence<0, 0, 0, 3>>{};
        static_assert(MultinomialCoeff2(std::get<0>(parts)) == 1);
        static_assert(MultinomialCoeff2(std::get<1>(parts)) == 3);
        static_assert(MultinomialCoeff2(std::get<2>(parts)) == 3);
        static_assert(MultinomialCoeff2(std::get<3>(parts)) == 3);
        static_assert(MultinomialCoeff2(std::get<4>(parts)) == 3);
        static_assert(MultinomialCoeff2(std::get<5>(parts)) == 6);
        static_assert(MultinomialCoeff2(std::get<6>(parts)) == 6);
        static_assert(MultinomialCoeff2(std::get<7>(parts)) == 3);
        static_assert(MultinomialCoeff2(std::get<8>(parts)) == 6);
        static_assert(MultinomialCoeff2(std::get<9>(parts)) == 3);
        static_assert(MultinomialCoeff2(std::get<10>(parts)) == 1);
        static_assert(MultinomialCoeff2(std::get<11>(parts)) == 3);
        static_assert(MultinomialCoeff2(std::get<12>(parts)) == 3);
        static_assert(MultinomialCoeff2(std::get<13>(parts)) == 3);
        static_assert(MultinomialCoeff2(std::get<14>(parts)) == 6);
        static_assert(MultinomialCoeff2(std::get<15>(parts)) == 3);
        static_assert(MultinomialCoeff2(std::get<16>(parts)) == 1);
        static_assert(MultinomialCoeff2(std::get<17>(parts)) == 3);
        static_assert(MultinomialCoeff2(std::get<18>(parts)) == 3);
        static_assert(MultinomialCoeff2(std::get<19>(parts)) == 1);
    }

    {
        constexpr auto multinomial_sequences = TrinomialSequencesMult<4, 6>();

        constexpr auto result = std::make_tuple(
            std::index_sequence<2, 2, 0>{}, std::index_sequence<1, 3, 0>{},
            std::index_sequence<0, 4, 0>{}, std::index_sequence<2, 1, 1>{},
            std::index_sequence<1, 2, 1>{}, std::index_sequence<0, 3, 1>{},
            std::index_sequence<2, 0, 2>{}, std::index_sequence<1, 1, 2>{},
            std::index_sequence<0, 2, 2>{}, std::index_sequence<1, 0, 3>{},
            std::index_sequence<0, 1, 3>{}, std::index_sequence<0, 0, 4>{});
        static_assert(
            std::is_same_v<decltype(multinomial_sequences), decltype(result)>);
    }

    return 0;
}
