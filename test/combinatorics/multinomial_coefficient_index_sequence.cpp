#include <combinatorics/multinomial_coefficient_index_sequence.hpp>

int main() {

    {
        constexpr auto res = adhoc4::MultinomialSequences<4, 4>();
        static_assert(std::tuple_size_v<decltype(res)> == 35);
        constexpr auto res_check = std::tuple<
            std::index_sequence<4, 0, 0, 0>, std::index_sequence<3, 1, 0, 0>,
            std::index_sequence<3, 0, 1, 0>, std::index_sequence<3, 0, 0, 1>,
            std::index_sequence<2, 2, 0, 0>, std::index_sequence<2, 1, 1, 0>,
            std::index_sequence<2, 1, 0, 1>, std::index_sequence<2, 0, 2, 0>,
            std::index_sequence<2, 0, 1, 1>, std::index_sequence<2, 0, 0, 2>,
            std::index_sequence<1, 3, 0, 0>, std::index_sequence<1, 2, 1, 0>,
            std::index_sequence<1, 2, 0, 1>, std::index_sequence<1, 1, 2, 0>,
            std::index_sequence<1, 1, 1, 1>, std::index_sequence<1, 1, 0, 2>,
            std::index_sequence<1, 0, 3, 0>, std::index_sequence<1, 0, 2, 1>,
            std::index_sequence<1, 0, 1, 2>, std::index_sequence<1, 0, 0, 3>,
            std::index_sequence<0, 4, 0, 0>, std::index_sequence<0, 3, 1, 0>,
            std::index_sequence<0, 3, 0, 1>, std::index_sequence<0, 2, 2, 0>,
            std::index_sequence<0, 2, 1, 1>, std::index_sequence<0, 2, 0, 2>,
            std::index_sequence<0, 1, 3, 0>, std::index_sequence<0, 1, 2, 1>,
            std::index_sequence<0, 1, 1, 2>, std::index_sequence<0, 1, 0, 3>,
            std::index_sequence<0, 0, 4, 0>, std::index_sequence<0, 0, 3, 1>,
            std::index_sequence<0, 0, 2, 2>, std::index_sequence<0, 0, 1, 3>,
            std::index_sequence<0, 0, 0, 4>>{};
        static_assert(std::is_same_v<decltype(res), decltype(res_check)>);
    }

    {
        constexpr auto res = adhoc4::MultinomialSequences<4, 5>();
        static_assert(std::tuple_size_v<decltype(res)> == 56);

        constexpr auto res_check = std::tuple<
            std::index_sequence<5, 0, 0, 0>, std::index_sequence<4, 1, 0, 0>,
            std::index_sequence<4, 0, 1, 0>, std::index_sequence<4, 0, 0, 1>,
            std::index_sequence<3, 2, 0, 0>, std::index_sequence<3, 1, 1, 0>,
            std::index_sequence<3, 1, 0, 1>, std::index_sequence<3, 0, 2, 0>,
            std::index_sequence<3, 0, 1, 1>, std::index_sequence<3, 0, 0, 2>,
            std::index_sequence<2, 3, 0, 0>, std::index_sequence<2, 2, 1, 0>,
            std::index_sequence<2, 2, 0, 1>, std::index_sequence<2, 1, 2, 0>,
            std::index_sequence<2, 1, 1, 1>, std::index_sequence<2, 1, 0, 2>,
            std::index_sequence<2, 0, 3, 0>, std::index_sequence<2, 0, 2, 1>,
            std::index_sequence<2, 0, 1, 2>, std::index_sequence<2, 0, 0, 3>,
            std::index_sequence<1, 4, 0, 0>, std::index_sequence<1, 3, 1, 0>,
            std::index_sequence<1, 3, 0, 1>, std::index_sequence<1, 2, 2, 0>,
            std::index_sequence<1, 2, 1, 1>, std::index_sequence<1, 2, 0, 2>,
            std::index_sequence<1, 1, 3, 0>, std::index_sequence<1, 1, 2, 1>,
            std::index_sequence<1, 1, 1, 2>, std::index_sequence<1, 1, 0, 3>,
            std::index_sequence<1, 0, 4, 0>, std::index_sequence<1, 0, 3, 1>,
            std::index_sequence<1, 0, 2, 2>, std::index_sequence<1, 0, 1, 3>,
            std::index_sequence<1, 0, 0, 4>, std::index_sequence<0, 5, 0, 0>,
            std::index_sequence<0, 4, 1, 0>, std::index_sequence<0, 4, 0, 1>,
            std::index_sequence<0, 3, 2, 0>, std::index_sequence<0, 3, 1, 1>,
            std::index_sequence<0, 3, 0, 2>, std::index_sequence<0, 2, 3, 0>,
            std::index_sequence<0, 2, 2, 1>, std::index_sequence<0, 2, 1, 2>,
            std::index_sequence<0, 2, 0, 3>, std::index_sequence<0, 1, 4, 0>,
            std::index_sequence<0, 1, 3, 1>, std::index_sequence<0, 1, 2, 2>,
            std::index_sequence<0, 1, 1, 3>, std::index_sequence<0, 1, 0, 4>,
            std::index_sequence<0, 0, 5, 0>, std::index_sequence<0, 0, 4, 1>,
            std::index_sequence<0, 0, 3, 2>, std::index_sequence<0, 0, 2, 3>,
            std::index_sequence<0, 0, 1, 4>, std::index_sequence<0, 0, 0, 5>>{};
        static_assert(std::is_same_v<decltype(res), decltype(res_check)>);
    }
    return 0;
}
