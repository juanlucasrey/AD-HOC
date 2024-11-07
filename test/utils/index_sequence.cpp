#include <utils/index_sequence.hpp>

#include <utility>

int main() {
    using namespace adhoc4;

    // MakeIndexSequence
    {
        {
            constexpr auto seq = make_index_sequence<3, 7>();
            constexpr auto res = std::index_sequence<3, 4, 5, 6, 7>{};

            static_assert(std::is_same_v<decltype(seq), decltype(res)>);
        }

        {
            constexpr auto seq = make_index_sequence<7, 3>();
            constexpr auto res = std::index_sequence<7, 6, 5, 4, 3>{};

            static_assert(std::is_same_v<decltype(seq), decltype(res)>);
        }
    }

    // Replace
    {
        constexpr auto bool_seq = make_bool_sequence<false, 4>();
        constexpr auto repl = replace<2, true>(bool_seq);
        constexpr auto res =
            std::integer_sequence<bool, false, false, true, false>{};
        static_assert(std::is_same_v<decltype(repl), decltype(res)>);
    }

    // Find
    {
        {
            constexpr auto bool_seq =
                std::integer_sequence<bool, false, false, true, false>{};
            static_assert(find<true>(bool_seq) == 2);
        }

        {
            constexpr auto bool_seq =
                std::integer_sequence<bool, false, false, false, false>{};
            static_assert(find<true>(bool_seq) == 4);
        }

        {
            constexpr auto bool_seq =
                std::integer_sequence<bool, true, false, true, false>{};
            static_assert(find<true, 0>(bool_seq) == 0);
            static_assert(find<true, 1>(bool_seq) == 2);
            static_assert(find<true, 2>(bool_seq) == 2);
            static_assert(find<true, 3>(bool_seq) == 4);
        }
    }
    return 0;
}
