#include <utils/index_sequence.hpp>
#include <utils/tuple.hpp>

#include <tuple>

int main() {
    using namespace adhoc4;

    // ReplaceFirst
    {
        constexpr std::tuple<double, int, char, int> temp;
        constexpr int v1 = 0;
        constexpr float v2 = 0;

        constexpr auto newtuple = replace_first(temp, v1, v2);
        constexpr std::tuple<double, float, char, int> res_check;
        static_assert(std::is_same_v<decltype(newtuple), decltype(res_check)>);

        constexpr auto newtuple2 = replace_first(temp, v2, v1);
        constexpr std::tuple<double, int, char, int> res_check2;
        static_assert(
            std::is_same_v<decltype(newtuple2), decltype(res_check2)>);
    }

    // Invert
    {
        constexpr std::tuple<double, int, char> temp;
        constexpr auto inverted_temp = invert(temp);
        constexpr std::tuple<char, int, double> res_check;
        static_assert(
            std::is_same_v<decltype(inverted_temp), decltype(res_check)>);
    }

    // Find
    {

        constexpr std::tuple<double, int, char> temp;
        constexpr auto loc2 = find<std::tuple<double, int, char>, int, 1>();
        static_assert(loc2 == 1);
    }

    // SubTuple
    {

        constexpr std::tuple<double, int, char> temp;
        constexpr auto sub = sub_tuple<1, 2>(temp);
        constexpr std::tuple<int, char> result;
        static_assert(std::is_same_v<decltype(result), decltype(sub)>);
    }

    return 0;
}
