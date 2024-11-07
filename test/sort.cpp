#include <calc_tree.hpp>
#include <sort.hpp>

int main() {
    using namespace adhoc4;

    ADHOC(x);
    ADHOC(y);
    ADHOC(z);
    auto res = x + (y * z);
    CalcTree ct(res);
    using PrimalNodes = decltype(ct)::ValuesTupleInverse;
    constexpr auto primal_nodes_inverted = PrimalNodes{};

    auto d1 = d(res);
    auto d2 = d(res) * d(x);
    auto d3 = d(res) * d(y);
    auto d4 = d(res) * d(z);
    auto d5 = d<2>(res);
    auto d6 = d<2>(res) * d<2>(x);

    constexpr auto restuple = sort_differential_operators(
        std::make_tuple(d4, d3, d2, d1, d5, d6), primal_nodes_inverted);

    constexpr auto res_check = std::make_tuple(d6, d5, d2, d3, d4, d1);
    static_assert(std::is_same_v<decltype(restuple), decltype(res_check)>);
    return 0;
}
