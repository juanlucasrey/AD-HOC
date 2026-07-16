#include <tape.hpp>

#include <test_simple_include.hpp>

#include <iostream>
#include <iterator>
#include <tuple>

namespace {

auto
generate_indices(std::size_t n_output, std::size_t lanes)
{

    adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tapeptr;
    tapeptr->configure(adhoc::Method::Bwd, 1, n_output, lanes);
    auto temp = tapeptr->subranges();

    std::vector<std::vector<std::size_t> > result;
    std::vector<std::size_t> result2;
    std::vector<std::vector<std::size_t> > result3;
    for (const auto& metricRange : tapeptr->subranges()) {
        result2.push_back(metricRange.global_index);
        result.emplace_back();
        result3.emplace_back();
        for (const auto& metricIndex : metricRange) {
            result.back().push_back(metricIndex.sub_index);
            result3.back().push_back(metricIndex.global_index);
        }
    }
    return std::make_tuple(std::move(result), std::move(result2), std::move(result3));
}

void
index_tests()
{

    {
        auto [result, result2, result3] = generate_indices(10, 3);
        std::vector<std::vector<std::size_t> > expected = { { 0, 1, 2 }, { 0, 1, 2 }, { 0, 1, 2 }, { 0 } };
        std::vector<std::size_t> expected2 = { 0, 3, 6, 9 };
        std::vector<std::vector<std::size_t> > expected3 = { { 0, 1, 2 }, { 3, 4, 5 }, { 6, 7, 8 }, { 9 } };
        EXPECT_EQUAL(result, expected);
        EXPECT_EQUAL(result2, expected2);
        EXPECT_EQUAL(result3, expected3);
    }

    {
        auto [result, result2, result3] = generate_indices(10, 4);
        std::vector<std::vector<std::size_t> > expected = { { 0, 1, 2, 3 }, { 0, 1, 2, 3 }, { 0, 1 } };
        std::vector<std::size_t> expected2 = { 0, 4, 8 };
        std::vector<std::vector<std::size_t> > expected3 = { { 0, 1, 2, 3 }, { 4, 5, 6, 7 }, { 8, 9 } };
        EXPECT_EQUAL(result, expected);
        EXPECT_EQUAL(result2, expected2);
        EXPECT_EQUAL(result3, expected3);
    }

    {
        auto [result, result2, result3] = generate_indices(10, 5);
        std::vector<std::vector<std::size_t> > expected = { { 0, 1, 2, 3, 4 }, { 0, 1, 2, 3, 4 } };
        std::vector<std::size_t> expected2 = { 0, 5 };
        std::vector<std::vector<std::size_t> > expected3 = { { 0, 1, 2, 3, 4 }, { 5, 6, 7, 8, 9 } };
        EXPECT_EQUAL(result, expected);
        EXPECT_EQUAL(result2, expected2);
        EXPECT_EQUAL(result3, expected3);
    }
}

} // namespace

auto
main() -> int
{
    index_tests();

    TEST_END;
}
