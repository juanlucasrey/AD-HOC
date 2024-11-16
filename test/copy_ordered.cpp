#include <backpropagator_tools.hpp>

#include "test_tools.hpp"

#include <tuple>
#include <type_traits>
#include <utility>

int main() {
    TEST_START;
    using namespace adhoc;

    // FilteredCopy
    {

        constexpr std::tuple<std::true_type, std::true_type, std::false_type,
                             std::true_type>
            temp;
        std::array<double, 4> arr{1, 2, 3, 4};
        std::array<double, 3> arr_to;

        detail::copy_filtered(arr, arr_to, temp);
        EXPECT_EQUAL(arr_to[0], 1);
        EXPECT_EQUAL(arr_to[1], 2);
        EXPECT_EQUAL(arr_to[2], 4);
    }

    // FilteredCopyInverse
    {

        {
            constexpr std::integer_sequence<bool, true, false, true, false>
                temp;
            std::array<double, 7> arr{1, 2, 3, 4, 5, 6, 7};
            std::array<double, 2> arr_to;

            detail::copy_filtered_inverted2<1, 5>(arr, arr_to, temp);

            EXPECT_EQUAL(arr_to[0], 5);
            EXPECT_EQUAL(arr_to[1], 3);
        }

        {
            constexpr std::tuple<std::true_type, std::true_type,
                                 std::false_type, std::true_type>
                temp;
            std::array<double, 4> arr{1, 2, 3, 4};
            std::array<double, 3> arr_to;

            detail::copy_filtered_inverted(arr, arr_to, temp);
            EXPECT_EQUAL(arr_to[0], 4);
            EXPECT_EQUAL(arr_to[1], 3);
            EXPECT_EQUAL(arr_to[2], 1);
        }

        {
            constexpr std::tuple<std::true_type, std::true_type> temp;
            std::array<double, 3> arr{1, 2, 3};
            std::array<double, 2> arr_to;

            detail::copy_filtered_inverted(arr, arr_to, temp);

            EXPECT_EQUAL(arr_to[0], 3);
            EXPECT_EQUAL(arr_to[1], 2);
        }

        {
            constexpr std::tuple<std::true_type, std::false_type> temp;
            std::array<double, 3> arr{1, 2, 3};
            std::array<double, 2> arr_to;
            arr_to.fill(0);

            detail::copy_filtered_inverted(arr, arr_to, temp);

            EXPECT_EQUAL(arr_to[0], 3);
            EXPECT_EQUAL(arr_to[1], 0);
        }
    }

    TEST_END;
}
