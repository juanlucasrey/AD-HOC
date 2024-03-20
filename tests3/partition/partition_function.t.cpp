#include "../../public3/partition/partition_function.hpp"

#include <gtest/gtest.h>

namespace adhoc3 {

TEST(Partition, PartitionFunction) {
    // https://oeis.org/A000041
    static_assert(partition_function(0) == 1);
    static_assert(partition_function(1) == 1);
    static_assert(partition_function(2) == 2);
    static_assert(partition_function(3) == 3);
    static_assert(partition_function(4) == 5);
    static_assert(partition_function(5) == 7);
    static_assert(partition_function(6) == 11);
    static_assert(partition_function(7) == 15);
    static_assert(partition_function(8) == 22);
    static_assert(partition_function(9) == 30);
}

} // namespace adhoc3
