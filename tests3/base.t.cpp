#include <base.hpp>

#include <gtest/gtest.h>

namespace adhoc3 {

TEST(Base, InvDer) {

    double val = 100.0;
    double res = 1.0 / val;

    {
        auto derivatives = inv_t<double>::d2<1>(res, val);

        EXPECT_NEAR(derivatives[0], -1.0 / (val * val), 1e-15);
    }

    {
        auto derivatives = inv_t<double>::d2<2>(res, val);

        EXPECT_NEAR(derivatives[0], -1.0 / (val * val), 1e-15);
        EXPECT_NEAR(derivatives[1], 2.0 / (val * val * val), 1e-15);
    }

    {
        auto derivatives = inv_t<double>::d2<3>(res, val);

        EXPECT_NEAR(derivatives[0], -1.0 / (val * val), 1e-15);
        EXPECT_NEAR(derivatives[1], 2.0 / (val * val * val), 1e-15);
        EXPECT_NEAR(derivatives[2], -6.0 / (val * val * val * val), 1e-14);
    }

    {
        auto derivatives = inv_t<double>::d2<4>(res, val);

        EXPECT_NEAR(derivatives[0], -1.0 / (val * val), 1e-15);
        EXPECT_NEAR(derivatives[1], 2.0 / (val * val * val), 1e-15);
        EXPECT_NEAR(derivatives[2], -6.0 / (val * val * val * val), 1e-14);
        EXPECT_NEAR(derivatives[3], 24.0 / (val * val * val * val * val),
                    1e-15);
    }
}

} // namespace adhoc3
