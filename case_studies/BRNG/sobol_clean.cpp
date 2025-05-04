#include "../../test/test_tools.hpp"
#include "test_tools_rng.hpp"

#include "sobol_engine.hpp"

#include <algorithm>
#include <iostream>
#include <random>
#include <vector>

auto mask(std::size_t w) -> std::uint32_t {
    if (w == 0) {
        return 0;
    }
    std::uint32_t ach = static_cast<std::uint32_t>(1U)
                        << static_cast<std::size_t>(w - 1U);

    return ach | (ach - static_cast<std::uint32_t>(1U));

    // equivalent to the following but shift does not overflow
    // return (static_cast<result_type>(1U) << w) - 1U;
}

int main() {
    TEST_START;

    auto test = [&_result]<adhoc::qrng_table t>() {
        std::uint32_t current_res = 0;
        std::vector<std::uint32_t> res;
        std::size_t count = 0;

        constexpr unsigned int max_dim = adhoc::detail::max_dimension(t);

        for (std::size_t i = 0; i < max_dim; i++) {
            auto const poly = adhoc::detail::poly(t, i);
            std::size_t const poly_d = adhoc::highestBitIndex(poly);

            for (std::size_t ii = 0; ii < poly_d; ++ii) {
                std::uint32_t const init = adhoc::detail::v_init<t>(i, ii);

                std::uint32_t const init_shift = init >> 1;

                if (ii) {
                    std::size_t bits = ii;
                    std::size_t shift = 0;
                    while (bits) {
                        std::size_t rembits = 32U - count;
                        std::size_t bits_here = std::min(rembits, bits);
                        std::uint32_t m = mask(bits_here);
                        std::uint32_t add = ((init_shift >> shift) & m)
                                            << count;

                        current_res |= add;
                        count += bits_here;
                        bits -= bits_here;
                        shift += bits_here;

                        if (count == 32) {
                            count = 0;
                            res.push_back(current_res);
                            current_res = 0;
                        }
                    }
                }

                // this->v[i][ii] = detail::v_init<t>(i, ii);
            }
        }

        if (count != 0) {
            res.push_back(current_res);
        }

        for (std::size_t ii = 0; ii < res.size(); ++ii) {
            std::cout << res[ii] << ", ";
        }
        std::cout << std::endl;

        std::size_t idx1 = 0;
        std::size_t idx2 = 0;
        std::uint32_t val = res[idx1++];
        for (std::size_t i = 0; i < max_dim; i++) {
            auto const poly = adhoc::detail::poly(t, i);
            std::size_t const poly_d = adhoc::highestBitIndex(poly);

            for (std::size_t ii = 0; ii < poly_d; ++ii) {
                std::uint32_t const init = adhoc::detail::v_init<t>(i, ii);

                std::uint32_t init2 = 0;
                {
                    std::size_t bits = ii;
                    std::size_t shift = 0;
                    std::size_t rembits = 32U - idx2;
                    while (bits) {
                        std::size_t bits_here = std::min(rembits, bits);
                        std::uint32_t m = mask(bits_here);
                        init2 |= (val & m) << shift;
                        val >>= bits_here;
                        idx2 += bits_here;
                        bits -= bits_here;
                        shift += bits_here;

                        if (idx2 == 32) {
                            idx2 = 0;
                            val = res[idx1++];
                        }
                        rembits = 32U - idx2;
                    }

                    init2 = (init2 << 1U) | 1U;
                }

                EXPECT_EQUAL(init, init2);
            }
        }
    };

    // test.template operator()<adhoc::qrng_table::joe_kuo_old_1111>();
    // test.template operator()<adhoc::qrng_table::joe_kuo_other_0_7600>();
    // test.template operator()<adhoc::qrng_table::joe_kuo_other_2_3900>();
    // test.template operator()<adhoc::qrng_table::joe_kuo_other_3_7300>();
    // test.template operator()<adhoc::qrng_table::joe_kuo_other_4_5600>();
    // test.template operator()<adhoc::qrng_table::new_joe_kuo_5_21201>();
    // test.template operator()<adhoc::qrng_table::new_joe_kuo_6_21201>();
    test.template operator()<adhoc::qrng_table::new_joe_kuo_7_21201>();

    TEST_END;
}
