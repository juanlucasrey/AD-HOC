#include "../../test/test_tools.hpp"
#include "test_tools_rng.hpp"

#include "sobol_engine.hpp"

#include <random>
#include <vector>

int main() {
    TEST_START;

    // first vals 1
    {
        std::vector<std::uint32_t> results{
            2147483648, 3221225472, 1073741824, 1610612736, 3758096384,
            2684354560, 536870912,  805306368,  2952790016, 4026531840,
            1879048192, 1342177280, 3489660928, 2415919104, 268435456,
            402653184,  2550136832, 3623878656, 1476395008, 2013265920,
            2147483648, 2147483648, 3221225472, 1073741824, 1073741824,
            3221225472, 1610612736, 1610612736, 3758096384, 3758096384,
            2684354560, 536870912,  536870912,  2684354560, 805306368,
            1342177280, 2952790016, 3489660928, 4026531840, 268435456,
            2147483648, 2147483648, 2147483648, 3221225472, 1073741824,
            3221225472, 1073741824, 3221225472, 1073741824, 1610612736,
            1610612736, 2684354560, 3758096384, 3758096384, 536870912,
            2684354560, 536870912,  1610612736, 536870912,  2684354560};

        std::size_t nsims = 20;
        std::size_t ndims = 4;
        std::size_t counter = 0;
        for (std::size_t i = 1; i < ndims; i++) {
            adhoc::sobol_engine<std::uint32_t, 32, true,
                                adhoc::qrng_table::joe_kuo_old_1111>
                rng1(i);
            adhoc::sobol_engine<std::uint32_t, 32, true,
                                adhoc::qrng_table::joe_kuo_other_0_7600>
                rng2(i);

            adhoc::sobol_engine<std::uint32_t, 32, true,
                                adhoc::qrng_table::joe_kuo_other_2_3900>
                rng3(i);

            adhoc::sobol_engine<std::uint32_t, 32, true,
                                adhoc::qrng_table::joe_kuo_other_3_7300>
                rng4(i);

            adhoc::sobol_engine<std::uint32_t, 32, true,
                                adhoc::qrng_table::joe_kuo_other_4_5600>
                rng5(i);

            for (std::size_t j = 0; j < nsims; j++) {
                auto val = results[counter++];
                EXPECT_EQUAL(rng1(), val);
                EXPECT_EQUAL(rng2(), val);
                EXPECT_EQUAL(rng3(), val);
                EXPECT_EQUAL(rng4(), val);
                EXPECT_EQUAL(rng5(), val);
            }
        }
    }

    // first vals 2
    {
        std::vector<std::uint32_t> results{
            2147483648, 3221225472, 1073741824, 1610612736, 3758096384,
            2684354560, 536870912,  805306368,  2952790016, 4026531840,
            1879048192, 1342177280, 3489660928, 2415919104, 268435456,
            402653184,  2550136832, 3623878656, 1476395008, 2013265920,
            2147483648, 2147483648, 3221225472, 1073741824, 1073741824,
            3221225472, 1610612736, 1610612736, 3758096384, 3758096384,
            2684354560, 536870912,  536870912,  2684354560, 805306368,
            1342177280, 2952790016, 3489660928, 4026531840, 268435456,
            2147483648, 2147483648, 2147483648, 3221225472, 1073741824,
            1073741824, 1073741824, 3221225472, 3221225472, 1610612736,
            1610612736, 2684354560, 3758096384, 3758096384, 536870912,
            2684354560, 536870912,  3758096384, 536870912,  2684354560};

        std::size_t nsims = 20;
        std::size_t ndims = 4;
        std::size_t counter = 0;
        for (std::size_t i = 1; i < ndims; i++) {
            adhoc::sobol_engine<std::uint32_t, 32, true,
                                adhoc::qrng_table::new_joe_kuo_5_21201>
                rng6(i);
            adhoc::sobol_engine<std::uint32_t, 32, true,
                                adhoc::qrng_table::new_joe_kuo_6_21201>
                rng7(i);
            adhoc::sobol_engine<std::uint32_t, 32, true,
                                adhoc::qrng_table::new_joe_kuo_7_21201>
                rng8(i);

            for (std::size_t j = 0; j < nsims; j++) {
                auto val = results[counter++];
                EXPECT_EQUAL(rng6(), val);
                EXPECT_EQUAL(rng7(), val);
                EXPECT_EQUAL(rng8(), val);
            }
        }
    }

    // all different values check
    {
        std::vector<std::uint32_t> results{4160749568, 1744830464, 134217728,
                                           402653184,  939524096,  2818572288,
                                           2281701376, 2550136832};

        std::size_t ndims = 15;
        adhoc::sobol_engine<std::uint32_t, 32, true,
                            adhoc::qrng_table::joe_kuo_old_1111>
            rng1(ndims);
        adhoc::sobol_engine<std::uint32_t, 32, true,
                            adhoc::qrng_table::joe_kuo_other_0_7600>
            rng2(ndims);

        adhoc::sobol_engine<std::uint32_t, 32, true,
                            adhoc::qrng_table::joe_kuo_other_2_3900>
            rng3(ndims);

        adhoc::sobol_engine<std::uint32_t, 32, true,
                            adhoc::qrng_table::joe_kuo_other_3_7300>
            rng4(ndims);

        adhoc::sobol_engine<std::uint32_t, 32, true,
                            adhoc::qrng_table::joe_kuo_other_4_5600>
            rng5(ndims);
        adhoc::sobol_engine<std::uint32_t, 32, true,
                            adhoc::qrng_table::new_joe_kuo_5_21201>
            rng6(ndims);
        adhoc::sobol_engine<std::uint32_t, 32, true,
                            adhoc::qrng_table::new_joe_kuo_6_21201>
            rng7(ndims);
        adhoc::sobol_engine<std::uint32_t, 32, true,
                            adhoc::qrng_table::new_joe_kuo_7_21201>
            rng8(ndims);

        std::size_t nsims = 359;
        for (std::size_t j = 0; j < nsims; j++) {
            rng1();
            rng2();
            rng3();
            rng4();
            rng5();
            rng6();
            rng7();
            rng8();
        }

        EXPECT_EQUAL(rng1(), results[0]);
        EXPECT_EQUAL(rng2(), results[1]);
        EXPECT_EQUAL(rng3(), results[2]);
        EXPECT_EQUAL(rng4(), results[3]);
        EXPECT_EQUAL(rng5(), results[4]);
        EXPECT_EQUAL(rng6(), results[5]);
        EXPECT_EQUAL(rng7(), results[6]);
        EXPECT_EQUAL(rng8(), results[7]);
    }

    // initial value check (skip first)
    {
        auto check_first = [&_result]<adhoc::qrng_table T>() {
            std::size_t ndims = 200;
            for (std::size_t i = 1; i < ndims; i++) {
                adhoc::sobol_engine<std::uint32_t, 32, true, T> rng(i);

                for (std::size_t j = 0; j < i; j++) {
                    auto val = rng();
                    EXPECT_EQUAL(val, 2147483648);
                }
            }
        };

        check_first.template operator()<adhoc::qrng_table::joe_kuo_old_1111>();
        check_first
            .template operator()<adhoc::qrng_table::joe_kuo_other_0_7600>();
        check_first
            .template operator()<adhoc::qrng_table::joe_kuo_other_2_3900>();
        check_first
            .template operator()<adhoc::qrng_table::joe_kuo_other_3_7300>();
        check_first
            .template operator()<adhoc::qrng_table::joe_kuo_other_4_5600>();
        check_first
            .template operator()<adhoc::qrng_table::new_joe_kuo_5_21201>();
        check_first
            .template operator()<adhoc::qrng_table::new_joe_kuo_6_21201>();
        check_first
            .template operator()<adhoc::qrng_table::new_joe_kuo_7_21201>();
    }

    // initial value check (no skip first)
    {
        auto check_first = [&_result]<adhoc::qrng_table T>() {
            std::size_t ndims = 200;
            for (std::size_t i = 1; i < ndims; i++) {
                adhoc::sobol_engine<std::uint32_t, 32, false, T> rng(i);

                for (std::size_t j = 0; j < i; j++) {
                    auto val = rng();
                    EXPECT_EQUAL(val, 0);
                }

                for (std::size_t j = 0; j < i; j++) {
                    auto val = rng();
                    EXPECT_EQUAL(val, 2147483648);
                }
            }
        };

        check_first.template operator()<adhoc::qrng_table::joe_kuo_old_1111>();
        check_first
            .template operator()<adhoc::qrng_table::joe_kuo_other_0_7600>();
        check_first
            .template operator()<adhoc::qrng_table::joe_kuo_other_2_3900>();
        check_first
            .template operator()<adhoc::qrng_table::joe_kuo_other_3_7300>();
        check_first
            .template operator()<adhoc::qrng_table::joe_kuo_other_4_5600>();
        check_first
            .template operator()<adhoc::qrng_table::new_joe_kuo_5_21201>();
        check_first
            .template operator()<adhoc::qrng_table::new_joe_kuo_6_21201>();
        check_first
            .template operator()<adhoc::qrng_table::new_joe_kuo_7_21201>();
    }

    // initial discard
    {
        auto check_discard = [&_result]<adhoc::qrng_table T, bool First>() {
            std::size_t ndims = 30;
            for (std::size_t i = 1; i < ndims; i++) {
                adhoc::sobol_engine<std::uint32_t, 32, true,
                                    adhoc::qrng_table::joe_kuo_old_1111>
                    gen(i);

                for (unsigned int j = 1; j < 1000; j++) {
                    gen();

                    adhoc::sobol_engine<std::uint32_t, 32, true,
                                        adhoc::qrng_table::joe_kuo_old_1111>
                        gen2(i);

                    gen2.discard(j);
                    EXPECT_EQUAL(gen, gen2);
                }
            }
        };

        check_discard
            .template operator()<adhoc::qrng_table::joe_kuo_old_1111, true>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_other_0_7600, true>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_other_2_3900, true>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_other_3_7300, true>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_other_4_5600, true>();
        check_discard.template
        operator()<adhoc::qrng_table::new_joe_kuo_5_21201, true>();
        check_discard.template
        operator()<adhoc::qrng_table::new_joe_kuo_6_21201, true>();
        check_discard.template
        operator()<adhoc::qrng_table::new_joe_kuo_7_21201, true>();

        check_discard
            .template operator()<adhoc::qrng_table::joe_kuo_old_1111, false>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_other_0_7600, false>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_other_2_3900, false>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_other_3_7300, false>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_other_4_5600, false>();
        check_discard.template
        operator()<adhoc::qrng_table::new_joe_kuo_5_21201, false>();
        check_discard.template
        operator()<adhoc::qrng_table::new_joe_kuo_6_21201, false>();
        check_discard.template
        operator()<adhoc::qrng_table::new_joe_kuo_7_21201, false>();
    }

    // full cycle discard
    {
        auto check_discard =
            [&_result]<adhoc::qrng_table T, bool First, std::size_t w = 32>() {
                std::size_t ndims = 200;
                for (std::size_t i = 1; i < ndims; i++) {
                    adhoc::sobol_engine<std::uint32_t, w, First, T> rng1(i);
                    adhoc::sobol_engine<std::uint32_t, w, First, T> rng2(i);

                    rng1.discard((static_cast<unsigned long long>(
                                      adhoc::sobol_engine<std::uint32_t, w,
                                                          First, T>::max()) +
                                  1) *
                                 i);
                    EXPECT_EQUAL(rng1, rng2);
                }
            };

        check_discard
            .template operator()<adhoc::qrng_table::joe_kuo_old_1111, false>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_old_1111, false, 31>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_old_1111, false, 30>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_old_1111, false, 23>();

        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_other_0_7600, false>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_other_2_3900, false>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_other_3_7300, false>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_other_4_5600, false>();
        check_discard.template
        operator()<adhoc::qrng_table::new_joe_kuo_5_21201, false>();
        check_discard.template
        operator()<adhoc::qrng_table::new_joe_kuo_6_21201, false>();
        check_discard.template
        operator()<adhoc::qrng_table::new_joe_kuo_7_21201, false>();

        check_discard
            .template operator()<adhoc::qrng_table::joe_kuo_old_1111, true>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_other_0_7600, true>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_other_2_3900, true>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_other_3_7300, true>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_other_4_5600, true>();
        check_discard.template
        operator()<adhoc::qrng_table::new_joe_kuo_5_21201, true>();
        check_discard.template
        operator()<adhoc::qrng_table::new_joe_kuo_6_21201, true>();
        check_discard.template
        operator()<adhoc::qrng_table::new_joe_kuo_7_21201, true>();
    }

    // right before full cycle discard
    {
        auto check_discard =
            [&_result]<adhoc::qrng_table T, bool First, std::size_t w = 32>() {
                std::size_t ndims = 200;
                for (std::size_t i = 180; i < ndims; i++) {
                    // std::size_t i = 180;
                    adhoc::sobol_engine<std::uint32_t, w, First, T> rng1(i);
                    adhoc::sobol_engine<std::uint32_t, w, First, T> rng2(i);

                    rng1.discard(((static_cast<unsigned long long>(
                                       adhoc::sobol_engine<std::uint32_t, w,
                                                           First, T>::max()) +
                                   1) *
                                  i) -
                                 10);

                    for (std::size_t j = 0; j < 10; j++) {
                        rng1();
                    }
                    EXPECT_EQUAL(rng1, rng2);
                }
            };

        check_discard
            .template operator()<adhoc::qrng_table::joe_kuo_old_1111, false>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_old_1111, false, 31>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_old_1111, false, 30>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_old_1111, false, 29>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_old_1111, false, 23>();

        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_other_0_7600, false>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_other_2_3900, false>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_other_3_7300, false>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_other_4_5600, false>();
        check_discard.template
        operator()<adhoc::qrng_table::new_joe_kuo_5_21201, false>();
        check_discard.template
        operator()<adhoc::qrng_table::new_joe_kuo_6_21201, false>();
        check_discard.template
        operator()<adhoc::qrng_table::new_joe_kuo_7_21201, false>();

        check_discard
            .template operator()<adhoc::qrng_table::joe_kuo_old_1111, true>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_other_0_7600, true>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_other_2_3900, true>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_other_3_7300, true>();
        check_discard.template
        operator()<adhoc::qrng_table::joe_kuo_other_4_5600, true>();
        check_discard.template
        operator()<adhoc::qrng_table::new_joe_kuo_5_21201, true>();
        check_discard.template
        operator()<adhoc::qrng_table::new_joe_kuo_6_21201, true>();
        check_discard.template
        operator()<adhoc::qrng_table::new_joe_kuo_7_21201, true>();
    }

    // 20 to 31 bits val check
    {
        auto full_range = []<class Generator>() {
            if constexpr (std::is_integral_v<std::invoke_result_t<Generator>>) {
                return static_cast<double>(Generator::max()) + 1.0;
            } else {
                return 1.0;
            }
        };

        auto generate_canonical =
            [&full_range]<class RealType, class Generator>(Generator &g) {
                constexpr auto fullrange =
                    full_range.template operator()<Generator>();

                constexpr RealType step =
                    (static_cast<RealType>(fullrange) -
                     (static_cast<RealType>(Generator::max()) -
                      static_cast<RealType>(Generator::min()))) /
                    static_cast<RealType>(2.0);

                constexpr RealType denominator =
                    ((static_cast<RealType>(1.0) * step *
                      static_cast<RealType>(2.0)) +
                     static_cast<RealType>(Generator::max()) -
                     static_cast<RealType>(Generator::min()));

                constexpr RealType one_over_range =
                    static_cast<RealType>(1.0) / denominator;

                constexpr RealType left_shift =
                    -static_cast<RealType>(Generator::min());
                return (static_cast<RealType>(g()) + left_shift) *
                       one_over_range;
            };

        auto check_uniform_vals =
            [&_result, &generate_canonical]<adhoc::qrng_table T, bool First,
                                            std::size_t w>() {
                std::size_t ndims = 20;
                adhoc::sobol_engine<std::uint32_t, w, true,
                                    adhoc::qrng_table::joe_kuo_old_1111>
                    rng(ndims);

                adhoc::sobol_engine<std::uint32_t, 32, true,
                                    adhoc::qrng_table::joe_kuo_old_1111>
                    rng32(ndims);

                for (std::size_t i = 0; i < 1000; i++) {
                    auto val1 =
                        generate_canonical.template operator()<float>(rng32);
                    auto val2 =
                        generate_canonical.template operator()<float>(rng);
                    EXPECT_EQUAL(val1, val2);
                }
            };

        auto check_uniform_vals_mult =
            [&_result, &check_uniform_vals]<adhoc::qrng_table T, bool First>() {
                check_uniform_vals.template operator()<T, First, 31>();
                check_uniform_vals.template operator()<T, First, 30>();
                check_uniform_vals.template operator()<T, First, 29>();
                check_uniform_vals.template operator()<T, First, 28>();
                check_uniform_vals.template operator()<T, First, 27>();
                check_uniform_vals.template operator()<T, First, 26>();
                check_uniform_vals.template operator()<T, First, 25>();
                check_uniform_vals.template operator()<T, First, 24>();
                check_uniform_vals.template operator()<T, First, 23>();
                check_uniform_vals.template operator()<T, First, 22>();
                check_uniform_vals.template operator()<T, First, 21>();
                check_uniform_vals.template operator()<T, First, 20>();
            };

        check_uniform_vals_mult
            .template operator()<adhoc::qrng_table::joe_kuo_old_1111, false>();
        check_uniform_vals_mult.template
        operator()<adhoc::qrng_table::joe_kuo_other_0_7600, false>();
        check_uniform_vals_mult.template
        operator()<adhoc::qrng_table::joe_kuo_other_2_3900, false>();
        check_uniform_vals_mult.template
        operator()<adhoc::qrng_table::joe_kuo_other_3_7300, false>();
        check_uniform_vals_mult.template
        operator()<adhoc::qrng_table::joe_kuo_other_4_5600, false>();
        check_uniform_vals_mult.template
        operator()<adhoc::qrng_table::new_joe_kuo_5_21201, false>();
        check_uniform_vals_mult.template
        operator()<adhoc::qrng_table::new_joe_kuo_6_21201, false>();
        check_uniform_vals_mult.template
        operator()<adhoc::qrng_table::new_joe_kuo_7_21201, false>();
    }

    // {
    //     std::size_t dims = 30;
    //     adhoc::sobol_engine<std::uint32_t, 32, true,
    //                         adhoc::qrng_table::new_joe_kuo_7_21201>
    //         rng1(dims);

    //     for (unsigned int j = 0; j < dims; j++) {
    //         std::cout << rng1() << std::endl;
    //     }
    //     std::cout << "end first" << std::endl;

    //     auto discardsize = 4294967294 * dims - 10;

    //     rng1.discard(discardsize);

    //     for (unsigned int j = 0; j < 20; j++) {
    //         auto val = rng1();
    //         // if (val == first) {
    //         //   std::cout << "first" << std::endl;
    //         // }
    //         std::cout << val << std::endl;
    //     }

    //     std::cout << std::endl;

    //     for (unsigned int j = 0; j < dims; j++) {
    //         std::cout << rng1() << std::endl;
    //     }

    //     for (unsigned int j = 0; j < dims; j++) {
    //         std::cout << rng1() << std::endl;
    //     }
    // }

    // {
    //     auto ndims = 1111;
    //     // adhoc::sobol_engine<std::uint32_t, 32, true,
    //     //                     adhoc::qrng_table::new_joe_kuo_5_21201>
    //     //     rng1(ndims);
    //     adhoc::sobol_engine<std::uint32_t, 32, true,
    //                         adhoc::qrng_table::joe_kuo_old_1111>
    //         rng2(ndims);

    //     // std::cout << "32" << std::endl;
    //     // rng1.discard(4294967294 * 3 / 2 - 10);
    //     // for (std::size_t i = 1; i < 20; i++) {
    //     //     auto val = rng1();
    //     //     std::bitset<32> valbitset(val);
    //     //     std::cout << valbitset << std::endl;
    //     // }
    //     // std::cout << std::endl;

    //     std::cout << "31" << std::endl;
    //     // auto disc = 4294967294 * ndims / 2 - 6;
    //     // auto disc = 4294967294 * ndims / 536870912 - 6;
    //     auto disc = 32768 * ndims - 6;
    //     std::cout << disc << std::endl;

    //     std::uint32_t mask = 1 << 20;
    //     std::bitset<32> maskbitset(mask);
    //     std::cout << maskbitset << std::endl;
    //     std::size_t i = 0;
    //     while (i < 3000000) {
    //         auto val = rng2();

    //         auto valmask = val & mask;
    //         if (valmask) {
    //             std::cout << i << ", ";
    //             std::cout << val << ", ";
    //             std::bitset<32> valbitset(val);
    //             std::cout << valbitset << ", ";
    //             break;
    //         }
    //         i++;
    //     }
    //     std::cout << std::endl;
    //     // for (std::size_t i = 0; i < 200000; i++) {
    //     // }

    //     // std::cout << disc << std::endl;
    //     // rng2.discard(6442450940);
    //     // rng2.discard(disc);
    //     // for (std::size_t i = 1; i < 20; i++) {
    //     //     auto val = rng2();
    //     //     std::bitset<32> valbitset(val);
    //     //     std::cout << valbitset << std::endl;
    //     // }
    //     // std::cout << std::endl;
    // }

    // std::cout << std::numeric_limits<unsigned int>::digits << std::endl;
    // auto const max = std::numeric_limits<std::uint32_t>::max();
    // std::cout << max << std::endl;

    // // std::cout << std::numeric_limits<std::uint32_t>::max() << std::endl;
    // adhoc::sobol_engine<std::uint32_t, 32, false> rng;

    // for (std::size_t i = 0; i < 20; ++i) {
    //     auto valui = rng();
    //     auto val = static_cast<double>(valui) / static_cast<double>(max);
    //     std::cout << valui << std::endl;
    //     std::cout << val << std::endl;
    // }

    // // rng();
    // // adhoc::sobol_engine<std::uint32_t> rng;
    // rng.discard(max - 30);
    // // auto val = static_cast<double>(rng()) / static_cast<double>(max);
    // // std::cout << val << std::endl; // 2147483648

    // for (std::size_t i = 0; i < 20; ++i) {
    //     auto valui = rng();
    //     auto val = static_cast<double>(valui) / static_cast<double>(max);
    //     std::cout << valui << std::endl;
    //     std::cout << val << std::endl;
    // }

    // std::uniform_real_distribution<double> temp;
    // std::cout << temp(rng) << std::endl;

    // auto max = std::numeric_limits<std::uint32_t>::max() - 2;
    // rng.discard(static_cast<uint64_t>(max));

    // std::cout << rng() << std::endl; // 2147483648
    // std::cout << rng() << std::endl; // 3221225472

    TEST_END;
}
