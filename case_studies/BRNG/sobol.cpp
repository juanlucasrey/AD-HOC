#include "../../test/test_tools.hpp"
#include "test_tools_rng.hpp"

#include "sobol_engine.hpp"

#include <random>
#include <vector>

int main() {
    TEST_START;

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

    {
        auto check_discard = [&_result]<adhoc::qrng_table T, bool First>() {
            std::size_t ndims = 200;
            for (std::size_t i = 1; i < ndims; i++) {
                adhoc::sobol_engine<std::uint32_t, 32, First, T> rng1(i);
                adhoc::sobol_engine<std::uint32_t, 32, First, T> rng2(i);

                rng1.discard(4294967295 * i);
                EXPECT_EQUAL(rng1, rng2);
            }
        };

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

    // {
    //     adhoc::sobol_engine<std::uint32_t, 32> rng1(3);
    //     adhoc::sobol_engine<std::uint32_t, 31> rng2(3);

    //     // std::cout << "32" << std::endl;
    //     // rng1.discard(4294967294 * 3 / 2 - 10);
    //     // for (std::size_t i = 1; i < 20; i++) {
    //     //     auto val = rng1();
    //     //     std::bitset<32> valbitset(val);
    //     //     std::cout << valbitset << std::endl;
    //     // }
    //     // std::cout << std::endl;

    //     std::cout << "31" << std::endl;
    //     // auto disc = 4294967294 * 3 / 2 - 2;
    //     // std::cout << disc << std::endl;
    //     rng2.discard(6442450940);
    //     for (std::size_t i = 1; i < 20; i++) {
    //         auto val = rng2();
    //         std::bitset<32> valbitset(val);
    //         std::cout << valbitset << std::endl;
    //     }
    //     std::cout << std::endl;
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
