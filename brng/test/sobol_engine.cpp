#include "../../test_simple/test_simple_include.hpp"
#include "test_tools_rng.hpp"

#include <distribution/uniform_distribution.hpp>
#include <sobol_engine.hpp>

#include <chrono>
#include <vector>

auto main() -> int {
    {
        static_assert(adhoc::sobol<>::min() == 1U);
        static_assert(adhoc::sobol<>::max() == 4294967295U);
    }

    // first dimension in sobol is centered every 2^n-1 simulations
    {
        std::size_t paths = 1000;
        adhoc::open<adhoc::sobol<>> gen(1U);
        double tolerance = 1e-15;
        double acc = 0;
        std::size_t acc_n = 0;

        for (std::size_t dim = 0; dim <= paths; dim++) {
            acc += gen();
            acc_n++;
            std::size_t simulated_paths_p1 = dim + 2;
            auto mean = acc / static_cast<double>(acc_n);
            if ((simulated_paths_p1 & (simulated_paths_p1 - 1)) == 0) {
                EXPECT_NEAR_ABS(mean, 0.5, 1e-15);
            } else {
                bool isfar =
                    mean > (0.5 + tolerance) || mean < (0.5 - tolerance);
                EXPECT_EQUAL(isfar, true);
            }
        }
    }

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
        for (unsigned int i = 1; i < ndims; i++) {
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
        for (unsigned int i = 1; i < ndims; i++) {
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

        unsigned int ndims = 15;
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
        auto check_first = []<adhoc::qrng_table T>() {
            std::size_t ndims = 200;
            for (unsigned int i = 1; i < ndims; i++) {
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
        auto check_first = []<adhoc::qrng_table T>() {
            std::size_t ndims = 200;
            for (unsigned int i = 1; i < ndims; i++) {
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
        auto check_discard = []<adhoc::qrng_table T, bool First>() {
            std::size_t ndims = 30;
            for (unsigned int i = 1; i < ndims; i++) {
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
            []<adhoc::qrng_table T, bool First, std::size_t w = 32>() {
                std::size_t ndims = 200;
                for (unsigned int i = 1; i < ndims; i++) {
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
            []<adhoc::qrng_table T, bool First, std::size_t w = 32>() {
                std::size_t ndims = 200;
                for (unsigned int i = 180; i < ndims; i++) {
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

        auto check_uniform_vals = []<adhoc::qrng_table T, bool First,
                                     class UIntType, std::size_t w>() {
            unsigned int ndims = 20;
            adhoc::open<
                adhoc::sobol_engine<UIntType, w, true,
                                    adhoc::qrng_table::joe_kuo_old_1111>,
                float>
                rng(ndims);

            adhoc::open<
                adhoc::sobol_engine<UIntType, 32, true,
                                    adhoc::qrng_table::joe_kuo_old_1111>,
                float>
                rng32(ndims);

            for (std::size_t i = 0; i < 1000; i++) {
                auto val1 = rng32();
                auto val2 = rng();
                EXPECT_EQUAL(val1, val2);
            }
        };

        auto check_uniform_vals_mult =
            [&check_uniform_vals]<adhoc::qrng_table T, bool First>() {
                check_uniform_vals
                    .template operator()<T, First, std::uint64_t, 34>();
                check_uniform_vals
                    .template operator()<T, First, std::uint64_t, 33>();

                check_uniform_vals
                    .template operator()<T, First, std::uint32_t, 31>();
                check_uniform_vals
                    .template operator()<T, First, std::uint32_t, 30>();
                check_uniform_vals
                    .template operator()<T, First, std::uint32_t, 29>();
                check_uniform_vals
                    .template operator()<T, First, std::uint32_t, 28>();
                check_uniform_vals
                    .template operator()<T, First, std::uint32_t, 27>();
                check_uniform_vals
                    .template operator()<T, First, std::uint32_t, 26>();
                check_uniform_vals
                    .template operator()<T, First, std::uint32_t, 25>();
                check_uniform_vals
                    .template operator()<T, First, std::uint32_t, 24>();
                check_uniform_vals
                    .template operator()<T, First, std::uint32_t, 23>();
                check_uniform_vals
                    .template operator()<T, First, std::uint32_t, 22>();
                check_uniform_vals
                    .template operator()<T, First, std::uint32_t, 21>();
                check_uniform_vals
                    .template operator()<T, First, std::uint32_t, 20>();
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

    // discard for w>32
    {
        unsigned long long z = 0;
        z--;
        unsigned int dim = 20;
        adhoc::sobol_engine<std::uint64_t, 64, true,
                            adhoc::qrng_table::new_joe_kuo_5_21201>
            rng1(dim);
        adhoc::sobol_engine<std::uint64_t, 64, true,
                            adhoc::qrng_table::new_joe_kuo_5_21201>
            rng2(dim);

        for (std::size_t i = 0; i < dim - 1; i++) {
            rng1();
            rng1.discard(z);
        }

        rng1();
        rng1.discard(z - 10);

        for (std::size_t i = 0; i < 10; i++) {
            rng1();
        }

        EXPECT_EQUAL(rng1, rng2);
    }

    // fwd and back
    {

        auto check_bwd = []<adhoc::qrng_table T, bool First, class UIntType,
                            std::size_t w>() {
            for (unsigned int i = 1; i < 20; i++) {
                unsigned int n_dims = i;
                adhoc::sobol_engine<UIntType, w, true,
                                    adhoc::qrng_table::joe_kuo_old_1111>
                    rng(n_dims);

                adhoc::sobol_engine<UIntType, w, true,
                                    adhoc::qrng_table::joe_kuo_old_1111>
                    rng2(n_dims);

                check_fwd_and_back(rng, 1000);
                EXPECT_EQUAL(rng, rng2);

                check_back_and_fwd(rng, 1000);
                EXPECT_EQUAL(rng, rng2);
            }
        };

        auto check_bwd_mult = [&check_bwd]<adhoc::qrng_table T, bool First>() {
            check_bwd.template operator()<T, First, std::uint64_t, 64>();
            check_bwd.template operator()<T, First, std::uint64_t, 52>();
            check_bwd.template operator()<T, First, std::uint32_t, 32>();
            check_bwd.template operator()<T, First, std::uint32_t, 31>();
            check_bwd.template operator()<T, First, std::uint32_t, 27>();
        };

        check_bwd_mult
            .template operator()<adhoc::qrng_table::joe_kuo_old_1111, false>();
        check_bwd_mult.template
        operator()<adhoc::qrng_table::joe_kuo_other_0_7600, false>();
        check_bwd_mult.template
        operator()<adhoc::qrng_table::joe_kuo_other_2_3900, false>();
        check_bwd_mult.template
        operator()<adhoc::qrng_table::joe_kuo_other_3_7300, false>();
        check_bwd_mult.template
        operator()<adhoc::qrng_table::joe_kuo_other_4_5600, false>();
        check_bwd_mult.template
        operator()<adhoc::qrng_table::new_joe_kuo_5_21201, false>();
        check_bwd_mult.template
        operator()<adhoc::qrng_table::new_joe_kuo_6_21201, false>();
        check_bwd_mult.template
        operator()<adhoc::qrng_table::new_joe_kuo_7_21201, false>();

        check_bwd_mult
            .template operator()<adhoc::qrng_table::joe_kuo_old_1111, true>();
        check_bwd_mult.template
        operator()<adhoc::qrng_table::joe_kuo_other_0_7600, true>();
        check_bwd_mult.template
        operator()<adhoc::qrng_table::joe_kuo_other_2_3900, true>();
        check_bwd_mult.template
        operator()<adhoc::qrng_table::joe_kuo_other_3_7300, true>();
        check_bwd_mult.template
        operator()<adhoc::qrng_table::joe_kuo_other_4_5600, true>();
        check_bwd_mult.template
        operator()<adhoc::qrng_table::new_joe_kuo_5_21201, true>();
        check_bwd_mult.template
        operator()<adhoc::qrng_table::new_joe_kuo_6_21201, true>();
        check_bwd_mult.template
        operator()<adhoc::qrng_table::new_joe_kuo_7_21201, true>();
    }

    std::size_t sims = 0;
    if (auto env_p = std::getenv("TIMING_SIMS")) {
        sims = std::atoi(env_p);
    }

    if (sims) {

        // int reps = 1;
        // if (auto env_p = std::getenv("REP_SIMS")) {
        //     reps = std::atoi(env_p);
        // }

        // we have to repeat multiple times because times are slow first
        // for (std::size_t rep = 0; rep < reps; rep++) {

        {
            std::uint64_t res1 = 0;
            adhoc::sobol<false> rng;
            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims; i++) {
                res1 += rng();
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();
            std::cout << "time not skipping first" << std::endl;
            std::cout << time << std::endl;
            std::cout << res1 << std::endl;
        }

        {
            std::uint64_t res1 = 0;
            adhoc::sobol<> rng;
            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims; i++) {
                res1 += rng();
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();
            std::cout << "time skipping first" << std::endl;
            std::cout << time << std::endl;
            std::cout << res1 << std::endl;
        }

        {
            std::uint64_t res1 = 0;
            adhoc::sobol<false> rng;
            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims; i++) {
                res1 += rng();
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();
            std::cout << "time not skipping first" << std::endl;
            std::cout << time << std::endl;
            std::cout << res1 << std::endl;
        }

        {
            std::uint64_t res1 = 0;
            adhoc::sobol<> rng;
            auto time1 = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i < sims; i++) {
                res1 += rng();
            }
            auto time2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                            time2 - time1)
                            .count();
            std::cout << "time skipping first" << std::endl;
            std::cout << time << std::endl;
            std::cout << res1 << std::endl;
        }
        // }
    }

    TEST_END;
}
