
#include <adhoc.hpp>
#include <base.hpp>

#include <functional>

#include "statistics.hpp"
#include "test_tools.hpp"

#define TIMING_TEST false

#if TIMING_TEST
#include <chrono>
#include <random>
#endif

namespace {

auto zeta_borwein_precision(double x, std::size_t i) {
    using namespace adhoc4;
    if (i == 3) {
        return detail::zeta_borwein<3>(x);
    }

    if (i == 4) {
        return detail::zeta_borwein<4>(x);
    }

    if (i == 5) {
        return detail::zeta_borwein<5>(x);
    }

    if (i == 6) {
        return detail::zeta_borwein<6>(x);
    }

    if (i == 7) {
        return detail::zeta_borwein<7>(x);
    }

    if (i == 8) {
        return detail::zeta_borwein<8>(x);
    }

    if (i == 9) {
        return detail::zeta_borwein<9>(x);
    }

    if (i == 10) {
        return detail::zeta_borwein<10>(x);
    }

    if (i == 11) {
        return detail::zeta_borwein<11>(x);
    }

    if (i == 12) {
        return detail::zeta_borwein<12>(x);
    }

    if (i == 13) {
        return detail::zeta_borwein<13>(x);
    }

    if (i == 14) {
        return detail::zeta_borwein<14>(x);
    }

    if (i == 15) {
        return detail::zeta_borwein<15>(x);
    }

    if (i == 16) {
        return detail::zeta_borwein<16>(x);
    }

    if (i == 17) {
        return detail::zeta_borwein<17>(x);
    }

    if (i == 18) {
        return detail::zeta_borwein<18>(x);
    }

    if (i == 19) {
        return detail::zeta_borwein<19>(x);
    }

    if (i == 20) {
        return detail::zeta_borwein<20>(x);
    }

    if (i == 21) {
        return detail::zeta_borwein<21>(x);
    }

    if (i == 22) {
        return detail::zeta_borwein<22>(x);
    }

    if (i == 23) {
        return detail::zeta_borwein<23>(x);
    }

    if (i == 24) {
        return detail::zeta_borwein<24>(x);
    }

    return detail::zeta_borwein<25>(x);
}

} // namespace

int main() {
    TEST_START;

    using namespace adhoc4;
    std::cout.precision(std::numeric_limits<double>::max_digits10);

    // ZetaPrecision
    {
        std::string name = "data/zeta.csv";
        std::vector<double> in;
        std::vector<std::array<double, 7>> out;

        get_values(name, in, out);

        std::vector<std::array<double, 3>> borwein_errors{
            {1e-2, 1e-2, 1e-1},    {1e-4, 1e-3, 1e-2},    {1e-4, 1e-4, 1e-3},
            {1e-6, 1e-5, 1e-4},    {1e-6, 1e-6, 1e-5},    {1e-8, 1e-7, 1e-5},
            {1e-7, 1e-7, 1e-6},    {1e-9, 1e-8, 1e-7},    {1e-9, 1e-9, 1e-8},
            {1e-10, 1e-10, 1e-8},  {1e-11, 1e-10, 1e-9},  {1e-12, 1e-11, 1e-10},
            {1e-12, 1e-12, 1e-11}, {1e-13, 1e-13, 1e-11}, {1e-14, 1e-14, 1e-12},
            {1e-15, 1e-14, 1e-13}, {1e-15, 1e-15, 1e-14}, {5e-16, 1e-15, 1e-14},
            {5e-16, 1e-15, 5e-14}, {1e-15, 1e-15, 1e-13}, {1e-16, 1e-15, 1e-14},
            {1e-15, 1e-15, 1e-14}};

        for (std::size_t i = 3; i < 25; i++) {
            std::function<double(double)> function = [i](double d) {
                return zeta_borwein_precision(d, i);
            };
            auto stats = statistics(in, out, function);
            EXPECT_LESS_THAN(std::get<0>(stats), borwein_errors[i - 3][0]);
            EXPECT_LESS_THAN(std::get<1>(stats), borwein_errors[i - 3][1]);
            EXPECT_LESS_THAN(std::get<2>(stats), borwein_errors[i - 3][2]);
        }

        std::function<double(double)> functionstd = [](double d) {
            return std::riemann_zeta(d);
        };

        auto stats = statistics(in, out, functionstd);
        EXPECT_LESS_THAN(std::get<0>(stats), 1e-15);
        EXPECT_LESS_THAN(std::get<1>(stats), 1e-14);
        EXPECT_LESS_THAN(std::get<2>(stats), 1e-13);

        std::vector<std::array<double, 3>> dererrors{{1e-15, 1e-14, 1e-12},
                                                     {1e-13, 1e-13, 1e-10},
                                                     {1e-11, 1e-11, 1e-8},
                                                     {1e-9, 1e-9, 1e-6}};

        for (std::size_t i = 1; i < 5; i++) {

            std::function<double(double)> function = [i](double d) {
                std::array<double, 10> results1;
                results1.fill(0.);
                riemann_zeta_t<double>::d<6>(d, d, results1);
                return results1[i - 1];
            };

            auto stats = statistics(in, out, function, i);
            EXPECT_LESS_THAN(std::get<0>(stats), dererrors[i - 1][0]);
            EXPECT_LESS_THAN(std::get<1>(stats), dererrors[i - 1][1]);
            EXPECT_LESS_THAN(std::get<2>(stats), dererrors[i - 1][2]);
        }
    }

    // ZetaPrecisionLargePositive
    {
        std::string name = "data/zeta_pos_large.csv";
        std::vector<double> in;
        std::vector<std::array<double, 7>> out;

        get_values(name, in, out);

        std::vector<std::array<double, 3>> borwein_errors{
            {1e-1, 1e-1, 1e-1},    {1e-2, 1e-2, 1e-2},    {1e-3, 1e-3, 1e-3},
            {1e-4, 1e-4, 1e-4},    {1e-5, 1e-5, 1e-5},    {1e-5, 1e-5, 1e-5},
            {1e-6, 1e-6, 1e-6},    {1e-7, 1e-7, 1e-7},    {1e-8, 1e-8, 1e-8},
            {1e-8, 1e-8, 1e-8},    {1e-9, 1e-9, 1e-9},    {1e-10, 1e-10, 1e-10},
            {1e-11, 1e-11, 1e-11}, {1e-11, 1e-11, 1e-11}, {1e-12, 1e-12, 1e-12},
            {1e-13, 1e-13, 1e-13}, {1e-14, 1e-14, 1e-14}, {1e-15, 1e-15, 1e-15},
            {1e-15, 1e-15, 1e-15}, {1e-16, 1e-16, 1e-16}, {1e-16, 1e-16, 1e-16},
            {1e-16, 1e-16, 1e-16}};

        for (std::size_t i = 3; i < 25; i++) {
            std::function<double(double)> function = [i](double d) {
                return zeta_borwein_precision(d, i);
            };
            auto stats = statistics(in, out, function);
            EXPECT_LESS_THAN(std::get<0>(stats), borwein_errors[i - 3][0]);
            EXPECT_LESS_THAN(std::get<1>(stats), borwein_errors[i - 3][1]);
            EXPECT_LESS_THAN(std::get<2>(stats), borwein_errors[i - 3][2]);
        }

        std::function<double(double)> functionstd = [](double d) {
            return std::riemann_zeta(d);
        };

        auto stats = statistics(in, out, functionstd);
        EXPECT_LESS_THAN(std::get<0>(stats), 1e-16);
        EXPECT_LESS_THAN(std::get<1>(stats), 1e-16);
        EXPECT_LESS_THAN(std::get<2>(stats), 1e-16);

        std::vector<std::array<double, 3>> dererrors{{1e-12, 1e-12, 1e-12},
                                                     {1e-12, 1e-12, 1e-12},
                                                     {1e-12, 1e-12, 1e-12},
                                                     {1e-12, 1e-12, 1e-12}};

        for (std::size_t i = 1; i < 5; i++) {
            std::function<double(double)> function = [i](double d) {
                std::array<double, 10> results1;
                results1.fill(0.);
                riemann_zeta_t<double>::d<6>(d, d, results1);
                return results1[i - 1];
            };

            auto stats = statistics(in, out, function, i);
            EXPECT_LESS_THAN(std::get<0>(stats), dererrors[i - 1][0]);
            EXPECT_LESS_THAN(std::get<1>(stats), dererrors[i - 1][1]);
            EXPECT_LESS_THAN(std::get<2>(stats), dererrors[i - 1][2]);
        }
    }

    // ZetaPrecisionLargeNegative
    {
        std::string name = "data/zeta_neg_large.csv";
        std::vector<double> in;
        std::vector<std::array<double, 7>> out;

        get_values(name, in, out);

        std::vector<std::array<double, 3>> borwein_errors{
            {1e-1, 1e-1, 1e-1},    {1e-2, 1e-2, 1e-2},    {1e-3, 1e-3, 1e-3},
            {1e-4, 1e-4, 1e-4},    {1e-5, 1e-5, 1e-5},    {1e-5, 1e-5, 1e-5},
            {1e-6, 1e-6, 1e-6},    {1e-7, 1e-7, 1e-7},    {1e-8, 1e-8, 1e-8},
            {1e-8, 1e-8, 1e-8},    {1e-9, 1e-9, 1e-9},    {1e-10, 1e-10, 1e-10},
            {1e-11, 1e-11, 1e-11}, {1e-11, 1e-11, 1e-11}, {1e-12, 1e-12, 1e-11},
            {1e-13, 1e-13, 1e-11}, {1e-14, 1e-13, 1e-11}, {1e-14, 1e-13, 1e-11},
            {1e-14, 1e-13, 1e-11}, {1e-14, 1e-13, 1e-11}, {1e-14, 1e-13, 1e-11},
            {1e-14, 1e-13, 1e-11}};

        for (std::size_t i = 3; i < 25; i++) {
            std::function<double(double)> function = [i](double d) {
                return zeta_borwein_precision(d, i);
            };
            auto stats = statistics(in, out, function);
            EXPECT_LESS_THAN(std::get<0>(stats), borwein_errors[i - 3][0]);
            EXPECT_LESS_THAN(std::get<1>(stats), borwein_errors[i - 3][1]);
            EXPECT_LESS_THAN(std::get<2>(stats), borwein_errors[i - 3][2]);
        }

        std::function<double(double)> functionstd = [](double d) {
            return std::riemann_zeta(d);
        };

        auto stats = statistics(in, out, functionstd);
        EXPECT_LESS_THAN(std::get<0>(stats), 1e-14);
        EXPECT_LESS_THAN(std::get<1>(stats), 1e-13);
        EXPECT_LESS_THAN(std::get<2>(stats), 1e-11);

        std::vector<std::array<double, 3>> dererrors{{1e-14, 1e-14, 1e-13},
                                                     {1e-14, 1e-14, 1e-12},
                                                     {1e-13, 1e-13, 1e-10},
                                                     {1e-14, 1e-14, 1e-11}};

        for (std::size_t i = 1; i < 5; i++) {
            std::function<double(double)> function = [i](double d) {
                std::array<double, 10> results1;
                results1.fill(0.);
                riemann_zeta_t<double>::d<6>(d, d, results1);
                return results1[i - 1];
            };

            auto stats = statistics(in, out, function, i);
            EXPECT_LESS_THAN(std::get<0>(stats), dererrors[i - 1][0]);
            EXPECT_LESS_THAN(std::get<1>(stats), dererrors[i - 1][1]);
            EXPECT_LESS_THAN(std::get<2>(stats), dererrors[i - 1][2]);
        }
    }

#if TIMING_TEST
    // ZetaBorweinTiming
    auto starttime = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();

    {
        std::uniform_real_distribution<double> d(1.1, 10.0);
        std::size_t iters = 100000;

        {
            std::mt19937 g(123);
            double acc = 0.;
            starttime = std::chrono::steady_clock::now();
            for (std::size_t i = 0; i < iters; i++) {
                double randomin = d(g);
                acc += std::riemann_zeta(randomin);
            }
            end = std::chrono::steady_clock::now();
            std::cout << "std "
                      << std::chrono::duration_cast<std::chrono::microseconds>(
                             end - starttime)
                             .count()
                      << " µs" << std::endl;

            acc /= static_cast<double>(iters);
            std::cout << acc << std::endl;
        }

        {
            std::mt19937 g(123);
            double acc = 0.;
            starttime = std::chrono::steady_clock::now();
            for (std::size_t i = 0; i < iters; i++) {
                double randomin = d(g);
                acc += detail::zeta_borwein<20>(randomin);
            }
            end = std::chrono::steady_clock::now();
            std::cout << "borwein 20 "
                      << std::chrono::duration_cast<std::chrono::microseconds>(
                             end - starttime)
                             .count()
                      << " µs" << std::endl;

            acc /= static_cast<double>(iters);
            std::cout << acc << std::endl;
        }
    }

    // ZetaBorweinTimingUnderOne
    {
        std::uniform_real_distribution<double> d(-10, 0.9);
        std::size_t iters = 100000;

        {
            std::mt19937 g(123);
            double acc = 0.;
            starttime = std::chrono::steady_clock::now();
            for (std::size_t i = 0; i < iters; i++) {
                double randomin = d(g);
                acc += std::riemann_zeta(randomin);
            }
            end = std::chrono::steady_clock::now();
            std::cout << "std "
                      << std::chrono::duration_cast<std::chrono::microseconds>(
                             end - starttime)
                             .count()
                      << " µs" << std::endl;
            acc /= static_cast<double>(iters);
            std::cout << acc << std::endl;
        }

        {
            std::mt19937 g(123);
            double acc = 0.;
            starttime = std::chrono::steady_clock::now();
            for (std::size_t i = 0; i < iters; i++) {
                double randomin = d(g);
                acc += detail::zeta_borwein<20>(randomin);
            }
            end = std::chrono::steady_clock::now();
            std::cout << "borwein 20 "
                      << std::chrono::duration_cast<std::chrono::microseconds>(
                             end - starttime)
                             .count()
                      << " µs" << std::endl;
            acc /= static_cast<double>(iters);
            std::cout << acc << std::endl;
        }
    }
#endif // TIMING_TEST

    TEST_END;
}
