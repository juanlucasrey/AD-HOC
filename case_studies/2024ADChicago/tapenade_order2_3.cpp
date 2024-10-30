#include "black_scholes_c.h"
#include "tapenade/c_version_bs_b_b.h"

#include <array>
#include <chrono>
#include <iostream>
#include <random>

int main() {
    std::mt19937 generator(123);
    std::uniform_real_distribution<double> stock_distr(90, 110.0),
        vol_distr(0.05, 0.3), time_distr(0.5, 1.5);

    std::size_t iters = 1;
    if (auto env_p = std::getenv("ITERATIONS")) {
        iters = std::stoul(env_p);
    }

    std::array<double, 15> results_average{};

    double S0, K0, v0, T0, Sb = 0, Kb = 0, vb = 0, Tb = 0;

    double call_priceb = 0.0;
    double Sb0 = 0, Kb0 = 0, vb0 = 0, Tb0 = 0;
    double call_pricebb;
    double Sbb = 0, Kbb = 0, vbb = 0, Tbb = 0;
    auto time1 = std::chrono::high_resolution_clock::now();

    for (std::size_t j = 0; j < iters; ++j) {
        S0 = stock_distr(generator);
        K0 = stock_distr(generator);
        v0 = vol_distr(generator);
        T0 = time_distr(generator);
        results_average[0] += call_price(S0, K0, v0, T0);

        // 1
        {
            call_priceb = 1;
            Sb = 1;
            Sbb = 1;
            call_price_b_b(S0, &Sb0, &Sb, &Sbb, K0, &Kb0, &Kb, &Kbb, v0, &vb0,
                           &vb, &vbb, T0, &Tb0, &Tb, &Tbb, call_priceb,
                           &call_pricebb);

            results_average[1] += call_pricebb;

            results_average[5] += Sb0;
            results_average[6] += Kb0;
            results_average[7] += vb0;
            results_average[8] += Tb0;

            call_pricebb = 0;
            Sb = 0;
            Kb = 0;
            vb = 0;
            Tb = 0;
            Sbb = 0;
            Kbb = 0;
            vbb = 0;
            Tbb = 0;
            Sb0 = 0;
            Kb0 = 0;
            vb0 = 0;
            Tb0 = 0;
        }

        // 2
        {
            call_priceb = 1;
            Kb = 1;
            Kbb = 1;
            call_price_b_b(S0, &Sb0, &Sb, &Sbb, K0, &Kb0, &Kb, &Kbb, v0, &vb0,
                           &vb, &vbb, T0, &Tb0, &Tb, &Tbb, call_priceb,
                           &call_pricebb);

            results_average[2] += call_pricebb;

            results_average[9] += Kb0;
            results_average[10] += vb0;
            results_average[11] += Tb0;

            call_pricebb = 0;
            Sb = 0;
            Kb = 0;
            vb = 0;
            Tb = 0;
            Sbb = 0;
            Kbb = 0;
            vbb = 0;
            Tbb = 0;
            Sb0 = 0;
            Kb0 = 0;
            vb0 = 0;
            Tb0 = 0;
        }

        // 3
        {
            call_priceb = 1;
            vb = 1;
            vbb = 1;
            call_price_b_b(S0, &Sb0, &Sb, &Sbb, K0, &Kb0, &Kb, &Kbb, v0, &vb0,
                           &vb, &vbb, T0, &Tb0, &Tb, &Tbb, call_priceb,
                           &call_pricebb);

            results_average[3] += call_pricebb;

            results_average[12] += vb0;
            results_average[13] += Tb0;

            call_pricebb = 0;
            Sb = 0;
            Kb = 0;
            vb = 0;
            Tb = 0;
            Sbb = 0;
            Kbb = 0;
            vbb = 0;
            Tbb = 0;
            Sb0 = 0;
            Kb0 = 0;
            vb0 = 0;
            Tb0 = 0;
        }

        // 4
        {
            call_priceb = 1;
            Tb = 1;
            Tbb = 1;
            call_price_b_b(S0, &Sb0, &Sb, &Sbb, K0, &Kb0, &Kb, &Kbb, v0, &vb0,
                           &vb, &vbb, T0, &Tb0, &Tb, &Tbb, call_priceb,
                           &call_pricebb);

            results_average[4] += call_pricebb;

            results_average[14] += Tb0;

            call_pricebb = 0;
            Sb = 0;
            Kb = 0;
            vb = 0;
            Tb = 0;
            Sbb = 0;
            Kbb = 0;
            vbb = 0;
            Tbb = 0;
            Sb0 = 0;
            Kb0 = 0;
            vb0 = 0;
            Tb0 = 0;
        }
    }

    auto time2 = std::chrono::high_resolution_clock::now();
    auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1)
            .count();

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    for (std::size_t i = 0; i < results_average.size(); ++i) {
        results_average[i] /= static_cast<double>(iters);
        std::cout << results_average[i] << std::endl;
    }

    std::cout << "iterations: " << iters << std::endl;
    std::cout << "Tapenade order 2 time (ms): " << time << std::endl;

    return 0;
}
