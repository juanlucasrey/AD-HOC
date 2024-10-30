#include "black_scholes_c.h"
#include "tapenade/c_version_bs_b_d.h"

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

    double Sdb = 0, Kdb = 0, vdb = 0, Tdb = 0;
    double Sd = 0, Kd = 0, vd = 0, Td = 0;
    double call_price_val = 0.;
    double call_priceb = 0.0;
    double call_price_db = 1.;

    double Sbd = 0, Kbd = 0, vbd = 0, Tbd = 0;
    double call_pricebd = 0;
    auto time1 = std::chrono::high_resolution_clock::now();

    for (std::size_t j = 0; j < iters; ++j) {
        S0 = stock_distr(generator);
        K0 = stock_distr(generator);
        v0 = vol_distr(generator);
        T0 = time_distr(generator);
        results_average[0] += call_price(S0, K0, v0, T0);

        // 1
        {
            Sd = 1.;
            call_priceb = 1.0;
            call_price_b_d(S0, Sd, &Sb, &Sbd, K0, Kd, &Kb, &Kbd, v0, vd, &vb,
                           &vbd, T0, Td, &Tb, &Tbd, call_priceb, call_pricebd);

            results_average[1] += Sb;
            results_average[2] += Kb;
            results_average[3] += vb;
            results_average[4] += Tb;

            results_average[5] += Sbd;
            results_average[6] += Kbd;
            results_average[7] += vbd;
            results_average[8] += Tbd;

            Sd = 0.;

            Sb = 0;
            Kb = 0;
            vb = 0;
            Tb = 0;
            Sbd = 0;
            Kbd = 0;
            vbd = 0;
            Tbd = 0;
        }

        // 2
        {
            Kd = 1.;
            call_priceb = 1.0;
            call_price_b_d(S0, Sd, &Sb, &Sbd, K0, Kd, &Kb, &Kbd, v0, vd, &vb,
                           &vbd, T0, Td, &Tb, &Tbd, call_priceb, call_pricebd);

            results_average[9] += Kbd;
            results_average[10] += vbd;
            results_average[11] += Tbd;

            Kd = 0.;

            Sb = 0;
            Kb = 0;
            vb = 0;
            Tb = 0;
            Sbd = 0;
            Kbd = 0;
            vbd = 0;
            Tbd = 0;
        }

        // 3
        {
            vd = 1.;
            call_priceb = 1.0;
            call_price_b_d(S0, Sd, &Sb, &Sbd, K0, Kd, &Kb, &Kbd, v0, vd, &vb,
                           &vbd, T0, Td, &Tb, &Tbd, call_priceb, call_pricebd);

            results_average[12] += vbd;
            results_average[13] += Tbd;

            vd = 0.;

            Sb = 0;
            Kb = 0;
            vb = 0;
            Tb = 0;
            Sbd = 0;
            Kbd = 0;
            vbd = 0;
            Tbd = 0;
        }

        // 4
        {
            Td = 1.;
            call_priceb = 1.0;
            call_price_b_d(S0, Sd, &Sb, &Sbd, K0, Kd, &Kb, &Kbd, v0, vd, &vb,
                           &vbd, T0, Td, &Tb, &Tbd, call_priceb, call_pricebd);

            results_average[14] += Tbd;

            Td = 0.;

            Sb = 0;
            Kb = 0;
            vb = 0;
            Tb = 0;
            Sbd = 0;
            Kbd = 0;
            vbd = 0;
            Tbd = 0;
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
