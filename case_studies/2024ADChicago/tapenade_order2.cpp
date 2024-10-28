#include "black_scholes_c.h"
#include "tapenade/c_version_bs_d_b.h"

#include <array>
#include <chrono>
#include <iostream>
#include <random>

int main() {
    std::mt19937 generator(123);
    std::uniform_real_distribution<double> stock_distr(90, 110.0);
    std::uniform_real_distribution<double> vol_distr(0.05, 0.3);
    std::uniform_real_distribution<double> time_distr(0.5, 1.5);

    using namespace std::chrono;
    std::chrono::time_point<std::chrono::high_resolution_clock> time1, time2;

    std::size_t iters = 1;
    if (auto env_p = std::getenv("ITERATIONS")) {
        iters = std::stoul(env_p);
    }

    std::array<double, 15> results_average;
    results_average.fill(0);

    double S0, K0, v0, T0, Sb = 0, Kb = 0, vb = 0, Tb = 0;

    double Sdb = 0, Kdb = 0, vdb = 0, Tdb = 0;
    double Sd = 0, Kd = 0, vd = 0, Td = 0;
    double call_price_val = 0.;
    double call_priceb = 0.0;
    double call_price_db = 1.;
    time1 = std::chrono::high_resolution_clock::now();

    for (std::size_t j = 0; j < iters; ++j) {
        S0 = stock_distr(generator);
        K0 = stock_distr(generator);
        v0 = vol_distr(generator);
        T0 = time_distr(generator);
        results_average[0] += call_price(S0, K0, v0, T0);

        // 1
        {
            Sd = 1.;
            call_price_d_b(S0, &Sb, Sd, &Sdb, K0, &Kb, Kd, &Kdb, v0, &vb, vd,
                           &vdb, T0, &Tb, Td, &Tdb, &call_price_val,
                           &call_priceb, call_price_db);

            results_average[1] += Sdb;
            results_average[2] += Kdb;
            results_average[3] += vdb;
            results_average[4] += Tdb;

            results_average[5] += Sb;
            results_average[6] += Kb;
            results_average[7] += vb;
            results_average[8] += Tb;

            Sd = 0.;

            Sb = 0;
            Kb = 0;
            vb = 0;
            Tb = 0;
            Sdb = 0;
            Kdb = 0;
            vdb = 0;
            Tdb = 0;
        }

        // 2
        {
            Kd = 1.;
            call_price_d_b(S0, &Sb, Sd, &Sdb, K0, &Kb, Kd, &Kdb, v0, &vb, vd,
                           &vdb, T0, &Tb, Td, &Tdb, &call_price_val,
                           &call_priceb, call_price_db);

            results_average[9] += Kb;
            results_average[10] += vb;
            results_average[11] += Tb;

            Kd = 0.;

            Sb = 0;
            Kb = 0;
            vb = 0;
            Tb = 0;
            Sdb = 0;
            Kdb = 0;
            vdb = 0;
            Tdb = 0;
        }

        // 3
        {
            vd = 1.;
            call_price_d_b(S0, &Sb, Sd, &Sdb, K0, &Kb, Kd, &Kdb, v0, &vb, vd,
                           &vdb, T0, &Tb, Td, &Tdb, &call_price_val,
                           &call_priceb, call_price_db);

            results_average[12] += vb;
            results_average[13] += Tb;

            vd = 0.;

            Sb = 0;
            Kb = 0;
            vb = 0;
            Tb = 0;
            Sdb = 0;
            Kdb = 0;
            vdb = 0;
            Tdb = 0;
        }

        // 4
        {
            Td = 1.;
            call_price_d_b(S0, &Sb, Sd, &Sdb, K0, &Kb, Kd, &Kdb, v0, &vb, vd,
                           &vdb, T0, &Tb, Td, &Tdb, &call_price_val,
                           &call_priceb, call_price_db);

            results_average[14] += Tb;

            Td = 0.;

            Sb = 0;
            Kb = 0;
            vb = 0;
            Tb = 0;
            Sdb = 0;
            Kdb = 0;
            vdb = 0;
            Tdb = 0;
        }
    }

    time2 = std::chrono::high_resolution_clock::now();
    auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1)
            .count();

    std::cout << "iterations: " << iters << std::endl;
    std::cout << "Tapenade order 2 time (ms): " << time << std::endl;

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    for (std::size_t i = 0; i < results_average.size(); ++i) {
        results_average[i] /= static_cast<double>(iters);
        std::cout << results_average[i] << std::endl;
    }

    return 0;
}
