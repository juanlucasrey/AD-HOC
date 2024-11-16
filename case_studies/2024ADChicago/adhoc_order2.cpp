#include <adhoc.hpp>
#include <backpropagator.hpp>
#include <calc_tree.hpp>
#include <differential_operator.hpp>

#include "black_scholes.hpp"

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

    using namespace adhoc;
    ADHOC(S);
    ADHOC(K);
    ADHOC(v);
    ADHOC(T);

    auto res = call_price(S, K, v, T);

    CalcTree ct(res);

    // order 1
    auto dS = d(S);
    auto dK = d(K);
    auto dv = d(v);
    auto dT = d(T);

    // order 2
    auto dSS = d<2>(S);
    auto dSK = d(S) * d(K);
    auto dSv = d(S) * d(v);
    auto dST = d(S) * d(T);
    auto dKK = d<2>(K);
    auto dKv = d(K) * d(v);
    auto dKT = d(K) * d(T);
    auto dvv = d<2>(v);
    auto dvT = d(v) * d(T);
    auto dTT = d<2>(T);

    auto dres = d(res);

    BackPropagator t(dS, dK, dv, dT, dSS, dSK, dSv, dST, dKK, dKv, dKT, dvv,
                     dvT, dTT, dres);

    auto time1 = std::chrono::high_resolution_clock::now();
    for (std::size_t j = 0; j < iters; ++j) {
        ct.set(S) = stock_distr(generator);
        ct.set(K) = stock_distr(generator);
        ct.set(v) = vol_distr(generator);
        ct.set(T) = time_distr(generator);
        ct.evaluate();
        t.set(dres) = 1.;
        t.backpropagate(ct);

        // average values in a single Taylor expansion
        results_average[0] += ct.get(res);
    }

    results_average[1] = t.get(dS);
    results_average[2] = t.get(dK);
    results_average[3] = t.get(dv);
    results_average[4] = t.get(dT);
    results_average[5] = t.get(dSS);
    results_average[6] = t.get(dSK);
    results_average[7] = t.get(dSv);
    results_average[8] = t.get(dST);
    results_average[9] = t.get(dKK);
    results_average[10] = t.get(dKv);
    results_average[11] = t.get(dKT);
    results_average[12] = t.get(dvv);
    results_average[13] = t.get(dvT);
    results_average[14] = t.get(dTT);

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
    std::cout << "AD-HOC order 2 time (ms): " << time << std::endl;

    return 0;
}
