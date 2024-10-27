#include "../../public4/adhoc.hpp"
#include "../../public4/backpropagator.hpp"
#include "../../public4/calc_tree.hpp"
#include "../../public4/differential_operator.hpp"

#include "black_scholes.hpp"

#include <chrono>
#include <iostream>
#include <random>
#include <tuple>

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

    std::array<double, 35> results_average;
    results_average.fill(0);

    using namespace adhoc4;
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
    auto dKK = d<2>(K);
    auto dvv = d<2>(v);
    auto dTT = d<2>(T);
    auto dSK = d(S) * d(K);
    auto dSv = d(S) * d(v);
    auto dST = d(S) * d(T);
    auto dKv = d(K) * d(v);
    auto dKT = d(K) * d(T);
    auto dvT = d(v) * d(T);

    // order 3
    auto dSSS = d<3>(S);
    auto dKKK = d<3>(K);
    auto dvvv = d<3>(v);
    auto dTTT = d<3>(T);

    auto dSSK = d<2>(S) * d(K);
    auto dSSv = d<2>(S) * d(v);
    auto dSST = d<2>(S) * d(T);
    auto dSKK = d(S) * d<2>(K);
    auto dSvv = d(S) * d<2>(v);
    auto dSTT = d(S) * d<2>(T);
    auto dKKv = d<2>(K) * d(v);
    auto dKKT = d<2>(K) * d(T);
    auto dKvv = d(K) * d<2>(v);
    auto dKTT = d(K) * d<2>(T);
    auto dvvT = d<2>(v) * d(T);
    auto dvTT = d(v) * d<2>(T);

    auto dSKv = d(S) * d(K) * d(v);
    auto dSKT = d(S) * d(K) * d(T);
    auto dSvT = d(S) * d(v) * d(T);
    auto dKvT = d(K) * d(v) * d(T);

    auto dres = d(res);

    BackPropagator t(dS, dK, dv, dT, dSS, dKK, dvv, dTT, dSK, dSv, dST, dKv,
                     dKT, dvT, dSSS, dKKK, dvvv, dTTT, dSSK, dSSv, dSST, dSKK,
                     dSvv, dSTT, dKKv, dKKT, dKvv, dKTT, dvvT, dvTT, dSKv, dSKT,
                     dSvT, dKvT, dres);

    time1 = std::chrono::high_resolution_clock::now();

    for (std::size_t j = 0; j < iters; ++j) {
        ct.set(S) = stock_distr(generator);
        ct.set(K) = stock_distr(generator);
        ct.set(v) = vol_distr(generator);
        ct.set(T) = time_distr(generator);
        ct.evaluate();
        t.reset();
        t.set(dres) = 1.;
        // t.backpropagate2<71>(ct);
        t.backpropagate(ct);

        // average to make sure compiler doesn't optimise calculations away
        results_average[0] += ct.get(res);
        results_average[1] += t.get(dS);
        results_average[2] += t.get(dK);
        results_average[3] += t.get(dv);
        results_average[4] += t.get(dT);
        results_average[5] += t.get(dSS);
        results_average[6] += t.get(dKK);
        results_average[7] += t.get(dvv);
        results_average[8] += t.get(dTT);
        results_average[9] += t.get(dSK);
        results_average[10] += t.get(dSv);
        results_average[11] += t.get(dST);
        results_average[12] += t.get(dKv);
        results_average[13] += t.get(dKT);
        results_average[14] += t.get(dvT);

        results_average[15] += t.get(dSSS);
        results_average[16] += t.get(dKKK);
        results_average[17] += t.get(dvvv);
        results_average[18] += t.get(dTTT);
        results_average[19] += t.get(dSSK);
        results_average[20] += t.get(dSSv);
        results_average[21] += t.get(dSST);
        results_average[22] += t.get(dSKK);
        results_average[23] += t.get(dSvv);
        results_average[24] += t.get(dSTT);
        results_average[25] += t.get(dKKv);
        results_average[26] += t.get(dKKT);
        results_average[27] += t.get(dKvv);
        results_average[28] += t.get(dKTT);
        results_average[29] += t.get(dvvT);
        results_average[30] += t.get(dvTT);
        results_average[31] += t.get(dSKv);
        results_average[32] += t.get(dSKT);
        results_average[33] += t.get(dSvT);
        results_average[34] += t.get(dKvT);
    }

    time2 = std::chrono::high_resolution_clock::now();
    auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1)
            .count();

    std::cout << "iterations: " << iters << std::endl;
    std::cout << "AD-HOC order 3 time (ms): " << time << std::endl;

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    for (std::size_t i = 0; i < results_average.size(); ++i) {
        results_average[i] /= static_cast<double>(iters);
        std::cout << results_average[i] << std::endl;
    }

    return 0;
}
