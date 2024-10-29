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

    std::array<double, 70> results_average{};

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
    auto dSK = d(S) * d(K);
    auto dSv = d(S) * d(v);
    auto dST = d(S) * d(T);
    auto dKK = d<2>(K);
    auto dKv = d(K) * d(v);
    auto dKT = d(K) * d(T);
    auto dvv = d<2>(v);
    auto dvT = d(v) * d(T);
    auto dTT = d<2>(T);

    // order 3
    auto dSSS = d<3>(S);
    auto dSSK = d<2>(S) * d(K);
    auto dSSv = d<2>(S) * d(v);
    auto dSST = d<2>(S) * d(T);
    auto dSKK = d(S) * d<2>(K);
    auto dSKv = d(S) * d(K) * d(v);
    auto dSKT = d(S) * d(K) * d(T);
    auto dSvv = d(S) * d<2>(v);
    auto dSvT = d(S) * d(v) * d(T);
    auto dSTT = d(S) * d<2>(T);
    auto dKKK = d<3>(K);
    auto dKKv = d<2>(K) * d(v);
    auto dKKT = d<2>(K) * d(T);
    auto dKvv = d(K) * d<2>(v);
    auto dKvT = d(K) * d(v) * d(T);
    auto dKTT = d(K) * d<2>(T);
    auto dvvv = d<3>(v);
    auto dvvT = d<2>(v) * d(T);
    auto dvTT = d(v) * d<2>(T);
    auto dTTT = d<3>(T);

    // order 4
    auto dSSSS = d<4>(S);
    auto dSSSK = d<3>(S) * d(K);
    auto dSSSv = d<3>(S) * d(v);
    auto dSSST = d<3>(S) * d(T);
    auto dSSKK = d<2>(S) * d<2>(K);
    auto dSSKv = d<2>(S) * d(K) * d(v);
    auto dSSKT = d<2>(S) * d(K) * d(T);
    auto dSSvv = d<2>(S) * d<2>(v);
    auto dSSvT = d<2>(S) * d(v) * d(T);
    auto dSSTT = d<2>(S) * d<2>(T);
    auto dSKKK = d(S) * d<3>(K);
    auto dSKKv = d(S) * d<2>(K) * d(v);
    auto dSKKT = d(S) * d<2>(K) * d(T);
    auto dSKvv = d(S) * d(K) * d<2>(v);
    auto dSKvT = d(S) * d(K) * d(v) * d(T);
    auto dSKTT = d(S) * d(K) * d<2>(T);
    auto dSvvv = d(S) * d<3>(v);
    auto dSvvT = d(S) * d<2>(v) * d(T);
    auto dSvTT = d(S) * d(v) * d<2>(T);
    auto dSTTT = d(S) * d<3>(T);
    auto dKKKK = d<4>(K);
    auto dKKKv = d<3>(K) * d(v);
    auto dKKKT = d<3>(K) * d(T);
    auto dKKvv = d<2>(K) * d<2>(v);
    auto dKKvT = d<2>(K) * d(v) * d(T);
    auto dKKTT = d<2>(K) * d<2>(T);
    auto dKvvv = d(K) * d<3>(v);
    auto dKvvT = d(K) * d<2>(v) * d(T);
    auto dKvTT = d(K) * d(v) * d<2>(T);
    auto dKTTT = d(K) * d<3>(T);
    auto dvvvv = d<4>(v);
    auto dvvvT = d<3>(v) * d(T);
    auto dvvTT = d<2>(v) * d<2>(T);
    auto dvTTT = d(v) * d<3>(T);
    auto dTTTT = d<4>(T);

    auto dres = d(res);

    BackPropagator t(dS, dK, dv, dT, dSS, dSK, dSv, dST, dKK, dKv, dKT, dvv,
                     dvT, dTT, dSSS, dSSK, dSSv, dSST, dSKK, dSKv, dSKT, dSvv,
                     dSvT, dSTT, dKKK, dKKv, dKKT, dKvv, dKvT, dKTT, dvvv, dvvT,
                     dvTT, dTTT, dSSSS, dSSSK, dSSSv, dSSST, dSSKK, dSSKv,
                     dSSKT, dSSvv, dSSvT, dSSTT, dSKKK, dSKKv, dSKKT, dSKvv,
                     dSKvT, dSKTT, dSvvv, dSvvT, dSvTT, dSTTT, dKKKK, dKKKv,
                     dKKKT, dKKvv, dKKvT, dKKTT, dKvvv, dKvvT, dKvTT, dKTTT,
                     dvvvv, dvvvT, dvvTT, dvTTT, dTTTT, dres);

    auto time1 = std::chrono::high_resolution_clock::now();
    for (std::size_t j = 0; j < iters; ++j) {
        ct.set(S) = stock_distr(generator);
        ct.set(K) = stock_distr(generator);
        ct.set(v) = vol_distr(generator);
        ct.set(T) = time_distr(generator);
        ct.evaluate();
        t.set(dres) = 1.;
        // t.backpropagate2<174>(ct);
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
    results_average[15] = t.get(dSSS);
    results_average[16] = t.get(dSSK);
    results_average[17] = t.get(dSSv);
    results_average[18] = t.get(dSST);
    results_average[19] = t.get(dSKK);
    results_average[20] = t.get(dSKv);
    results_average[21] = t.get(dSKT);
    results_average[22] = t.get(dSvv);
    results_average[23] = t.get(dSvT);
    results_average[24] = t.get(dSTT);
    results_average[25] = t.get(dKKK);
    results_average[26] = t.get(dKKv);
    results_average[27] = t.get(dKKT);
    results_average[28] = t.get(dKvv);
    results_average[29] = t.get(dKvT);
    results_average[30] = t.get(dKTT);
    results_average[31] = t.get(dvvv);
    results_average[32] = t.get(dvvT);
    results_average[33] = t.get(dvTT);
    results_average[34] = t.get(dTTT);
    results_average[35] = t.get(dSSSS);
    results_average[36] = t.get(dSSSK);
    results_average[37] = t.get(dSSSv);
    results_average[38] = t.get(dSSST);
    results_average[39] = t.get(dSSKK);
    results_average[40] = t.get(dSSKv);
    results_average[41] = t.get(dSSKT);
    results_average[42] = t.get(dSSvv);
    results_average[43] = t.get(dSSvT);
    results_average[44] = t.get(dSSTT);
    results_average[45] = t.get(dSKKK);
    results_average[46] = t.get(dSKKv);
    results_average[47] = t.get(dSKKT);
    results_average[48] = t.get(dSKvv);
    results_average[49] = t.get(dSKvT);
    results_average[50] = t.get(dSKTT);
    results_average[51] = t.get(dSvvv);
    results_average[52] = t.get(dSvvT);
    results_average[53] = t.get(dSvTT);
    results_average[54] = t.get(dSTTT);
    results_average[55] = t.get(dKKKK);
    results_average[56] = t.get(dKKKv);
    results_average[57] = t.get(dKKKT);
    results_average[58] = t.get(dKKvv);
    results_average[59] = t.get(dKKvT);
    results_average[60] = t.get(dKKTT);
    results_average[61] = t.get(dKvvv);
    results_average[62] = t.get(dKvvT);
    results_average[63] = t.get(dKvTT);
    results_average[64] = t.get(dKTTT);
    results_average[65] = t.get(dvvvv);
    results_average[66] = t.get(dvvvT);
    results_average[67] = t.get(dvvTT);
    results_average[68] = t.get(dvTTT);
    results_average[69] = t.get(dTTTT);

    auto time2 = std::chrono::high_resolution_clock::now();
    auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1)
            .count();

    std::cout << "iterations: " << iters << std::endl;
    std::cout << "AD-HOC order 4 time (ms): " << time << std::endl;

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    for (std::size_t i = 0; i < results_average.size(); ++i) {
        results_average[i] /= static_cast<double>(iters);
        std::cout << results_average[i] << std::endl;
    }

    return 0;
}
