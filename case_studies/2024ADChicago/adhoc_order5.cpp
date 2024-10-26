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

    std::array<double, 126> results_average;
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

    // order 5
    auto dSSSSS = d<5>(S);
    auto dSSSSK = d<4>(S) * d(K);
    auto dSSSSv = d<4>(S) * d(v);
    auto dSSSST = d<4>(S) * d(T);
    auto dSSSKK = d<3>(S) * d<2>(K);
    auto dSSSKv = d<3>(S) * d(K) * d(v);
    auto dSSSKT = d<3>(S) * d(K) * d(T);
    auto dSSSvv = d<3>(S) * d<2>(v);
    auto dSSSvT = d<3>(S) * d(v) * d(T);
    auto dSSSTT = d<3>(S) * d<2>(T);
    auto dSSKKK = d<2>(S) * d<3>(K);
    auto dSSKKv = d<2>(S) * d<2>(K) * d(v);
    auto dSSKKT = d<2>(S) * d<2>(K) * d(T);
    auto dSSKvv = d<2>(S) * d(K) * d<2>(v);
    auto dSSKvT = d<2>(S) * d(K) * d(v) * d(T);
    auto dSSKTT = d<2>(S) * d(K) * d<2>(T);
    auto dSSvvv = d<2>(S) * d<3>(v);
    auto dSSvvT = d<2>(S) * d<2>(v) * d(T);
    auto dSSvTT = d<2>(S) * d(v) * d<2>(T);
    auto dSSTTT = d<2>(S) * d<3>(T);
    auto dSKKKK = d(S) * d<4>(K);
    auto dSKKKv = d(S) * d<3>(K) * d(v);
    auto dSKKKT = d(S) * d<3>(K) * d(T);
    auto dSKKvv = d(S) * d<2>(K) * d<2>(v);
    auto dSKKvT = d(S) * d<2>(K) * d(v) * d(T);
    auto dSKKTT = d(S) * d<2>(K) * d<2>(T);
    auto dSKvvv = d(S) * d(K) * d<3>(v);
    auto dSKvvT = d(S) * d(K) * d<2>(v) * d(T);
    auto dSKvTT = d(S) * d(K) * d(v) * d<2>(T);
    auto dSKTTT = d(S) * d(K) * d<3>(T);
    auto dSvvvv = d(S) * d<4>(v);
    auto dSvvvT = d(S) * d<3>(v) * d(T);
    auto dSvvTT = d(S) * d<2>(v) * d<2>(T);
    auto dSvTTT = d(S) * d(v) * d<3>(T);
    auto dSTTTT = d(S) * d<4>(T);
    auto dKKKKK = d<5>(K);
    auto dKKKKv = d<4>(K) * d(v);
    auto dKKKKT = d<4>(K) * d(T);
    auto dKKKvv = d<3>(K) * d<2>(v);
    auto dKKKvT = d<3>(K) * d(v) * d(T);
    auto dKKKTT = d<3>(K) * d<2>(T);
    auto dKKvvv = d<2>(K) * d<3>(v);
    auto dKKvvT = d<2>(K) * d<2>(v) * d(T);
    auto dKKvTT = d<2>(K) * d(v) * d<2>(T);
    auto dKKTTT = d<2>(K) * d<3>(T);
    auto dKvvvv = d(K) * d<4>(v);
    auto dKvvvT = d(K) * d<3>(v) * d(T);
    auto dKvvTT = d(K) * d<2>(v) * d<2>(T);
    auto dKvTTT = d(K) * d(v) * d<3>(T);
    auto dKTTTT = d(K) * d<4>(T);
    auto dvvvvv = d<5>(v);
    auto dvvvvT = d<4>(v) * d(T);
    auto dvvvTT = d<3>(v) * d<2>(T);
    auto dvvTTT = d<2>(v) * d<3>(T);
    auto dvTTTT = d(v) * d<4>(T);
    auto dTTTTT = d<5>(T);

    auto dres = d(res);

    BackPropagator t(
        dS, dK, dv, dT, dSS, dKK, dvv, dTT, dSK, dSv, dST, dKv, dKT, dvT, dSSS,
        dKKK, dvvv, dTTT, dSSK, dSSv, dSST, dSKK, dSvv, dSTT, dKKv, dKKT, dKvv,
        dKTT, dvvT, dvTT, dSKv, dSKT, dSvT, dKvT, dSSSS, dSSSK, dSSSv, dSSST,
        dSSKK, dSSKv, dSSKT, dSSvv, dSSvT, dSSTT, dSKKK, dSKKv, dSKKT, dSKvv,
        dSKvT, dSKTT, dSvvv, dSvvT, dSvTT, dSTTT, dKKKK, dKKKv, dKKKT, dKKvv,
        dKKvT, dKKTT, dKvvv, dKvvT, dKvTT, dKTTT, dvvvv, dvvvT, dvvTT, dvTTT,
        dTTTT, dSSSSS, dSSSSK, dSSSSv, dSSSST, dSSSKK, dSSSKv, dSSSKT, dSSSvv,
        dSSSvT, dSSSTT, dSSKKK, dSSKKv, dSSKKT, dSSKvv, dSSKvT, dSSKTT, dSSvvv,
        dSSvvT, dSSvTT, dSSTTT, dSKKKK, dSKKKv, dSKKKT, dSKKvv, dSKKvT, dSKKTT,
        dSKvvv, dSKvvT, dSKvTT, dSKTTT, dSvvvv, dSvvvT, dSvvTT, dSvTTT, dSTTTT,
        dKKKKK, dKKKKv, dKKKKT, dKKKvv, dKKKvT, dKKKTT, dKKvvv, dKKvvT, dKKvTT,
        dKKTTT, dKvvvv, dKvvvT, dKvvTT, dKvTTT, dKTTTT, dvvvvv, dvvvvT, dvvvTT,
        dvvTTT, dvTTTT, dTTTTT, dres);

    time1 = std::chrono::high_resolution_clock::now();

    for (std::size_t j = 0; j < iters; ++j) {
        ct.set(S) = stock_distr(generator);
        ct.set(K) = stock_distr(generator);
        ct.set(v) = vol_distr(generator);
        ct.set(T) = time_distr(generator);
        ct.evaluate();
        t.reset();
        t.set(dres) = 1.;
        t.backpropagate2<368>(ct);
        // t.backpropagate(ct);

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
        results_average[35] += t.get(dSSSS);
        results_average[36] += t.get(dSSSK);
        results_average[37] += t.get(dSSSv);
        results_average[38] += t.get(dSSST);
        results_average[39] += t.get(dSSKK);
        results_average[40] += t.get(dSSKv);
        results_average[41] += t.get(dSSKT);
        results_average[42] += t.get(dSSvv);
        results_average[43] += t.get(dSSvT);
        results_average[44] += t.get(dSSTT);
        results_average[45] += t.get(dSKKK);
        results_average[46] += t.get(dSKKv);
        results_average[47] += t.get(dSKKT);
        results_average[48] += t.get(dSKvv);
        results_average[49] += t.get(dSKvT);
        results_average[50] += t.get(dSKTT);
        results_average[51] += t.get(dSvvv);
        results_average[52] += t.get(dSvvT);
        results_average[53] += t.get(dSvTT);
        results_average[54] += t.get(dSTTT);
        results_average[55] += t.get(dKKKK);
        results_average[56] += t.get(dKKKv);
        results_average[57] += t.get(dKKKT);
        results_average[58] += t.get(dKKvv);
        results_average[59] += t.get(dKKvT);
        results_average[60] += t.get(dKKTT);
        results_average[61] += t.get(dKvvv);
        results_average[62] += t.get(dKvvT);
        results_average[63] += t.get(dKvTT);
        results_average[64] += t.get(dKTTT);
        results_average[65] += t.get(dvvvv);
        results_average[66] += t.get(dvvvT);
        results_average[67] += t.get(dvvTT);
        results_average[68] += t.get(dvTTT);
        results_average[69] += t.get(dTTTT);
        results_average[70] += t.get(dSSSSS);
        results_average[71] += t.get(dSSSSK);
        results_average[72] += t.get(dSSSSv);
        results_average[73] += t.get(dSSSST);
        results_average[74] += t.get(dSSSKK);
        results_average[75] += t.get(dSSSKv);
        results_average[76] += t.get(dSSSKT);
        results_average[77] += t.get(dSSSvv);
        results_average[78] += t.get(dSSSvT);
        results_average[79] += t.get(dSSSTT);
        results_average[80] += t.get(dSSKKK);
        results_average[81] += t.get(dSSKKv);
        results_average[82] += t.get(dSSKKT);
        results_average[83] += t.get(dSSKvv);
        results_average[84] += t.get(dSSKvT);
        results_average[85] += t.get(dSSKTT);
        results_average[86] += t.get(dSSvvv);
        results_average[87] += t.get(dSSvvT);
        results_average[88] += t.get(dSSvTT);
        results_average[89] += t.get(dSSTTT);
        results_average[90] += t.get(dSKKKK);
        results_average[91] += t.get(dSKKKv);
        results_average[92] += t.get(dSKKKT);
        results_average[93] += t.get(dSKKvv);
        results_average[94] += t.get(dSKKvT);
        results_average[95] += t.get(dSKKTT);
        results_average[96] += t.get(dSKvvv);
        results_average[97] += t.get(dSKvvT);
        results_average[98] += t.get(dSKvTT);
        results_average[99] += t.get(dSKTTT);
        results_average[100] += t.get(dSvvvv);
        results_average[101] += t.get(dSvvvT);
        results_average[102] += t.get(dSvvTT);
        results_average[103] += t.get(dSvTTT);
        results_average[104] += t.get(dSTTTT);
        results_average[105] += t.get(dKKKKK);
        results_average[106] += t.get(dKKKKv);
        results_average[107] += t.get(dKKKKT);
        results_average[108] += t.get(dKKKvv);
        results_average[109] += t.get(dKKKvT);
        results_average[110] += t.get(dKKKTT);
        results_average[111] += t.get(dKKvvv);
        results_average[112] += t.get(dKKvvT);
        results_average[113] += t.get(dKKvTT);
        results_average[114] += t.get(dKKTTT);
        results_average[115] += t.get(dKvvvv);
        results_average[116] += t.get(dKvvvT);
        results_average[117] += t.get(dKvvTT);
        results_average[118] += t.get(dKvTTT);
        results_average[119] += t.get(dKTTTT);
        results_average[120] += t.get(dvvvvv);
        results_average[121] += t.get(dvvvvT);
        results_average[122] += t.get(dvvvTT);
        results_average[123] += t.get(dvvTTT);
        results_average[124] += t.get(dvTTTT);
        results_average[125] += t.get(dTTTTT);
    }

    time2 = std::chrono::high_resolution_clock::now();
    auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1)
            .count();

    std::cout << "iterations: " << iters << std::endl;
    std::cout << "adhoc order 5 time (ms): " << time << std::endl;

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    for (std::size_t i = 0; i < results_average.size(); ++i) {
        results_average[i] /= static_cast<double>(iters);
        std::cout << results_average[i] << std::endl;
    }

    return 0;
}
