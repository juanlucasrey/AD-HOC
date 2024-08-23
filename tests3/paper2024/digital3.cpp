#include "digital3.hpp"

#include <calc_tree.hpp>
#include <differential_operator/differential_operator.hpp>
#include <differential_operator/select_root_derivatives.hpp>
#include <init.hpp>
#include <tape2.hpp>

#include "../call_price.hpp"

namespace adhoc3 {

void digital_order3(double S_val, double K_val, double v_val, double T_val,
                    std::array<double, 35> &results) {

    auto [S, K, v, T, r] = Init<5>();

    auto res = digital_price(S, K, v, T);

    CalcTree ct(res);

    ct.set(S) = S_val;
    ct.set(K) = K_val;
    ct.set(v) = v_val;
    ct.set(T) = T_val;
    ct.evaluate();

    results[0] = ct.val(res);

    // order 1
    auto dS = d<3>(S);
    auto dK = d<3>(K);
    auto dv = d<3>(v);
    auto dT = d<3>(T);

    // order 2
    auto dSS = d<2>(S) * d(S);
    auto dKK = d<2>(K) * d(K);
    auto dvv = d<2>(v) * d(v);
    auto dTT = d<2>(T) * d(T);
    auto dSK = d<2>(S) * d(K);
    auto dSv = d<2>(S) * d(v);
    auto dST = d<2>(S) * d(T);
    auto dKv = d<2>(K) * d(v);
    auto dKT = d<2>(K) * d(T);
    auto dvT = d<2>(v) * d(T);

    // order 3
    auto dSSS = pow<3>(d(S));
    auto dKKK = pow<3>(d(K));
    auto dvvv = pow<3>(d(v));
    auto dTTT = pow<3>(d(T));

    auto dSSK = pow<2>(d(S)) * d(K);
    auto dSSv = pow<2>(d(S)) * d(v);
    auto dSST = pow<2>(d(S)) * d(T);
    auto dSKK = d(S) * pow<2>(d(K));
    auto dSvv = d(S) * pow<2>(d(v));
    auto dSTT = d(S) * pow<2>(d(T));
    auto dKKv = pow<2>(d(K)) * d(v);
    auto dKKT = pow<2>(d(K)) * d(T);
    auto dKvv = d(K) * pow<2>(d(v));
    auto dKTT = d(K) * pow<2>(d(T));
    auto dvvT = pow<2>(d(v)) * d(T);
    auto dvTT = d(v) * pow<2>(d(T));

    auto dSKv = d(S) * d(K) * d(v);
    auto dSKT = d(S) * d(K) * d(T);
    auto dSvT = d(S) * d(v) * d(T);
    auto dKvT = d(K) * d(v) * d(T);

    auto dr3 = d<3>(res);

    auto t =
        Tape2(dS, dK, dv, dT, dSS, dKK, dvv, dTT, dSK, dSv, dST, dKv, dKT, dvT,
              dSSS, dKKK, dvvv, dTTT, dSSK, dSSv, dSST, dSKK, dSvv, dSTT, dKKv,
              dKKT, dKvv, dKTT, dvvT, dvTT, dSKv, dSKT, dSvT, dKvT, dr3);

    t.set(dr3) = 1.;
    t.backpropagate(ct);

    results[1] = t.get_d(dS);
    results[2] = t.get_d(dK);
    results[3] = t.get_d(dv);
    results[4] = t.get_d(dT);
    results[5] = t.get_d(dSS);
    results[6] = t.get_d(dKK);
    results[7] = t.get_d(dvv);
    results[8] = t.get_d(dTT);
    results[9] = t.get_d(dSK);
    results[10] = t.get_d(dSv);
    results[11] = t.get_d(dST);
    results[12] = t.get_d(dKv);
    results[13] = t.get_d(dKT);
    results[14] = t.get_d(dvT);
    results[15] = t.get_d(dSSS);
    results[16] = t.get_d(dKKK);
    results[17] = t.get_d(dvvv);
    results[18] = t.get_d(dTTT);
    results[19] = t.get_d(dSSK);
    results[20] = t.get_d(dSSv);
    results[21] = t.get_d(dSST);
    results[22] = t.get_d(dSKK);
    results[23] = t.get_d(dSvv);
    results[24] = t.get_d(dSTT);
    results[25] = t.get_d(dKKv);
    results[26] = t.get_d(dKKT);
    results[27] = t.get_d(dKvv);
    results[28] = t.get_d(dKTT);
    results[29] = t.get_d(dvvT);
    results[30] = t.get_d(dvTT);
    results[31] = t.get_d(dSKv);
    results[32] = t.get_d(dSKT);
    results[33] = t.get_d(dSvT);
    results[34] = t.get_d(dKvT);
}

} // namespace adhoc3
