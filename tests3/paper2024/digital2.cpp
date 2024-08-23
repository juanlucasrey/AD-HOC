#include "digital2.hpp"

#include <calc_tree.hpp>
#include <differential_operator/differential_operator.hpp>
#include <differential_operator/select_root_derivatives.hpp>
#include <init.hpp>
#include <tape2.hpp>

#include "../call_price.hpp"

namespace adhoc3 {

void digital_order2(double S_val, double K_val, double v_val, double T_val,
                    std::array<double, 15> &results) {

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
    auto dS = d<2>(S);
    auto dK = d<2>(K);
    auto dv = d<2>(v);
    auto dT = d<2>(T);

    // order 2
    auto dSS = d(S) * d(S);
    auto dKK = d(K) * d(K);
    auto dvv = d(v) * d(v);
    auto dTT = d(T) * d(T);
    auto dSK = d(S) * d(K);
    auto dSv = d(S) * d(v);
    auto dST = d(S) * d(T);
    auto dKv = d(K) * d(v);
    auto dKT = d(K) * d(T);
    auto dvT = d(v) * d(T);

    auto dr2 = d<2>(res);

    auto t = Tape2(dS, dK, dv, dT, dSS, dKK, dvv, dTT, dSK, dSv, dST, dKv, dKT,
                   dvT, dr2);

    t.set(dr2) = 1.;
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
}

} // namespace adhoc3
