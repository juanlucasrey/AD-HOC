#include "digital1.hpp"

#include <calc_tree.hpp>
#include <differential_operator/differential_operator.hpp>
#include <differential_operator/select_root_derivatives.hpp>
#include <init.hpp>
#include <tape2.hpp>

#include "../call_price.hpp"

namespace adhoc3 {

void digital_order1(double S_val, double K_val, double v_val, double T_val,
                    std::array<double, 5> &results) {

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
    auto dS = d(S);
    auto dK = d(K);
    auto dv = d(v);
    auto dT = d(T);

    auto dr2 = d(res);

    auto t = Tape2(dS, dK, dv, dT, dr2);

    t.set(dr2) = 1.;
    t.backpropagate(ct);

    results[1] = t.get_d(dS);
    results[2] = t.get_d(dK);
    results[3] = t.get_d(dv);
    results[4] = t.get_d(dT);
}

} // namespace adhoc3
