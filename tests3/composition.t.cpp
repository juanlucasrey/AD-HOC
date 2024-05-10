#include <adhoc.hpp>
#include <differential_operator/differential_operator.hpp>
#include <init.hpp>
#include <tape2.hpp>

#include "call_price.hpp"

#include <gtest/gtest.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <random>

namespace adhoc3 {

// TEST(AD, BlackScholes_compositionOrder1) {
//     {
//         auto [x1, x2] = Init<2>();
//         auto y = log(log(x1 * x2) * x2);
//         CalcTree t(y);
//         t.set(x1) = 0.5;
//         t.set(x2) = 0.33;
//         t.evaluate();

//         Tape2 bp(d(x1), d(x2), d(y));
//         bp.set(d(y)) = 1.;
//         bp.backpropagate(t);

//         std::cout << bp.get(d(x1)) << std::endl;
//         std::cout << bp.get(d(x2)) << std::endl;
//     }
//     {
//         auto [x1, x2] = Init<2>();
//         auto y1 = log(x1 * x2);

//         CalcTree t1(y1);
//         t1.set(x1) = 0.5;
//         t1.set(x2) = 0.33;
//         t1.evaluate();

//         CalcTree t2(y1);
//         t2.set(x1) = t1.val(y1);
//         t2.set(x2) = 0.33;
//         t2.evaluate();

//         Tape2 bp(d(x1), d(x2), d(y1));
//         bp.set(d(y1)) = 1.;
//         bp.backpropagate(t2);

//         Tape2 bp2(d(x1), d(x2), d(y1));
//         bp2.set(d(y1)) = bp.get(d(x1));
//         bp2.set(d(x2)) = bp.get(d(x2));
//         bp2.backpropagate(t1);

//         std::cout << bp2.get(d(x1)) << std::endl;
//         std::cout << bp2.get(d(x2)) << std::endl;
//     }
// }

// TEST(AD, BlackScholes_compositionOrder2) {
//     double dx1_2 = 0;
//     double dx2_2 = 0;
//     double dx1_x2 = 0;
//     {
//         auto [x1, x2] = Init<2>();
//         auto y = erfc(erfc(x1 * x2) * x2);
//         CalcTree t(y);
//         t.set(x1) = 0.5;
//         t.set(x2) = 0.33;
//         t.evaluate();

//         Tape2 bp(pow<2>(d(x1)), pow<2>(d(x2)), d(x1) * d(x2), d<2>(y));
//         bp.set(d<2>(y)) = 1.;
//         bp.backpropagate(t);

//         std::cout << bp.get(pow<2>(d(x1))) << std::endl;
//         std::cout << bp.get(pow<2>(d(x2))) << std::endl;
//         std::cout << bp.get(d(x1) * d(x2)) << std::endl;

//         dx1_2 = bp.get(pow<2>(d(x1)));
//         dx2_2 = bp.get(pow<2>(d(x2)));
//         dx1_x2 = bp.get(d(x1) * d(x2));
//     }
//     {
//         auto [x1, x2] = Init<2>();
//         auto y1 = erfc(x1 * x2);

//         CalcTree t1(y1);
//         t1.set(x1) = 0.5;
//         t1.set(x2) = 0.33;
//         t1.evaluate();

//         CalcTree t2(y1);
//         t2.set(x1) = t1.val(y1);
//         t2.set(x2) = 0.33;
//         t2.evaluate();

//         Tape2 bp(d<2>(y1), pow<2>(d(x1)), pow<2>(d(x2)), d(x1) * d(x2),
//                  d<2>(x1) /* , d<2>(x2) */);
//         bp.set(d<2>(y1)) = 1.;
//         bp.backpropagate(t2);

//         Tape2 bp2(pow<2>(d(y1)), d(y1) * d(x2), d<2>(y1),
//                   d<2>(x1) /* , d<2>(x2) */, pow<2>(d(x1)), pow<2>(d(x2)),
//                   d(x1) * d(x2));

//         bp2.set(pow<2>(d(y1))) = bp.get(pow<2>(d(x1)));
//         bp2.set(pow<2>(d(x2))) = bp.get(pow<2>(d(x2)));
//         bp2.set(d(y1) * d(x2)) = bp.get(d(x1) * d(x2));
//         bp2.set(d<2>(y1)) = bp.get(d<2>(x1));
//         // bp2.set(d<2>(x2)) = bp.get(d<2>(x2));

//         bp2.backpropagate(t1);
//         std::cout << bp2.get(pow<2>(d(x1))) << std::endl;
//         std::cout << bp2.get(pow<2>(d(x2))) << std::endl;
//         std::cout << bp2.get(d(x1) * d(x2)) << std::endl;

//         EXPECT_NEAR(dx1_2, bp2.get(pow<2>(d(x1))), 1e-15);
//         EXPECT_NEAR(dx2_2, bp2.get(pow<2>(d(x2))), 1e-15);
//         EXPECT_NEAR(dx1_x2, bp2.get(d(x1) * d(x2)), 1e-15);
//     }
// }

// TEST(AD, BlackScholes_compositionOrder1And2) {
//     std::cout.precision(std::numeric_limits<double>::max_digits10);

//     double dx1_2 = 0;
//     double dx2_2 = 0;
//     double dx1_x2 = 0;
//     {
//         auto [x1, x2] = Init<2>();
//         auto y = erfc(erfc(x1 * x2) * x2);
//         CalcTree t(y);
//         t.set(x1) = 0.5;
//         t.set(x2) = 0.33;
//         t.evaluate();

//         auto dx1 = d(x1);
//         auto dx2 = d(x2);
//         auto dy = d(y);

//         auto d2x1 = pow<2>(d(x1));
//         auto d2x2 = pow<2>(d(x2));
//         auto d2x1x2 = d(x1) * d(x2);
//         auto d2y = d<2>(y);

//         Tape2 bp(dx1, dx2, dy, d2x1, d2x2, d2x1x2, d<2>(x1), d<2>(x2), d2y);
//         bp.set(dy) = 1.;
//         bp.set(d2y) = 1.;
//         bp.backpropagate(t);

//         std::cout << "final" << std::endl;
//         std::cout << bp.get(d(x1)) << std::endl;
//         std::cout << bp.get(d(x2)) << std::endl;

//         std::cout << bp.get(d<2>(x1)) << std::endl;
//         std::cout << bp.get(d<2>(x2)) << std::endl;
//         std::cout << bp.get(pow<2>(d(x1))) << std::endl;
//         std::cout << bp.get(pow<2>(d(x2))) << std::endl;
//         std::cout << bp.get(d(x1) * d(x2)) << std::endl;

//         dx1_2 = bp.get(d2x1);
//         dx2_2 = bp.get(d2x2);
//         dx1_x2 = bp.get(d2x1x2);
//     }
//     {
//         auto [x1, x2] = Init<2>();
//         auto y1 = erfc(x1 * x2);

//         CalcTree t1(y1);
//         t1.set(x1) = 0.5;
//         t1.set(x2) = 0.33;
//         t1.evaluate();

//         CalcTree t2(y1);
//         t2.set(x1) = t1.val(y1);
//         t2.set(x2) = 0.33;
//         t2.evaluate();

//         Tape2 bp(d(x1), d(x2), d(y1), pow<2>(d(x1)), pow<2>(d(x2)),
//                  d(x1) * d(x2), d<2>(x1), d<2>(x2), d<2>(y1));
//         bp.set(d(y1)) = 1.;
//         bp.set(d<2>(y1)) = 1.;
//         bp.backpropagate(t2);

//         std::cout << bp.get(d(x1)) << std::endl;
//         std::cout << bp.get(d(x2)) << std::endl;

//         std::cout << bp.get(d<2>(x1)) << std::endl;
//         std::cout << bp.get(d<2>(x2)) << std::endl;
//         std::cout << bp.get(pow<2>(d(x1))) << std::endl;
//         std::cout << bp.get(pow<2>(d(x2))) << std::endl;
//         std::cout << bp.get(d(x1) * d(x2)) << std::endl;

//         auto dyx2 = d(y1) * d(x2);
//         Tape2 bp2(d(x1), d(x2), d(y1), pow<2>(d(x1)), pow<2>(d(x2)),
//                   d(x1) * d(x2), d<2>(y1), dyx2, pow<2>(d(y1)), d<2>(x1),
//                   d<2>(x2));
//         bp2.set(d(y1)) = bp.get(d(x1));
//         bp2.set(d(x2)) = bp.get(d(x2));

//         bp2.set(d<2>(y1)) = bp.get(d<2>(x1));
//         bp2.set(d<2>(x2)) = bp.get(d<2>(x2));
//         bp2.set(pow<2>(d(y1))) = bp.get(pow<2>(d(x1)));
//         bp2.set(pow<2>(d(x2))) = bp.get(pow<2>(d(x2)));
//         bp2.set(d(y1) * d(x2)) = bp.get(d(x1) * d(x2));
//         bp2.set(d<2>(x1)) = bp.get(d<2>(y1));

//         bp2.backpropagate(t1);

//         std::cout << "final" << std::endl;
//         std::cout << bp2.get(d(x1)) << std::endl;
//         std::cout << bp2.get(d(x2)) << std::endl;

//         std::cout << bp2.get(d<2>(x1)) << std::endl; // wrong??
//         std::cout << bp2.get(d<2>(x2)) << std::endl;
//         std::cout << bp2.get(pow<2>(d(x1))) << std::endl;
//         std::cout << bp2.get(pow<2>(d(x2))) << std::endl;
//         std::cout << bp2.get(d(x1) * d(x2)) << std::endl;

//         EXPECT_NEAR(dx1_2, bp2.get(pow<2>(d(x1))), 1e-15);
//         EXPECT_NEAR(dx2_2, bp2.get(pow<2>(d(x2))), 1e-15);
//         EXPECT_NEAR(dx1_x2, bp2.get(d(x1) * d(x2)), 1e-15);
//     }
// }

// // template <class... Ids, std::size_t... Orders, std::size_t... Powers>
// // constexpr auto
// // multiplicity2(std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...> /* in
// */)
// //     -> std::size_t {
// //     return MultinomialCoeff2(std::index_sequence<Powers * Orders...>{});
// // }

// TEST(AD, BlackScholes_compositionOrder3) {
//     std::cout.precision(std::numeric_limits<double>::max_digits10);

//     auto [x1, x2] = Init<2>();
//     auto y = erfc(erfc(x1 * x2) * x2);
//     CalcTree t(y);
//     t.set(x1) = 0.5;
//     t.set(x2) = 0.33;
//     t.evaluate();

//     // auto dx1 = d(x1);
//     // auto dx2 = d(x1);
//     // auto dy = d(y);

//     // auto d2x1 = pow<2>(d(x1));
//     // auto d2x2 = pow<2>(d(x2));
//     // auto d2x1x2 = d(x1) * d(x2);
//     // auto d2y = d<2>(y);

//     Tape2 bp(d(x1), d(x2), d(y), pow<2>(d(x1)), pow<2>(d(x2)), d(x1) * d(x2),
//              d<2>(x1), d<2>(x2), d<2>(y), pow<3>(d(x1)), pow<3>(d(x2)),
//              pow<2>(d(x1)) * d(x2), d(x1) * pow<2>(d(x2)), d<3>(x1),
//              d<3>(x2), d<2>(x1) * d(x1), d<2>(x1) * d(x2), d<2>(x2) * d(x1),
//              d<2>(x2) * d(x2), d<3>(y), pow<4>(d(x1)), pow<3>(d(x1)) * d(x2),
//              pow<2>(d(x1)) * pow<2>(d(x2)), d(x1) * pow<3>(d(x2)),
//              pow<4>(d(x2)), d<4>(x1), d<3>(x1) * d(x1), d<3>(x1) * d(x2),
//              pow<2>(d<2>(x1)), d<2>(x1) * pow<2>(d(x1)),
//              d<2>(x1) * d(x1) * d(x2), d<2>(x1) * pow<2>(d(x2)), d<4>(y));
//     bp.set(d(y)) = 1.;
//     bp.set(d<2>(y)) = 1.;
//     bp.set(d<3>(y)) = 1.;
//     bp.set(d<4>(y)) = 1.;
//     bp.backpropagate(t);

//     std::cout << "1" << std::endl;
//     std::cout << bp.get_d(d(x1)) << std::endl;
//     std::cout << bp.get_d(d(x2)) << std::endl;

//     std::cout << "2" << std::endl;
//     std::cout << bp.get_d(pow<2>(d(x1))) << std::endl;
//     std::cout << bp.get_d(pow<2>(d(x2))) << std::endl;
//     std::cout << bp.get_d(d(x1) * d(x2)) << std::endl;

//     std::cout << "2r" << std::endl;
//     std::cout << bp.get_d(d<2>(x1)) << std::endl;
//     std::cout << bp.get_d(d<2>(x2)) << std::endl;

//     std::cout << "3" << std::endl;
//     std::cout << bp.get_d(pow<3>(d(x1))) << std::endl;
//     std::cout << bp.get_d(pow<3>(d(x2))) << std::endl;
//     std::cout << bp.get_d(pow<2>(d(x1)) * d(x2)) << std::endl;
//     std::cout << bp.get_d(d(x1) * pow<2>(d(x2))) << std::endl;

//     std::cout << "3r" << std::endl;
//     std::cout << bp.get_d(d<3>(x1)) << std::endl;
//     std::cout << bp.get_d(d<3>(x2)) << std::endl;
//     std::cout << bp.get_d(d<2>(x1) * d(x1)) << std::endl;
//     std::cout << bp.get_d(d<2>(x1) * d(x2)) << std::endl;
//     std::cout << bp.get_d(d<2>(x2) * d(x1)) << std::endl;
//     std::cout << bp.get_d(d<2>(x2) * d(x2)) << std::endl;

//     std::cout << "4" << std::endl;
//     std::cout << bp.get_d(pow<4>(d(x1))) << std::endl;
//     std::cout << bp.get_d(pow<3>(d(x1)) * d(x2)) << std::endl;
//     std::cout << bp.get_d(pow<2>(d(x1)) * pow<2>(d(x2))) << std::endl;
//     std::cout << bp.get_d(d(x1) * pow<3>(d(x2))) << std::endl;
//     std::cout << bp.get_d(pow<4>(d(x2))) << std::endl;

//     std::cout << "4r" << std::endl;
//     std::cout << bp.get_d(d<4>(x1)) << std::endl;
//     std::cout << bp.get_d(d<3>(x1) * d(x1)) << std::endl;
//     std::cout << bp.get_d(d<3>(x1) * d(x2)) << std::endl;
//     std::cout << bp.get_d(pow<2>(d<2>(x1))) << std::endl;
//     std::cout << bp.get(pow<2>(d<2>(x1))) << std::endl;
//     std::cout << bp.get_d(d<2>(x1) * pow<2>(d(x1))) << std::endl;
//     std::cout << bp.get_d(d<2>(x1) * d(x1) * d(x2)) << std::endl;
//     std::cout << bp.get_d(d<2>(x1) * pow<2>(d(x2))) << std::endl;
// }

TEST(AD, CrossDeriv) {
    std::cout.precision(std::numeric_limits<double>::max_digits10);

    auto [x1, x2] = Init<2>();
    auto y = erfc(x1 * x2);
    CalcTree t(y);
    t.set(x1) = 0.5;
    t.set(x2) = 0.33;
    t.evaluate();

    // auto dx1 = d(x1);
    // auto dx2 = d(x1);
    // auto dy = d(y);

    // auto d2x1 = pow<2>(d(x1));
    // auto d2x2 = pow<2>(d(x2));
    // auto d2x1x2 = d(x1) * d(x2);
    // auto d2y = d<2>(y);

    Tape2 bp(d(x1), d(x2), d(y), pow<2>(d(x1)), pow<2>(d(x2)), d(x1) * d(x2),
             d<2>(x1), d<2>(x2), d<2>(y), pow<3>(d(x1)), pow<3>(d(x2)),
             pow<2>(d(x1)) * d(x2), d(x1) * pow<2>(d(x2)), d<3>(x1), d<3>(x2),
             d<2>(x1) * d(x1), d<2>(x1) * d(x2), d<2>(x2) * d(x1),
             d<2>(x2) * d(x2), d<3>(y), pow<4>(d(x1)), pow<3>(d(x1)) * d(x2),
             pow<2>(d(x1)) * pow<2>(d(x2)), d(x1) * pow<3>(d(x2)),
             pow<4>(d(x2)), d<4>(x1), d<3>(x1) * d(x1), d<3>(x1) * d(x2),
             pow<2>(d<2>(x1)), d<2>(x1) * pow<2>(d(x1)),
             d<2>(x1) * d(x1) * d(x2), d<2>(x1) * pow<2>(d(x2)),
             d<2>(x1) * d<2>(x2), d<4>(y));
    bp.set(d(y)) = 1.;
    bp.set(d<2>(y)) = 1.;
    bp.set(d<3>(y)) = 1.;
    bp.set(d<4>(y)) = 1.;
    bp.backpropagate(t);

    std::cout << "1" << std::endl;
    std::cout << bp.get(d(x1)) << std::endl;
    std::cout << bp.get(d(x2)) << std::endl;

    std::cout << "2" << std::endl;
    std::cout << bp.get(pow<2>(d(x1))) << std::endl;
    std::cout << bp.get(pow<2>(d(x2))) << std::endl;
    std::cout << bp.get(d(x1) * d(x2)) << std::endl;

    std::cout << "2r" << std::endl;
    std::cout << bp.get(d<2>(x1)) << std::endl;
    std::cout << bp.get(d<2>(x2)) << std::endl;

    std::cout << "3" << std::endl;
    std::cout << bp.get(pow<3>(d(x1))) << std::endl;
    std::cout << bp.get(pow<3>(d(x2))) << std::endl;
    std::cout << bp.get(pow<2>(d(x1)) * d(x2)) << std::endl;
    std::cout << bp.get(d(x1) * pow<2>(d(x2))) << std::endl;

    std::cout << "3r" << std::endl;
    std::cout << bp.get(d<3>(x1)) << std::endl;
    std::cout << bp.get(d<3>(x2)) << std::endl;
    std::cout << bp.get(d<2>(x1) * d(x1)) << std::endl;
    std::cout << bp.get(d<2>(x1) * d(x2)) << std::endl;
    std::cout << bp.get(d<2>(x2) * d(x1)) << std::endl;
    std::cout << bp.get(d<2>(x2) * d(x2)) << std::endl;

    std::cout << "4" << std::endl;
    std::cout << bp.get(pow<4>(d(x1))) << std::endl;
    std::cout << bp.get(pow<3>(d(x1)) * d(x2)) << std::endl;
    std::cout << bp.get(pow<2>(d(x1)) * pow<2>(d(x2))) << std::endl;
    std::cout << bp.get(d(x1) * pow<3>(d(x2))) << std::endl;
    std::cout << bp.get(pow<4>(d(x2))) << std::endl;

    std::cout << "4r" << std::endl;
    std::cout << bp.get(d<4>(x1)) << std::endl;
    std::cout << bp.get(d<3>(x1) * d(x1)) << std::endl;
    std::cout << bp.get(d<3>(x1) * d(x2)) << std::endl;
    std::cout << bp.get(pow<2>(d<2>(x1))) << std::endl;
    std::cout << bp.get(d<2>(x1) * pow<2>(d(x1))) << std::endl;
    std::cout << bp.get(d<2>(x1) * d(x1) * d(x2)) << std::endl;
    std::cout << bp.get(d<2>(x1) * pow<2>(d(x2))) << std::endl;

    std::cout << bp.get(d<2>(x1) * d<2>(x2)) << std::endl;
}

} // namespace adhoc3
