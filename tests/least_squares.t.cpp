#define DCO_EXT_NO_EIGEN
#define DCO_AUTO_SUPPORT
#define DCO_DISABLE_AVX2_WARNING
#define DCO_CHUNK_TAPE
#include "../benchmark/packages/dco/dco.hpp"

#include <Eigen/Dense>
#include <gtest/gtest.h>

#include <limits>
#include <random>

namespace BloombergLP::qfd::test {

namespace {

template <class D> using EigenArray = typename Eigen::Array<D, -1, 1>;
template <class D> using EigenArrayXX = typename Eigen::Array<D, -1, -1>;
template <class D> using EigenMatrix = typename Eigen::Matrix<D, -1, -1>;
template <class D> using EigenVector = typename Eigen::Matrix<D, -1, 1>;

constexpr std::size_t numRegressors = 3;
// constexpr double pvTarget = 0.234135;
constexpr double pvTarget = 0.18556876913162901;
// constexpr double dBeta1dSpotTarget = 0.060057;
constexpr double dBeta1dSpotTarget = 0.054958418524778369;
// constexpr double dBeta1dVolTarget = 0.029752;
constexpr double dBeta1dVolTarget = 0.011737403136111536;
// constexpr double dPVdSpotTarget = 0.0102026;
constexpr double dPVdSpotTarget = 0.0086732402860533551;
// constexpr double dPVdVolTarget = 0.00454295;
constexpr double dPVdVolTarget = 0.0027125387140718568;

template <class D> inline auto cdf(D const &x) {
    const double one_over_root_two = 0.70710678118654757;
    using std::erfc;
    return 0.5 * erfc(x * -one_over_root_two);
};

template <class D>
static void Calculate_d1_d2(D const &forward, D const &strike,
                            D const &total_vol, D &d1, D &d2) {
    using std::log;
    D log_term = log(forward / strike) / total_vol;
    D half_vol = total_vol * 0.5;
    d1 = log_term + half_vol;
    d2 = log_term - half_vol;
}

template <class D>
D CallBSPV(D const &fDF, D const &fK, D const &fFwd, D const &fTotalVol) {
    if (fK <= 0.0) {
        // No more an option, but a model-independent forward price
        return fDF * (fFwd - fK);
    }

    D d1;
    D d2;
    Calculate_d1_d2(fFwd, fK, fTotalVol, d1, d2);
    return fDF * (fFwd * cdf(d1) - fK * cdf(d2));
}

template <typename D>
D futureSpot(const D &todaysSpot, const D &vol, std::mt19937 &generator,
             std::normal_distribution<double> &distribution) {
    const double W = distribution(generator);
    return todaysSpot * exp(vol * W - 0.5 * vol * vol);
}

template <typename D> D price(const D &fwd, const D &vol) {
    const double time = 1.25;
    const D totalVol = vol * std::sqrt(time);
    const double df = 0.98;
    const double strike = 1.0;
    return CallBSPV(D(df), D(strike), fwd, totalVol);
}

template <typename D>
EigenVector<D> solveForBeta(const EigenMatrix<D> &XTX,
                            const EigenVector<D> &XTy) {
    return (XTX).colPivHouseholderQr().solve(XTy);
}

template <typename B, typename EAO> void adjointBeta(EAO *eao) {
    /********************************************************/
    /* Implicitly differentiate w.r.t regression input r,   */
    /*   then multiply by betaAdj:                          */
    /*    XTX * beta = XTy                                  */
    /*    (dXTXdr * beta + XTX *dBetadr - dXTydr)*betaAdj   */
    /*    Since rAdj = dBetadr * betaAdj, we get an implict */
    /*    solution for rAdj:                                */
    /*    rAdj = XTXInv*(dXTydr - dXTXdr * beta)*betaAdj    */
    /*    and similarly for yAdj.                           */
    /********************************************************/
    const std::size_t numPresimPaths = eao->template read_data<std::size_t>();
    const EigenMatrix<B> XPrimal = eao->template read_data<EigenMatrix<B>>();
    const EigenVector<B> YPrimal = eao->template read_data<EigenVector<B>>();
    const EigenMatrix<B> XTXPrimal = eao->template read_data<EigenMatrix<B>>();
    const EigenVector<B> betaPrimal = eao->template read_data<EigenVector<B>>();
    const EigenMatrix<B> XTXInv = XTXPrimal.inverse();
    EigenVector<B> betaAdj(numRegressors);

    EigenMatrix<double> XAdj =
        EigenMatrix<double>::Zero(numPresimPaths, numRegressors);
    EigenVector<double> YAdj = EigenVector<double>::Zero(numPresimPaths);

    for (std::size_t i = 0; i < numRegressors; ++i) {
        betaAdj(i) = eao->get_output_adjoint();
    }

    for (std::size_t path = 0; path < numPresimPaths; ++path) {
        EigenVector<B> dXTydy = EigenVector<B>::Zero(numRegressors);

        for (std::size_t r = 0; r < numRegressors; ++r) {
            EigenVector<B> dXTydr = EigenVector<B>::Zero(numRegressors);
            EigenMatrix<B> dXTXdr =
                EigenMatrix<B>::Zero(numRegressors, numRegressors);

            for (std::size_t i = 0; i < numRegressors; ++i) {
                for (std::size_t j = 0; j < numRegressors; ++j) {
                    if ((r == i) && (r == j)) {
                        dXTXdr(i, j) = 2 * XPrimal(path, i);
                    } else if (r == i) {
                        dXTXdr(i, j) = XPrimal(path, j);
                    } else if (r == j) {
                        dXTXdr(i, j) = XPrimal(path, i);
                    }
                }
            }

            dXTydr(r) = YPrimal(path);
            dXTydy(r) = XPrimal(path, r);
            EigenVector<B> symbolicAdjoint =
                XTXInv * (dXTydr - dXTXdr * betaPrimal);
            XAdj(path, r) = symbolicAdjoint.dot(betaAdj);
            eao->increment_input_adjoint(XAdj(path, r));
        }

        YAdj(path) = (XTXInv * dXTydy).dot(betaAdj);
        eao->increment_input_adjoint(YAdj(path));
    }
}

template <typename D>
std::tuple<D, D, D, D>
preSimPathXy(const D &spot, const D &vol, std::mt19937 &generator,
             std::normal_distribution<double> &distribution) {
    const D fwd = futureSpot(spot, vol, generator, distribution);
    const D y = price(fwd, vol);
    return {1.0, fwd, fwd * fwd, y};
}

template <typename B, typename EAO> void adjointPerPathPreSimulation(EAO *eao) {
    using dco_mode = dco::ga1v<B, 1>;
    using D = typename dco_mode::type;
    using tape_type = typename dco_mode::tape_t;

    const std::size_t tapeSymbolicStart =
        dco::size_of(dco_mode::global_tape, tape_type::size_of_default);
    std::cout << "TapeSymbolicStart: " << tapeSymbolicStart << "\n";

    const D eurusdAdj = eao->template read_data<B>();
    const D volAdj = eao->template read_data<B>();
    const std::size_t numPresimPaths = eao->template read_data<std::size_t>();
    const EigenVector<B> y = eao->template read_data<EigenVector<B>>();
    const EigenMatrix<B> X = eao->template read_data<EigenMatrix<B>>();
    const EigenVector<B> beta = eao->template read_data<EigenVector<B>>();
    const EigenMatrix<B> XTXInv = eao->template read_data<EigenMatrix<B>>();
    auto tapePosition1 = dco_mode::global_tape->get_position();
    dco_mode::global_tape->register_variable(eurusdAdj);
    dco_mode::global_tape->register_variable(volAdj);

    EigenVector<B> betaAdj(numRegressors);

    for (std::size_t i = 0; i < numRegressors; ++i) {
        betaAdj(i) = eao->get_output_adjoint();
    }

    // compute adjoints via Mike Giles 2008: matrix derivatives
    const EigenVector<B> XTyAdjoint = XTXInv.transpose() * betaAdj;
    const EigenMatrix<B> XTXAdjoint = -XTyAdjoint * beta.transpose();
    const EigenVector<B> yAdj = X * XTyAdjoint;
    // contribututions to XAdj from XTXAdjoint (2) and XTyAdjoint (1)
    const EigenMatrix<B> XAdj = X * XTXAdjoint + X * XTXAdjoint.transpose() +
                                y * XTyAdjoint.transpose();

    // re-run pre simulation
    std::mt19937 generator(1u);
    std::normal_distribution<double> distribution(0., 1.);

    const std::size_t tapeSymbolicBeforePreSim =
        dco::size_of(dco_mode::global_tape, tape_type::size_of_default);
    EXPECT_EQ(tapeSymbolicBeforePreSim, tapeSymbolicStart + 32);
    std::size_t maxTape = 0.;

    for (std::size_t path = 0; path < numPresimPaths; ++path) {
        auto pos = dco_mode::global_tape->get_position();
        std::vector<D> Xrow(numRegressors);
        D y = 0.;
        std::tie(Xrow[0], Xrow[1], Xrow[2], y) =
            preSimPathXy(eurusdAdj, volAdj, generator, distribution);

        for (std::size_t j = 0; j < numRegressors; ++j) {
            dco_mode::global_tape->register_output_variable(Xrow[j]);
            // seed outputs
            dco::derivative(Xrow[j]) = XAdj(path, j);
        }

        dco_mode::global_tape->register_output_variable(y);
        dco::derivative(y) = yAdj(path);
        const std::size_t tapeOnePath =
            dco::size_of(dco_mode::global_tape, tape_type::size_of_default);
        maxTape = std::max(maxTape, tapeOnePath);
        dco_mode::global_tape->interpret_adjoint_and_reset_to(pos);
    }

    const std::size_t tapeSymbolicAfterPreSim =
        dco::size_of(dco_mode::global_tape, tape_type::size_of_default);
    EXPECT_LT(tapeSymbolicAfterPreSim, 2 * tapeSymbolicStart);
    EXPECT_LT(maxTape, 1200);

    eao->increment_input_adjoint(dco::derivative(eurusdAdj));
    eao->increment_input_adjoint(dco::derivative(volAdj));
    dco_mode::global_tape->reset_to(tapePosition1);
}

template <typename D>
EigenVector<D> preSimulation(const D &todaysSpot, const D &vol,
                             std::size_t numPresimPaths,
                             std::mt19937 &generator,
                             std::normal_distribution<double> &distribution) {
    EigenMatrix<D> XTX = EigenMatrix<D>::Zero(numRegressors, numRegressors);
    EigenVector<D> XTy = EigenVector<D>::Zero(numRegressors);
    const double weight = 1.0 / static_cast<double>(numPresimPaths);
    std::vector<D> regressors(numRegressors);
    D y = 0.;

    for (std::size_t path = 0; path < numPresimPaths; ++path) {
        std::tie(regressors[0], regressors[1], regressors[2], y) =
            preSimPathXy(todaysSpot, vol, generator, distribution);

        for (std::size_t i = 0; i < numRegressors; ++i) {
            for (std::size_t j = 0; j < numRegressors; ++j) {
                // online update of smaller matrix
                XTX(i, j) += regressors[i] * regressors[j] * weight;
            }
            XTy(i, 0) += regressors[i] * y * weight;
        }
    }

    return solveForBeta(XTX, XTy);
}

template <typename D>
EigenVector<D> preSimulationSymbolicAdjoint(
    const D &todaysSpot, const D &vol, std::size_t numPresimPaths,
    std::mt19937 &generator, std::normal_distribution<double> &distribution) {
    using dco_mode = typename dco::mode<D>;
    using EAO = typename dco_mode::external_adjoint_object_t;
    EAO *eao = dco_mode::global_tape->template create_callback_object<EAO>();
    EigenMatrix<double> X =
        EigenMatrix<double>::Zero(numPresimPaths, numRegressors);
    EigenVector<double> Y = EigenVector<double>::Zero(numPresimPaths);
    EigenMatrix<double> XTX =
        EigenMatrix<double>::Zero(numRegressors, numRegressors);
    EigenVector<double> XTy = EigenVector<double>::Zero(numRegressors);
    std::vector<D> regressors(numRegressors);
    D y = 0.;

    for (std::size_t path = 0; path < numPresimPaths; ++path) {
        std::tie(regressors[0], regressors[1], regressors[2], y) =
            preSimPathXy(todaysSpot, vol, generator, distribution);
        X(path, 0) = eao->register_input(regressors[0]);
        X(path, 1) = eao->register_input(regressors[1]);
        X(path, 2) = eao->register_input(regressors[2]);
        Y(path) = eao->register_input(y);

        for (std::size_t i = 0; i < numRegressors; ++i) {
            for (std::size_t j = 0; j < numRegressors; ++j) {
                // online update of smaller matrix
                XTX(i, j) += X(path, i) * X(path, j);
            }
            XTy(i, 0) += X(path, i) * Y(path);
        }
    }

    // primal solution
    auto betaPrimal = solveForBeta(XTX, XTy);
    EigenVector<D> beta(numRegressors);

    for (std::size_t i = 0; i < numRegressors; ++i) {
        beta(i) = eao->register_output(betaPrimal(i));
    }

    eao->write_data(numPresimPaths);
    eao->write_data(X);
    eao->write_data(Y);
    eao->write_data(XTX);
    eao->write_data(betaPrimal);
    dco_mode::global_tape->template insert_callback<EAO>(
        adjointBeta<double, EAO>, eao);
    return beta;
}

template <typename D>
EigenVector<D>
preSimulationFancy(const D &spot, const D &vol, std::size_t numPresimPaths,
                   std::mt19937 &generator,
                   std::normal_distribution<double> &distribution) {
    using dco_mode = typename dco::mode<D>;
    using EAO = typename dco_mode::external_adjoint_object_t;

    EAO *eao = dco_mode::global_tape->template create_callback_object<EAO>();
    eao->register_input(spot);
    eao->register_input(vol);
    const double spotPrimal = dco::value(spot);
    const double volPrimal = dco::value(vol);
    eao->write_data(spotPrimal);
    eao->write_data(volPrimal);
    eao->write_data(numPresimPaths);
    // do not record heavy computation on tape now
    dco_mode::global_tape
        ->switch_to_passive(); // probably overkill, as below is double-only
    EigenMatrix<double> X =
        EigenMatrix<double>::Zero(numPresimPaths, numRegressors);
    EigenVector<double> y = EigenVector<double>::Zero(numPresimPaths);
    EigenMatrix<double> XTX =
        EigenMatrix<double>::Zero(numRegressors, numRegressors);
    EigenVector<double> XTy = EigenVector<double>::Zero(numRegressors);

    for (std::size_t path = 0; path < numPresimPaths; ++path) {

        std::tie(X(path, 0), X(path, 1), X(path, 2), y(path)) =
            preSimPathXy(spotPrimal, volPrimal, generator, distribution);

        for (std::size_t i = 0; i < numRegressors; ++i) {
            for (std::size_t j = 0; j < numRegressors; ++j) {
                // online update of smaller matrix
                XTX(i, j) += X(path, i) * X(path, j);
            }
            XTy(i, 0) += X(path, i) * y(path);
        }
    }

    EigenVector<double> betaPrimal = solveForBeta(XTX, XTy);
    EigenMatrix<double> XTXInv = XTX.inverse();
    dco_mode::global_tape->switch_to_active(); // turn recording back on
    eao->write_data(y);
    eao->write_data(X);
    eao->write_data(betaPrimal);
    eao->write_data(XTXInv);
    // eao->write_data(generator);
    // eao->write_data(distribution);
    EigenVector<D> beta(numRegressors);

    for (std::size_t i = 0; i < numRegressors; ++i) {
        beta(i) = eao->register_output(betaPrimal(i));
    }

    dco_mode::global_tape->template insert_callback<EAO>(
        adjointPerPathPreSimulation<double, EAO>, eao);
    return beta;
}

template <typename D>
double
mainSimulationWithReplay(const D &spot, const D &vol,
                         const EigenVector<D> &beta, std::size_t numPaths,
                         std::mt19937 &generator,
                         std::normal_distribution<double> &distribution) {
    using dco_mode = typename dco::mode<D>;
    EigenVector<D> regressors(numRegressors);
    const double weight = 1.0 / static_cast<double>(numPaths);
    double pv = 0.;

    for (std::size_t path = 0; path < numPaths; ++path) {
        auto pos = dco_mode::global_tape->get_position();
        const D fwd = futureSpot(spot, vol, generator, distribution);
        regressors << 1.0, fwd, fwd * fwd;
        const D approx = beta.dot(regressors) * weight;
        dco_mode::global_tape->register_output_variable(approx);
        dco::derivative(approx)[0] = 1.0;
        pv += dco::passive_value(approx);
        dco_mode::global_tape->interpret_adjoint_and_reset_to(pos);
    }

    return pv;
}

template <typename B, typename EAO> void adjointCheckpoint(EAO *eao) {
    using dco_mode = dco::ga1v<B, 1>;
    using D = typename dco_mode::type;
    using tape_type = typename dco_mode::tape_t;

    const D eurusdAdj = eao->template read_data<B>();
    const D volAdj = eao->template read_data<B>();
    const std::size_t numPresimPaths = eao->template read_data<std::size_t>();

    const std::size_t tapeBefore =
        dco::size_of(dco_mode::global_tape, tape_type::size_of_default);

    auto tapePosition1 = dco_mode::global_tape->get_position();
    dco_mode::global_tape->register_variable(eurusdAdj);
    dco_mode::global_tape->register_variable(volAdj);
    auto tapePosition2 = dco_mode::global_tape->get_position();
    // re-run pre simulation
    std::mt19937 generator(1u);
    std::normal_distribution<double> distribution(0., 1.);

    auto betaAdj = preSimulation(eurusdAdj, volAdj, numPresimPaths, generator,
                                 distribution);

    for (std::size_t i = 0; i < numRegressors; ++i) {
        dco::derivative(betaAdj(i)) = eao->get_output_adjoint();
    }

    const std::size_t tapeAfter =
        dco::size_of(dco_mode::global_tape, tape_type::size_of_default);
    EXPECT_GT(tapeAfter, 380 * tapeBefore);
    dco_mode::global_tape->interpret_adjoint_and_reset_to(tapePosition2);
    eao->increment_input_adjoint(dco::derivative(eurusdAdj));
    eao->increment_input_adjoint(dco::derivative(volAdj));
    dco_mode::global_tape->reset_to(tapePosition1);
}
} // namespace

TEST(LeastSquares, RegressionErrors) {

    std::mt19937 generator(1u);
    std::normal_distribution<double> distribution(0., 1.);

    using D = double;
    const D eurusd = 1.2;
    const D vol = 0.15;
    const auto beta = preSimulation(eurusd, vol, 100, generator, distribution);

    // Main simulation
    const std::size_t numPaths = 10;
    EigenVector<D> regressors(numRegressors);
    double pv = 0.;

    for (std::size_t path = 0; path < numPaths; ++path) {
        const D fwd = futureSpot(eurusd, vol, generator, distribution);
        regressors << 1.0, fwd, fwd * fwd;
        const D pathApprox = beta.dot(regressors);
        const D pathPv = price(fwd, vol);
        const double error = dco::passive_value(pathApprox - pathPv);
        EXPECT_LT(std::abs(error), 2e-2);
        pv += pathApprox / numPaths;
    }

    EXPECT_NEAR(pvTarget, pv, 1e-6);
}

TEST(LeastSquares, TargetDerivativesFromTangent) {
    std::mt19937 generator(1u);
    std::normal_distribution<double> distribution(0., 1.);

    using D = dco::gt1v<double, 8>::type;
    D h0 = 0.;
    D h1 = 0.;
    dco::derivative(h0)[0] = 0.01;
    dco::derivative(h1)[1] = 0.01;

    // augmentation
    const D eurusd = 1.2 * (1. + h0); // multiplicative shift
    const D vol = 0.15 + h1;          // additive shift

    // Pre simulation
    const auto beta = preSimulation(eurusd, vol, 100, generator, distribution);
    const D beta1 = beta(1);
    EXPECT_NEAR(dBeta1dSpotTarget, dco::derivative(beta1)[0], 1e-6);
    EXPECT_NEAR(dBeta1dVolTarget, dco::derivative(beta1)[1], 1e-6);

    // Main simulation
    const std::size_t numPaths = 10;
    EigenVector<D> regressors(numRegressors);
    D pv = 0.;

    for (std::size_t path = 0; path < numPaths; ++path) {
        const D fwd = futureSpot(eurusd, vol, generator, distribution);
        regressors << 1.0, fwd, fwd * fwd;
        const D pathPv = beta.dot(regressors);
        pv += pathPv / numPaths;
    }

    EXPECT_NEAR(pvTarget, dco::passive_value(pv), 1e-6);
    const double dPVdSpot = dco::derivative(pv)[0];
    const double dPVdVol = dco::derivative(pv)[1];
    EXPECT_NEAR(dPVdSpotTarget, dPVdSpot, 1e-6);
    EXPECT_NEAR(dPVdVolTarget, dPVdVol, 1e-6);
}

TEST(LeastSquares, Preaccumulation) {
    std::mt19937 generator(1u);
    std::normal_distribution<double> distribution(0., 1.);

    using dco_mode = dco::ga1v<double, 1>;
    using D = dco_mode::type;
    using tape_type = dco_mode::tape_t;
    using preacc = dco_mode::jacobian_preaccumulator_t;

    dco::tape_options o;
    o.set_chunk_size_in_kbyte(100);
    dco_mode::global_tape = tape_type::create(o);
    D h0 = 0.;
    D h1 = 0.;
    dco_mode::global_tape->register_variable(h0);
    dco_mode::global_tape->register_variable(h1);

    // augmentation
    const D eurusd = 1.2 * (1. + h0 * 0.01); // multiplicative shift
    const D vol = 0.15 + h1 * 0.01;          // additive shift
    const std::size_t tapeStart =
        dco::size_of(dco_mode::global_tape, tape_type::size_of_default);

    // Preaccumulate Pre simulation
    preacc jp(dco::tape(eurusd)); // record inputs
    jp.start();                   // start recording of intermediates
    auto beta =
        preSimulation(eurusd, vol, 100, generator, distribution); // compute
    const std::size_t tapePreSim =
        dco::size_of(dco_mode::global_tape, tape_type::size_of_default);
    EXPECT_GT(tapePreSim,
              1700 * tapeStart); // tape grows a lot during pre-simulation
    jp.register_output(beta(0)); // record output
    jp.register_output(beta(1)); // record output
    jp.register_output(beta(2)); // record output
    jp.finish();                 // preaccumulate
    const std::size_t tapePreAccum =
        dco::size_of(dco_mode::global_tape, tape_type::size_of_default);
    EXPECT_LT(tapePreAccum,
              tapePreSim /
                  3); // tape reduces by about 1/3 after preaccumulation

    // Main simulation
    const double pv = mainSimulationWithReplay(eurusd, vol, beta, 10, generator,
                                               distribution);

    // check results
    EXPECT_NEAR(pvTarget, pv, 1e-6);
    dco_mode::global_tape->interpret_adjoint();
    EXPECT_NEAR(dPVdSpotTarget, dco::derivative(h0)[0], 1e-6);
    EXPECT_NEAR(dPVdVolTarget, dco::derivative(h1)[0], 1e-6);
    tape_type::remove(dco_mode::global_tape);
}

TEST(LeastSquares, SymbolicAdjoints) {
    std::mt19937 generator(1u);
    std::normal_distribution<double> distribution(0., 1.);

    using dco_mode = dco::ga1v<double, 1>;
    using D = dco_mode::type;
    using tape_type = dco_mode::tape_t;

    dco::tape_options o;
    o.set_chunk_size_in_kbyte(100);
    dco_mode::global_tape = tape_type::create(o);
    D h0 = 0.;
    D h1 = 0.;
    dco_mode::global_tape->register_variable(h0);
    dco_mode::global_tape->register_variable(h1);

    // augmentation
    const D eurusd = 1.2 * (1. + h0 * 0.01); // multiplicative shift
    const D vol = 0.15 + h1 * 0.01;          // additive shift
    const std::size_t tapeStart =
        dco::size_of(dco_mode::global_tape, tape_type::size_of_default);

    auto beta =
        preSimulationSymbolicAdjoint(eurusd, vol, 100, generator, distribution);
    const std::size_t tapePreSim =
        dco::size_of(dco_mode::global_tape, tape_type::size_of_default);
    EXPECT_LT(tapePreSim,
              875 *
                  tapeStart); // tape remains pretty small during pre-simulation

    // Main simulation
    const double pv = mainSimulationWithReplay(eurusd, vol, beta, 10, generator,
                                               distribution);

    // check results
    EXPECT_NEAR(pvTarget, pv, 1e-6);
    dco_mode::global_tape->interpret_adjoint();
    EXPECT_NEAR(dPVdSpotTarget, dco::derivative(h0)[0], 1e-6);
    EXPECT_NEAR(dPVdVolTarget, dco::derivative(h1)[0], 1e-6);
    tape_type::remove(dco_mode::global_tape);
}

TEST(LeastSquares, Checkpointing) {
    std::mt19937 generator(1u);
    std::normal_distribution<double> distribution(0., 1.);

    using dco_mode = dco::ga1v<double, 1>;
    using D = dco_mode::type;
    using tape_type = dco_mode::tape_t;
    using EAO = dco_mode::external_adjoint_object_t;

    dco::tape_options o;
    o.set_chunk_size_in_kbyte(100);
    dco_mode::global_tape = tape_type::create(o);
    D h0 = 0.;
    D h1 = 0.;
    dco_mode::global_tape->register_variable(h0);
    dco_mode::global_tape->register_variable(h1);
    const std::size_t numPresimPaths = 100;

    // augmentation
    const D eurusd = 1.2 * (1. + h0 * 0.01); // multiplicative shift
    const D vol = 0.15 + h1 * 0.01;          // additive shift
    const std::size_t tapeStart =
        dco::size_of(dco_mode::global_tape, tape_type::size_of_default);

    // create checkpoint
    EAO *eao = dco_mode::global_tape->template create_callback_object<EAO>();
    eao->register_input(eurusd);
    eao->register_input(vol);
    eao->write_data(dco::value(eurusd));
    eao->write_data(dco::value(vol));
    eao->write_data(numPresimPaths);
    dco_mode::global_tape->switch_to_passive();
    // do not record heavy computation on tape now
    auto beta =
        preSimulation(eurusd, vol, numPresimPaths, generator, distribution);
    dco_mode::global_tape->switch_to_active();
    const std::size_t tapePreSim =
        dco::size_of(dco_mode::global_tape, tape_type::size_of_default);
    EXPECT_EQ(tapePreSim, tapeStart); // tape is not used in presimulation

    for (std::size_t i = 0; i < numRegressors; ++i) {
        beta(i) = eao->register_output(dco::value(beta(i)));
    }

    dco_mode::global_tape->template insert_callback<EAO>(
        adjointCheckpoint<double, EAO>, eao);

    // Main simulation
    const double pv = mainSimulationWithReplay(eurusd, vol, beta, 10, generator,
                                               distribution);

    // check results
    EXPECT_NEAR(pvTarget, pv, 1e-6);
    dco_mode::global_tape->interpret_adjoint();
    EXPECT_NEAR(dPVdSpotTarget, dco::derivative(h0)[0], 1e-6);
    EXPECT_NEAR(dPVdVolTarget, dco::derivative(h1)[0], 1e-6);
    tape_type::remove(dco_mode::global_tape);
}

TEST(LeastSquares, CheckpointingWithSymbolicAdjoints) {
    std::mt19937 generator(1u);
    std::normal_distribution<double> distribution(0., 1.);

    using dco_mode = dco::ga1v<double, 1>;
    using D = dco_mode::type;
    using tape_type = dco_mode::tape_t;

    dco::tape_options o;
    o.set_chunk_size_in_kbyte(100);
    dco_mode::global_tape = tape_type::create(o);
    D h0 = 0.;
    D h1 = 0.;
    dco_mode::global_tape->register_variable(h0);
    dco_mode::global_tape->register_variable(h1);
    const std::size_t numPresimPaths = 100;

    // augmentation
    const D eurusd = 1.2 * (1. + h0 * 0.01); // multiplicative shift
    const D vol = 0.15 + h1 * 0.01;          // additive shift
    const std::size_t tapeStart =
        dco::size_of(dco_mode::global_tape, tape_type::size_of_default);

    // fancy checkpointing
    auto beta = preSimulationFancy(eurusd, vol, numPresimPaths, generator,
                                   distribution);
    const std::size_t tapePreSim =
        dco::size_of(dco_mode::global_tape, tape_type::size_of_default);
    EXPECT_LE(tapePreSim, 3 * tapeStart); // a tiny bit of data is written on
                                          // the tape during presimulation

    // Main simulation
    const double pv = mainSimulationWithReplay(eurusd, vol, beta, 10, generator,
                                               distribution);

    // check results
    EXPECT_NEAR(pvTarget, pv, 1e-6);
    dco_mode::global_tape->interpret_adjoint();
    EXPECT_NEAR(dPVdSpotTarget, dco::derivative(h0)[0], 1e-6);
    EXPECT_NEAR(dPVdVolTarget, dco::derivative(h1)[0], 1e-6);
    tape_type::remove(dco_mode::global_tape);
}

} // namespace BloombergLP::qfd::test

namespace adhoc {}
