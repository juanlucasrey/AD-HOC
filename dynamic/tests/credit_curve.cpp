#include <test_simple_include.hpp>

#include "credit/dates.hpp"
#include "credit/fi_instruments.hpp"
#include "credit/mdspan.hpp"
#include "credit/yield_curve.hpp"
#include "data/load.hpp"
#include <tape.hpp>

#include <chrono>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <vector>

using namespace std::chrono;

std::vector<double> cdsSpreadVals;
std::vector<double> survivalProbVals;
std::vector<double> expectedRates;
std::vector<double> expectedFirstDerivativesVals;
adhoc::mdspan<const double, 2> expectedFirstDerivatives;
std::vector<double> expectedSecondDerivativesVals;
adhoc::mdspan<const double, 3> expectedSecondDerivatives;
std::vector<double> expectedSecondDerivativesADHOCVals;
adhoc::mdspan<const double, 3> expectedSecondDerivativesADHOC;
std::chrono::year_month_day pricingDate;
std::chrono::year_month_day settleDate;
std::vector<std::chrono::year_month_day> discountCurveDates;
std::vector<double> discountCurveParams;

std::vector<double> expectedFirstDerivativesDFsVals;
adhoc::mdspan<const double, 2> expectedFirstDerivativesDFs;

template<class Float = double>
auto
get_yc(std::vector<Float> discountCurveParams_in) -> adhoc::yield_curve<Float>
{
    return adhoc::yield_curve<Float>(pricingDate, discountCurveParams_in, discountCurveDates);
}

struct CDS {
    std::vector<std::chrono::year_month_day> schedule;
    adhoc::DayCountConvention dcc;
    std::map<std::chrono::year_month_day, std::chrono::year_month_day> paymentDateMapper;
};

auto
get_cds() -> std::vector<CDS>
{
    const adhoc::DayCountConvention dcc_v = adhoc::DayCountConvention::e_ACTUAL_360;
    std::vector<CDS> cdsInstruments;
    cdsInstruments.push_back(CDS{ .schedule = { 2016y / 6 / 20, 2016y / 9 / 20, 2016y / 12 / 20 }, .dcc = dcc_v });
    cdsInstruments.push_back(CDS{
      .schedule = { 2016y / 6 / 20, 2016y / 9 / 20, 2016y / 12 / 20, 2017y / 3 / 20, 2017y / 6 / 20 }, .dcc = dcc_v });
    cdsInstruments.push_back(CDS{ .schedule = { 2016y / 6 / 20,
                                                2016y / 9 / 20,
                                                2016y / 12 / 20,
                                                2017y / 3 / 20,
                                                2017y / 6 / 20,
                                                2017y / 9 / 20,
                                                2017y / 12 / 20,
                                                2018y / 3 / 20,
                                                2018y / 6 / 20 },
                                  .dcc = dcc_v });
    cdsInstruments.push_back(CDS{ .schedule = { 2016y / 6 / 20,
                                                2016y / 9 / 20,
                                                2016y / 12 / 20,
                                                2017y / 3 / 20,
                                                2017y / 6 / 20,
                                                2017y / 9 / 20,
                                                2017y / 12 / 20,
                                                2018y / 3 / 20,
                                                2018y / 6 / 20,
                                                2018y / 9 / 20,
                                                2018y / 12 / 20,
                                                2019y / 3 / 20,
                                                2019y / 6 / 20 },
                                  .dcc = dcc_v });
    cdsInstruments.push_back(CDS{ .schedule = { 2016y / 6 / 20,
                                                2016y / 9 / 20,
                                                2016y / 12 / 20,
                                                2017y / 3 / 20,
                                                2017y / 6 / 20,
                                                2017y / 9 / 20,
                                                2017y / 12 / 20,
                                                2018y / 3 / 20,
                                                2018y / 6 / 20,
                                                2018y / 9 / 20,
                                                2018y / 12 / 20,
                                                2019y / 3 / 20,
                                                2019y / 6 / 20,
                                                2019y / 9 / 20,
                                                2019y / 12 / 20,
                                                2020y / 3 / 20,
                                                2020y / 6 / 20 },
                                  .dcc = dcc_v });
    cdsInstruments.back().paymentDateMapper = { { 2020y / 6 / 20, 2020y / 6 / 22 } };
    cdsInstruments.push_back(CDS{
      .schedule = { 2016y / 6 / 20,  2016y / 9 / 20, 2016y / 12 / 20, 2017y / 3 / 20, 2017y / 6 / 20,  2017y / 9 / 20,
                    2017y / 12 / 20, 2018y / 3 / 20, 2018y / 6 / 20,  2018y / 9 / 20, 2018y / 12 / 20, 2019y / 3 / 20,
                    2019y / 6 / 20,  2019y / 9 / 20, 2019y / 12 / 20, 2020y / 3 / 20, 2020y / 6 / 22,  2020y / 9 / 21,
                    2020y / 12 / 21, 2021y / 3 / 22, 2021y / 6 / 20 },
      .dcc = dcc_v });
    cdsInstruments.back().paymentDateMapper = { { 2021y / 6 / 20, 2021y / 6 / 21 } };
    cdsInstruments.push_back(CDS{
      .schedule = { 2016y / 6 / 20,  2016y / 9 / 20, 2016y / 12 / 20, 2017y / 3 / 20, 2017y / 6 / 20,  2017y / 9 / 20,
                    2017y / 12 / 20, 2018y / 3 / 20, 2018y / 6 / 20,  2018y / 9 / 20, 2018y / 12 / 20, 2019y / 3 / 20,
                    2019y / 6 / 20,  2019y / 9 / 20, 2019y / 12 / 20, 2020y / 3 / 20, 2020y / 6 / 22,  2020y / 9 / 21,
                    2020y / 12 / 21, 2021y / 3 / 22, 2021y / 6 / 21,  2021y / 9 / 20, 2021y / 12 / 20, 2022y / 3 / 21,
                    2022y / 6 / 20,  2022y / 9 / 20, 2022y / 12 / 20, 2023y / 3 / 20, 2023y / 6 / 20 },
      .dcc = dcc_v });
    cdsInstruments.push_back(CDS{
      .schedule = { 2016y / 6 / 20,  2016y / 9 / 20, 2016y / 12 / 20, 2017y / 3 / 20, 2017y / 6 / 20,  2017y / 9 / 20,
                    2017y / 12 / 20, 2018y / 3 / 20, 2018y / 6 / 20,  2018y / 9 / 20, 2018y / 12 / 20, 2019y / 3 / 20,
                    2019y / 6 / 20,  2019y / 9 / 20, 2019y / 12 / 20, 2020y / 3 / 20, 2020y / 6 / 22,  2020y / 9 / 21,
                    2020y / 12 / 21, 2021y / 3 / 22, 2021y / 6 / 21,  2021y / 9 / 20, 2021y / 12 / 20, 2022y / 3 / 21,
                    2022y / 6 / 20,  2022y / 9 / 20, 2022y / 12 / 20, 2023y / 3 / 20, 2023y / 6 / 20,  2023y / 9 / 20,
                    2023y / 12 / 20, 2024y / 3 / 20, 2024y / 6 / 20,  2024y / 9 / 20, 2024y / 12 / 20, 2025y / 3 / 20,
                    2025y / 6 / 20,  2025y / 9 / 22, 2025y / 12 / 22, 2026y / 3 / 20, 2026y / 6 / 20 },
      .dcc = dcc_v });
    cdsInstruments.back().paymentDateMapper = { { 2026y / 6 / 20, 2026y / 6 / 22 } };

    return cdsInstruments;
}

template<class D, class F>
void
zbrac(F func, D& x1, D& f1, D& x2, D& f2)
{
    constexpr double GOLD_FACTOR = 1.6;
    constexpr int NTRY_ZBRAC = 50;

    func(x1, f1);
    func(x2, f2);

    for (int numIter = 0; numIter < NTRY_ZBRAC; ++numIter) {
        if (f1 * f2 <= 0) { // The equality here means we consider to have bracketed the
                            // root when the function evaluates to zero at the roots too
            return;
        }

        using std::abs;
        if (abs(f1) < abs(f2)) {
            D tmp = x1;
            x1 += GOLD_FACTOR * (x1 - x2);
            x2 = tmp;
            func(x1, f1);
        }
        else {
            D tmp = x2;
            x2 += GOLD_FACTOR * (x2 - x1);
            x1 = tmp;
            func(x2, f2);
        }
    }

    throw std::runtime_error("zbrac: Unable to bracket root");
}

using namespace adhoc;

template<class Float = double>
class credit_curve final {
  public:
    template<class Float2>
    explicit credit_curve(credit_curve<Float2> const& other)
      : asofdate(other.asofdate)
      , params(other.params.begin(), other.params.end())
      , hazard(other.hazard.begin(), other.hazard.end())
      , dates(other.dates.begin(), other.dates.end())
    {
    }

    explicit credit_curve(std::chrono::year_month_day const& asofdate_in)
      : asofdate(asofdate_in)
    {
        dates.push_back(asofdate_in);
        params.push_back(1.);
    }
    explicit credit_curve(std::chrono::year_month_day const& asofdate_in,
                          std::vector<Float> const& hazard_in,
                          std::vector<std::chrono::year_month_day> const& segments_in)
      : asofdate(asofdate_in)
    {
        dates.reserve(segments_in.size());
        dates.insert(end(dates), begin(segments_in), end(segments_in));

        hazard.reserve(hazard_in.size());
        hazard.insert(end(hazard), begin(hazard_in), end(hazard_in));

        params.reserve(hazard_in.size() + 1);
        params.push_back(1.);
        for (std::size_t i = 0; i < hazard_in.size(); ++i) {
            auto const daycountFactor =
              GetCoverage(adhoc::DayCountConvention::e_ACTUAL_365_FIXED, dates[i], dates[i + 1]);
            auto const result = params[i] * exp(-hazard_in[i] * daycountFactor);
            params.push_back(result);
        }
    }

    auto getDiscountFactor(std::chrono::year_month_day const& date) const -> Float
    {
        auto lower_bound_it = std::lower_bound(this->dates.begin(), this->dates.end(), date);
        auto distance = std::distance(this->dates.begin(), lower_bound_it);

        if (lower_bound_it != this->dates.end() && *lower_bound_it == date) {
            return (this->params[static_cast<std::size_t>(distance)]);
        }

        auto date_index = distance - 1;
        auto hr_index = distance - 1;

        if (hr_index == static_cast<int>(this->hazard.size())) {
            hr_index--;
        }

        auto const& df = this->params[static_cast<std::size_t>(date_index)];
        auto const& hr = this->hazard[static_cast<std::size_t>(hr_index)];
        auto const& prev_node_date = this->dates[static_cast<std::size_t>(date_index)];

        auto dcc = adhoc::DayCountConvention::e_ACTUAL_365_FIXED;
        double const extra_dt = GetCoverage(dcc, prev_node_date, date);

        using std::exp;
        return df * exp(-hr * extra_dt);
    }

    auto setLastHazardRate(Float const& value, adhoc::DayCountConvention const& dcc) -> void
    {
        if (dates.size() < 2) {
            throw std::runtime_error("credit_curve::setLastHazardRate: not enough points to set "
                                     "hazard rate");
        }
        auto const d = dates.size() - 1;

        auto const daycountFactor = GetCoverage(dcc, dates[d - 1], dates[d]);
        hazard.back() = value;
        auto result = params[d - 1] * exp(-value * daycountFactor);
        params.back() = result;
    }

    auto setLastPoint(Float const& value) -> void { params.back() = value; }

    auto addLastPoint(std::chrono::year_month_day const& date, Float const& value) -> void
    {
        // make sure date is after last date
        if (date <= dates.back()) {
            throw std::runtime_error("credit_curve::addLastPoint: date must be after last date");
        }
        dates.push_back(date);
        hazard.push_back(0.0);
        params.push_back(value);
    }

    auto getParams() /* const */ -> std::vector<Float> /* const */& { return this->params; }

    auto getHazard() /* const */ -> std::vector<Float> /* const */& { return this->hazard; }

    auto getDates() const -> std::vector<std::chrono::year_month_day> const& { return this->dates; }

  private:
    std::chrono::year_month_day asofdate;
    std::vector<Float> params{};
    std::vector<Float> hazard{};
    std::vector<std::chrono::year_month_day> dates{};

    template<class Float2>
    friend class credit_curve;
};

template<class D, class D2>
auto
protectionPV(credit_curve<D> const& sc,
             adhoc::yield_curve<D2> const& dc,
             std::chrono::year_month_day const& periodStart,
             std::chrono::year_month_day const& periodEnd,
             std::chrono::year_month_day const& periodEndPaymentDate) -> D
{
    auto const scStart = sc.getDiscountFactor(periodStart);
    auto const scEnd = sc.getDiscountFactor(periodEnd);
    auto const dfStart = dc.getDiscountFactor(periodStart);
    auto const dfEnd = dc.getDiscountFactor(periodEndPaymentDate);
    auto ht = log(scStart / scEnd);
    auto ft = log(dfStart / dfEnd);
    auto const lt = ht + ft;
    return (dfStart * scStart - dfEnd * scEnd) * ht / lt;
}

template<class D, class D2>
auto
accruedPV(credit_curve<D> const& sc,
          adhoc::yield_curve<D2> const& dc,
          std::chrono::year_month_day const& periodStart,
          std::chrono::year_month_day const& periodEnd,
          std::chrono::year_month_day const& periodEndPaymentDate,
          double const& accruedFraction,
          double const& riskyFraction) -> D
{
    auto const scStart = sc.getDiscountFactor(periodStart);
    auto const scEnd = sc.getDiscountFactor(periodEnd);
    auto const dfStart = dc.getDiscountFactor(periodStart);
    auto const dfEnd = dc.getDiscountFactor(periodEndPaymentDate);
    auto ht = log(scStart / scEnd);
    auto ft = log(dfStart / dfEnd);
    auto const lt = ht + ft;
    return (dfStart * scStart - dfEnd * scEnd * (1 + lt)) * riskyFraction * ht / (lt * lt) +
           (dfStart * scStart - dfEnd * scEnd) * accruedFraction * ht / lt;
}

template<class D, class D2>
auto
PremiumLegPv01(credit_curve<D> const& sc,
               adhoc::yield_curve<D2> const& dc,
               adhoc::DayCountConvention const& dcc,
               std::chrono::year_month_day const& pricingDate,
               std::chrono::year_month_day const& maturityDate,
               std::chrono::year_month_day const& periodStart,
               std::chrono::year_month_day const& periodEnd,
               std::chrono::year_month_day const& periodEndPaymentDate,
               D& pv)
{
    auto plusOneIfMaturity = periodEnd == maturityDate ? periodEnd + 1 : periodEnd;
    double const daycount = GetCoverage(dcc, periodStart, plusOneIfMaturity);

    D sp = sc.getDiscountFactor(periodEnd);
    auto const df = dc.getDiscountFactor(periodEndPaymentDate);
    D riskyCoupon = daycount * sp * df;
    pv += riskyCoupon;

    double accruedFraction(0);
    if (periodStart <= pricingDate) {
        accruedFraction = GetCoverage(dcc, periodStart, pricingDate + 1);
    }

    double riskyFraction = daycount - accruedFraction;

    D acc_pv = accruedPV(
      sc, dc, std::max(pricingDate, periodStart), periodEnd, periodEndPaymentDate, accruedFraction, riskyFraction);
    pv += acc_pv;
}

template<class D, class D2>
void
set_rate_and_price(credit_curve<D>& sc,
                   adhoc::yield_curve<D2> const& yc,
                   CDS const& cds,
                   D const& cds_quote,
                   std::chrono::year_month_day const& pricingDate,
                   D const& x,
                   D& f)
{
    sc.setLastHazardRate(x, adhoc::DayCountConvention::e_ACTUAL_365_FIXED);
    std::vector<std::chrono::year_month_day> risky_schedule;
    risky_schedule.push_back(pricingDate);
    for (auto const& date : cds.schedule) {
        if (date > pricingDate) {
            risky_schedule.push_back(date);
        }
    }

    double notional = 100.;
    double recovery_Rate = 0.4;
    D protectionPV_value = 0;
    for (std::size_t i = 0; i < risky_schedule.size() - 1; ++i) {
        std::chrono::year_month_day periodEndPaymentDate = risky_schedule[i + 1];
        auto it = cds.paymentDateMapper.find(risky_schedule[i + 1]);
        if (it != cds.paymentDateMapper.end()) {
            periodEndPaymentDate = it->second;
        }

        auto pvi = protectionPV(sc, yc, risky_schedule[i], risky_schedule[i + 1], periodEndPaymentDate);
        protectionPV_value += pvi;
    }
    protectionPV_value = notional * protectionPV_value * (1.0 - recovery_Rate);

    std::vector<std::chrono::year_month_day> const& cds_schedule2 = cds.schedule;

    D pv01 = 0;
    for (std::size_t i = 0; i < cds_schedule2.size() - 1; ++i) {
        std::chrono::year_month_day periodEndPaymentDate = cds_schedule2[i + 1];
        auto it = cds.paymentDateMapper.find(cds_schedule2[i + 1]);
        if (it != cds.paymentDateMapper.end()) {
            periodEndPaymentDate = it->second;
        }

        PremiumLegPv01(sc,
                       yc,
                       cds.dcc,
                       pricingDate,
                       cds_schedule2.back(),
                       cds_schedule2[i],
                       cds_schedule2[i + 1],
                       periodEndPaymentDate,
                       pv01);
    }

    pv01 *= notional;

    auto premiumPV = -pv01 * cds_quote;
    auto dirtyPV = protectionPV_value + premiumPV;
    auto const dfFromSettle = yc.getDiscountFactor(settleDate);
    auto cashOnSettle = dirtyPV / dfFromSettle;
    double accruedFraction = GetCoverage(cds.dcc, cds_schedule2.front(), pricingDate + 1);
    double accrued_pv01 = accruedFraction * notional;
    auto accruedInterest = -accrued_pv01 * cds_quote;
    D rebateInterestDiscounted = 0;

    auto upfront = cashOnSettle - accruedInterest - rebateInterestDiscounted;

    double quoteUpfrontValue = 0.0;
    f = upfront - quoteUpfrontValue;
}

auto
calibrateQuote(credit_curve<double>& sc,
               adhoc::yield_curve<double> const& yc,
               CDS const& cds,
               double cds_quote,
               std::chrono::year_month_day const& pricingDate) -> double
{
    std::chrono::year_month_day newpoint = cds.schedule.back();
    auto it = cds.paymentDateMapper.find(newpoint);
    if (it != cds.paymentDateMapper.end()) {
        newpoint = it->second;
    }
    sc.addLastPoint(newpoint, 1.0);

    double x1 = 1.0E-8;
    double f1 = 0;
    double x2 = 0.01;
    double f2 = 0;

    auto func = [&](double x, double& f) { set_rate_and_price(sc, yc, cds, cds_quote, pricingDate, x, f); };

    zbrac(func, x1, f1, x2, f2);

    double b = x2;
    double fb = f2;
    double a = x1;
    double fa = f1;

    double currentGuess = 0.5 * (a + b);

    using std::abs;
    double previousChangeInCurrentGuess = abs(b - a);

    double changeInCurrentGuess = previousChangeInCurrentGuess;
    double f, df = 0.0;
    double boundWithNegF, boundWithPosF = 0.0;

    if (fa > 0.0) // switch a and b as we fa < 0.0 and fb > 0.0
    {
        boundWithNegF = b;
        boundWithPosF = a;
    }
    else {
        boundWithNegF = a;
        boundWithPosF = b;
    }

    auto func_d = [&](double x, double& f, double& df) {
        using adhoc_mode = adhoc::opcode<double>;
        using D = adhoc_mode::type;
        D inputD = x;

        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tapeptr;
        auto& tape = *tapeptr;
        tape.register_variable(inputD);
        D outputD = 0.0;

        D cds_quote2 = cds_quote;
        credit_curve<D> sc2(sc);

        set_rate_and_price(sc2, yc, cds, cds_quote2, pricingDate, inputD, outputD);

        tape.register_output_variable(outputD);

        tape.set_derivative(outputD, 1.0);
        tape.backpropagate();
        df = tape.get_derivative(inputD);
        f = adhoc::passive_value(outputD);
    };

    func_d(currentGuess, f, df);

    double newtonRestartValue = 0.0;

    int const index = 1;
    std::vector<double> previousNewtonGuesses(static_cast<std::size_t>(index), 0.0);
    std::size_t indexToSetNext = 0;
    auto const maxIndex = static_cast<std::size_t>(index - 1);
    int const minNewtonStepsToReturn = 2; // Necessary to give correct AD Greeks up to order 2

    double tolX = 1e-9;
    double tolY = 1e-15;
    int numNewtonStepsSinceLastBisection = 0;
    int const maxTotalSteps = 102;
    bool allowBisectionSteps = true;

    const auto fMeetsTolY = [tolY](const double& f) { return std::abs(f) < tolY; };
    for (int numIterations = 1; numIterations < maxTotalSteps; numIterations++) {
        using std::abs;
        if (allowBisectionSteps &&
            (((currentGuess - boundWithNegF) * df - f) * ((currentGuess - boundWithPosF) * df - f) >
               0.0 // Next Newton guess would exceed bounds
             || (std::abs(2 * f) > std::abs(previousChangeInCurrentGuess * df)))) // Convergence isn't fast enough
        {                                                                         // Bisect
            numNewtonStepsSinceLastBisection = 0;
            previousChangeInCurrentGuess = changeInCurrentGuess;
            changeInCurrentGuess = 0.5 * (boundWithPosF - boundWithNegF);
            currentGuess = boundWithNegF + changeInCurrentGuess;

            if (std::abs(changeInCurrentGuess) < tolX) {
                allowBisectionSteps = false;
            }
        }
        else {
            // Newton step
            numNewtonStepsSinceLastBisection++;
            previousNewtonGuesses[indexToSetNext] = currentGuess;
            indexToSetNext = indexToSetNext == maxIndex ? 0 : indexToSetNext + 1;
            previousChangeInCurrentGuess = changeInCurrentGuess;
            changeInCurrentGuess = f / df;
            currentGuess -= f / df;

            if ((std::abs(changeInCurrentGuess) < tolX || fMeetsTolY(f)) &&
                numNewtonStepsSinceLastBisection > minNewtonStepsToReturn) {
                b = currentGuess;

                newtonRestartValue = previousNewtonGuesses[indexToSetNext];
                break;
            }
        }

        func_d(currentGuess, f, df);
        if (f < 0.0) // Update bounds so root \in [boundWithNegF, boundWithPosF]
        {
            boundWithNegF = currentGuess;
        }
        else {
            boundWithPosF = currentGuess;
        }
    }

    return currentGuess;
}

// template<DerivativeMethod dm = DerivativeMethod::Bwd>
auto
createCDSCurve(std::vector<CDS> const& cds,
               adhoc::yield_curve<double> const& yc,
               std::vector<double> const& cds_quote,
               std::chrono::year_month_day const& pricing_date) -> credit_curve<double>
{

    auto sc = credit_curve(pricing_date);
    for (std::size_t i = 0; i < cds.size(); ++i) {
        double calibratedRates = calibrateQuote(sc, yc, cds[i], cds_quote[i], pricing_date);
        sc.setLastHazardRate(calibratedRates, adhoc::DayCountConvention::e_ACTUAL_365_FIXED);
    }

    return sc;
}

namespace {
auto
passive_value(double value) -> double
{
    return value;
}
} // namespace

template<class D, class DYC, bool adhocOrder1 = true, bool Experimental = false>
auto
pureNewtonSolve(credit_curve<D>& sc,
                adhoc::yield_curve<DYC> const& yc,
                CDS const& cds,
                D cds_quote,
                D& currentGuess,
                std::chrono::year_month_day const& pricingDate) -> void
{

    auto func_d2 = [&](double x, double& df) {
        using D2 = adhoc::opcode<double>::type;

        auto& hazard = sc.getHazard();
        std::vector<D2> clonedHazardRates;
        clonedHazardRates.reserve(hazard.size());
        for (const auto& rate : hazard) {
            clonedHazardRates.emplace_back(passive_value(rate));
        }
        auto& dates = sc.getDates();
        credit_curve<D2> sc2(pricingDate, clonedHazardRates, dates);

        auto& paramsyc = yc.getParams();
        std::vector<DYC> clonedParamsYC;
        clonedParamsYC.reserve(paramsyc.size());
        for (const auto& rate : paramsyc) {
            clonedParamsYC.emplace_back(passive_value(rate));
        }
        adhoc::yield_curve<DYC> yc2(yc.getAsOfDate(), clonedParamsYC, yc.getDates());

        D2 inputD = x;
        D2 outputD = 0.0;
        D2 cds_quote_d2 = passive_value(cds_quote);

        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tapeptr;
        auto& tape = *tapeptr;
        tape.register_variable(inputD);

        set_rate_and_price(sc2, yc2, cds, cds_quote_d2, pricingDate, inputD, outputD);

        tape.register_output_variable(outputD);

        tape.set_derivative(outputD, 1.0);
        tape.backpropagate();
        df = tape.get_derivative(inputD);
        // f = adhoc::passive_value(outputD);
    };

    auto func_d3 = [&](double x, double& df, double& df2) {
        using D2 = adhoc::opcode<double>::type;

        auto& hazard = sc.getHazard();
        std::vector<D2> clonedHazardRates;
        clonedHazardRates.reserve(hazard.size());
        for (const auto& rate : hazard) {
            clonedHazardRates.emplace_back(passive_value(rate));
        }
        auto& dates = sc.getDates();
        credit_curve<D2> sc_d2(pricingDate, clonedHazardRates, dates);

        auto& paramsyc = yc.getParams();
        std::vector<DYC> clonedParamsYC;
        clonedParamsYC.reserve(paramsyc.size());
        for (const auto& rate : paramsyc) {
            clonedParamsYC.emplace_back(passive_value(rate));
        }
        adhoc::yield_curve<DYC> yc2(yc.getAsOfDate(), clonedParamsYC, yc.getDates());

        D2 inputD = x;
        D2 outputD = 0.0;
        D2 cds_quote_d2 = passive_value(cds_quote);

        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tapeptr;
        auto& tape = *tapeptr;
        tape.configure(adhoc::Method::SecondOrderSimple, 1, 1);

        tape.register_variable(inputD);

        set_rate_and_price(sc_d2, yc2, cds, cds_quote_d2, pricingDate, inputD, outputD);

        tape.register_output_variable(outputD);

        tape.set_derivative(outputD, 1.0);
        tape.backpropagate();
        tape.get_derivative(inputD);
        df = tape.get_derivative(inputD);
        df2 = tape.get_derivative(inputD, inputD);

        // f = adhoc::passive_value(outputD);
    };

    auto func = [&](D const& x, D& f) { set_rate_and_price(sc, yc, cds, cds_quote, pricingDate, x, f); };

    const double originalNewtonRoot = passive_value(currentGuess);

    if constexpr (std::is_same_v<D, adhoc::opcode<double>::type>) {
        if constexpr (Experimental) {
            auto& tape = cds_quote.get_tape();
            // auto tape = sc.getHazard().front().get_tape();
            tape.start_implicit_function();
            tape.register_implicit_variable(currentGuess);
            D f = 0;
            set_rate_and_price(sc, yc, cds, cds_quote, pricingDate, currentGuess, f);
            tape.register_constraint_variable(f);
            tape.end_implicit_function();
        }
        else {
            if constexpr (adhocOrder1) {
                double df = 0;
                func_d2(originalNewtonRoot, df);
                D f = 0.0;
                func(currentGuess, f);
                currentGuess -= f / df;
            }
            else {
                double df = 0, df2 = 0;
                func_d3(originalNewtonRoot, df, df2);

                D f = 0.0;
                func(currentGuess, f);
                currentGuess -= f / df;
                func(currentGuess, f);

                const D deltaInput = currentGuess - passive_value(currentGuess);
                D df_bis = df + deltaInput * df2;
                currentGuess -= f / df_bis;
            }
        }
    }
    else {
        double df = 0;
        func_d2(originalNewtonRoot, df);
        D f = 0.0;
        func(currentGuess, f);
        currentGuess -= f / df;
    }

    if constexpr (!Experimental) {
        passive_value(currentGuess) = originalNewtonRoot;
    }

    sc.setLastHazardRate(currentGuess, adhoc::DayCountConvention::e_ACTUAL_365_FIXED);
}

void
normal_calib()
{
    std::vector<CDS> cdsInstruments = get_cds();
    const auto yc = get_yc(discountCurveParams);

    std::vector<double> cdsSpreads(cdsSpreadVals.begin(), cdsSpreadVals.end());

    auto sc = createCDSCurve(cdsInstruments, yc, cdsSpreads, pricingDate);
    auto const& rates = sc.getHazard();
    for (std::size_t i = 0; i < rates.size(); ++i) {
        EXPECT_NEAR_ABS(expectedRates[i], rates[i], 1e-15);
    }

    auto const& params = sc.getParams();
    for (std::size_t i = 0; i < rates.size(); ++i) {
        EXPECT_NEAR_ABS(survivalProbVals[i], params[i], 1e-15);
    }
}

void
fd1()
{
    double epsilon = 1e-6;
    std::vector<CDS> cdsInstruments = get_cds();
    const auto yc = get_yc(discountCurveParams);
    std::vector<double> cdsSpreads(cdsSpreadVals.begin(), cdsSpreadVals.end());

    auto sc = createCDSCurve(cdsInstruments, yc, cdsSpreads, pricingDate);

    for (std::size_t j = 0; j < cdsSpreads.size(); ++j) {
        double originalCDSSpread = cdsSpreads[j];
        cdsSpreads[j] = originalCDSSpread + epsilon;
        auto sc_up = createCDSCurve(cdsInstruments, yc, cdsSpreads, pricingDate);

        cdsSpreads[j] = originalCDSSpread - epsilon;
        auto sc_dn = createCDSCurve(cdsInstruments, yc, cdsSpreads, pricingDate);

        cdsSpreads[j] = originalCDSSpread;

        auto const& sp_up = sc_up.getParams();
        auto const& sp_dn = sc_dn.getParams();
        for (std::size_t i = 1; i < sp_up.size(); ++i) {
            const auto fdSensitivity = (sp_up[i] - sp_dn[i]) / (2.0 * epsilon);
            EXPECT_NEAR_ABS(expectedFirstDerivatives(i - 1, j), fdSensitivity, 2e-7);
        }
    }
}

void
create_curve_adhoc()
{
    using adhoc_mode = adhoc::opcode<double>;
    using D = adhoc_mode::type;
    adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tapeptr;
    auto& tape = *tapeptr;

    std::vector<CDS> cdsInstruments = get_cds();
    const auto yc = get_yc(discountCurveParams);
    std::vector<double> cdsSpreads_double(cdsSpreadVals.begin(), cdsSpreadVals.end());

    auto sc_double = createCDSCurve(cdsInstruments, yc, cdsSpreads_double, pricingDate);

    auto const& rates = sc_double.getHazard();
    for (std::size_t i = 0; i < rates.size(); ++i) {
        EXPECT_NEAR_ABS(expectedRates[i], rates[i], 1e-15);
    }

    auto const& params = sc_double.getParams();
    for (std::size_t i = 0; i < rates.size(); ++i) {
        EXPECT_NEAR_ABS(survivalProbVals[i], params[i], 1e-15);
    }
}

void
bwd_adhoc_df()
{
    using adhoc_mode = adhoc::opcode<double>;
    using D = adhoc_mode::type;
    adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tapeptr;
    auto& tape = *tapeptr;

    std::vector<CDS> cdsInstruments = get_cds();
    std::vector<D> discountCurveParamsD{ discountCurveParams.begin(), discountCurveParams.end() };
    tape.register_variable(discountCurveParamsD);

    const auto yc = get_yc(discountCurveParams);
    const auto ycD = get_yc(discountCurveParamsD);
    std::vector<double> cdsSpreads_double(cdsSpreadVals.begin(), cdsSpreadVals.end());

    auto sc_double = createCDSCurve(cdsInstruments, yc, cdsSpreads_double, pricingDate);
    auto const& calibratedRates = sc_double.getHazard();

    std::vector<D> cdsSpreads(cdsSpreadVals.begin(), cdsSpreadVals.end());
    // const auto initialTapePosition = tape.get_position();
    tape.register_variable(cdsSpreads);

    auto sc = credit_curve<D>(pricingDate);

    for (std::size_t j = 0; j < cdsInstruments.size(); ++j) {
        D calibratedRate = calibratedRates[j];
        std::chrono::year_month_day newpoint = cdsInstruments[j].schedule.back();
        auto it = cdsInstruments[j].paymentDateMapper.find(newpoint);
        if (it != cdsInstruments[j].paymentDateMapper.end()) {
            newpoint = it->second;
        }
        sc.addLastPoint(newpoint, 1.0);

        pureNewtonSolve<D>(sc, ycD, cdsInstruments[j], cdsSpreads[j], calibratedRate, pricingDate);
    }

    auto& y = sc.getParams();
    tape.register_output_variable(y);

    for (std::size_t j = 1; j < y.size(); ++j) {
        EXPECT_NEAR_ABS(survivalProbVals[j], y[j].get_value(), 2e-15);

        tape.set_derivative(y[j], 1.0);
        tape.backpropagate();

        for (std::size_t i = 0; i < cdsSpreads.size(); ++i) {
            auto const der = tape.get_derivative(cdsSpreads[i]);
            auto const expexted_der = expectedFirstDerivatives(j - 1, i);
            EXPECT_NEAR_ABS(expexted_der, der, 1e-10);
        }

        for (std::size_t i = 0; i < discountCurveParamsD.size(); ++i) {
            auto const der = tape.get_derivative(discountCurveParamsD[i]);
            auto const expexted_der = expectedFirstDerivativesDFs(j - 1, i);
            EXPECT_NEAR_ABS(expexted_der, der, 1e-10);
        }

        tape.set_derivative(y[j], 0.0);
        tape.zero_adjoints();
    }
}

void
bwd_adhoc()
{
    using adhoc_mode = adhoc::opcode<double>;
    using D = adhoc_mode::type;
    adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tapeptr;
    auto& tape = *tapeptr;

    std::vector<CDS> cdsInstruments = get_cds();
    const auto yc = get_yc(discountCurveParams);
    std::vector<double> cdsSpreads_double(cdsSpreadVals.begin(), cdsSpreadVals.end());

    auto sc_double = createCDSCurve(cdsInstruments, yc, cdsSpreads_double, pricingDate);
    auto const& calibratedRates = sc_double.getHazard();

    std::vector<D> cdsSpreads(cdsSpreadVals.begin(), cdsSpreadVals.end());
    // const auto initialTapePosition = tape.get_position();
    tape.register_variable(cdsSpreads);

    auto sc = credit_curve<D>(pricingDate);

    for (std::size_t j = 0; j < cdsInstruments.size(); ++j) {
        D calibratedRate = calibratedRates[j];
        std::chrono::year_month_day newpoint = cdsInstruments[j].schedule.back();
        auto it = cdsInstruments[j].paymentDateMapper.find(newpoint);
        if (it != cdsInstruments[j].paymentDateMapper.end()) {
            newpoint = it->second;
        }
        sc.addLastPoint(newpoint, 1.0);

        pureNewtonSolve<D>(sc, yc, cdsInstruments[j], cdsSpreads[j], calibratedRate, pricingDate);
    }

    auto& y = sc.getParams();
    tape.register_output_variable(y);

    for (std::size_t j = 1; j < y.size(); ++j) {
        EXPECT_NEAR_ABS(survivalProbVals[j], y[j].get_value(), 2e-15);

        tape.set_derivative(y[j], 1.0);
        tape.backpropagate();

        for (std::size_t i = 0; i < cdsSpreads.size(); ++i) {
            auto const der = tape.get_derivative(cdsSpreads[i]);
            auto const expexted_der = expectedFirstDerivatives(j - 1, i);
            EXPECT_NEAR_ABS(expexted_der, der, 1e-10);
        }

        tape.set_derivative(y[j], 0.0);
        tape.zero_adjoints();
    }
}

void
bwd_adhoc2()
{
    using adhoc_mode = adhoc::opcode<double>;
    using D = adhoc_mode::type;
    adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tapeptr;
    auto& tape = *tapeptr;
    tape.configure(adhoc::Method::SecondOrderSimple, 1, 1);

    std::vector<CDS> cdsInstruments = get_cds();
    const auto yc = get_yc(discountCurveParams);
    std::vector<double> cdsSpreads_double(cdsSpreadVals.begin(), cdsSpreadVals.end());

    auto sc_double = createCDSCurve(cdsInstruments, yc, cdsSpreads_double, pricingDate);
    auto const& calibratedRates = sc_double.getHazard();

    std::vector<D> cdsSpreads(cdsSpreadVals.begin(), cdsSpreadVals.end());
    tape.register_variable(cdsSpreads);

    auto sc = credit_curve<D>(pricingDate);

    for (std::size_t j = 0; j < cdsInstruments.size(); ++j) {
        D calibratedRate = calibratedRates[j];
        std::chrono::year_month_day newpoint = cdsInstruments[j].schedule.back();
        auto it = cdsInstruments[j].paymentDateMapper.find(newpoint);
        if (it != cdsInstruments[j].paymentDateMapper.end()) {
            newpoint = it->second;
        }
        sc.addLastPoint(newpoint, 1.0);

        pureNewtonSolve<D, double, false>(sc, yc, cdsInstruments[j], cdsSpreads[j], calibratedRate, pricingDate);
    }

    auto& y = sc.getParams();
    tape.register_output_variable(y);

    for (std::size_t j = 1; j < y.size(); ++j) {
        EXPECT_NEAR_ABS(survivalProbVals[j], y[j].get_value(), 2e-15);

        tape.set_derivative(y[j], 1.0);
        tape.backpropagate();

        for (std::size_t i = 0; i < cdsSpreads.size(); ++i) {
            auto const der = tape.get_derivative(cdsSpreads[i]);
            auto const expexted_der = expectedFirstDerivatives(j - 1, i);
            EXPECT_NEAR_ABS(expexted_der, der, 1e-10);
        }

        for (std::size_t i = 0; i < cdsSpreads.size(); ++i) {
            for (std::size_t k = 0; k < cdsSpreads.size(); ++k) {
                auto const der = tape.get_derivative(cdsSpreads[i], cdsSpreads[k]);
                auto const expexted_der = expectedSecondDerivativesADHOC(k, i, j - 1);
                EXPECT_NEAR_ABS(expexted_der, der, 1e-10);
            }
        }

        tape.set_derivative(y[j], 0.0);
        tape.zero_adjoints();
    }
}

void
bwd_adhoc2_payoff()
{

    // common data
    std::vector<CDS> cdsInstruments = get_cds();
    const auto yc = get_yc(discountCurveParams);

    auto sc_double = createCDSCurve(cdsInstruments, yc, cdsSpreadVals, pricingDate);
    auto const& calibratedRates = sc_double.getHazard();

    // adhoc results
    std::vector<double> first_derivs_adhoc(cdsSpreadVals.size());
    std::vector<double> second_derivs_vals_adhoc(cdsSpreadVals.size() * cdsSpreadVals.size());
    adhoc::mdspan<double, 2> second_derivs_adhoc(
      second_derivs_vals_adhoc.data(), cdsSpreadVals.size(), cdsSpreadVals.size());

    {
        using adhoc_mode = adhoc::opcode<double>;
        using D = adhoc_mode::type;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tapeptr;
        auto& tape = *tapeptr;
        tape.configure(adhoc::Method::SecondOrderSimple, 1, 1);

        std::vector<D> cdsSpreads(cdsSpreadVals.begin(), cdsSpreadVals.end());
        tape.register_variable(cdsSpreads);

        auto sc = credit_curve<D>(pricingDate);

        for (std::size_t j = 0; j < cdsInstruments.size(); ++j) {
            D calibratedRate = calibratedRates[j];
            std::chrono::year_month_day newpoint = cdsInstruments[j].schedule.back();
            auto it = cdsInstruments[j].paymentDateMapper.find(newpoint);
            if (it != cdsInstruments[j].paymentDateMapper.end()) {
                newpoint = it->second;
            }
            sc.addLastPoint(newpoint, 1.0);

            pureNewtonSolve<D, double, false>(sc, yc, cdsInstruments[j], cdsSpreads[j], calibratedRate, pricingDate);
        }

        // write a random payoff using sc's getDiscountFactor
        std::vector<std::chrono::year_month_day> payoff_dates{ 2024y / 6 / 15, 2025y / 6 / 15, 2026y / 6 / 15 };
        std::vector<D> payoff_amounts{ 1000.0, 1000.0, 1000.0 };
        D total_payoff = 0.0;
        for (std::size_t i = 0; i < payoff_dates.size(); ++i) {
            D df = sc.getDiscountFactor(payoff_dates[i]);
            total_payoff += df * payoff_amounts[i];
        }

        tape.register_output_variable(total_payoff);

        tape.set_derivative(total_payoff, 1.0);
        tape.backpropagate();
        for (std::size_t i = 0; i < cdsSpreads.size(); ++i) {
            first_derivs_adhoc[i] = tape.get_derivative(cdsSpreads[i]);
        }

        for (std::size_t i = 0; i < cdsSpreads.size(); ++i) {
            for (std::size_t k = 0; k < cdsSpreads.size(); ++k) {
                second_derivs_adhoc(i, k) = tape.get_derivative(cdsSpreads[i], cdsSpreads[k]);
            }
        }
    }

    std::vector<double> first_derivs_res(cdsSpreadVals.size());
    std::vector<double> second_derivs_vals_res(cdsSpreadVals.size() * cdsSpreadVals.size());
    adhoc::mdspan<double, 2> second_derivs_res(
      second_derivs_vals_res.data(), cdsSpreadVals.size(), cdsSpreadVals.size());

    // compare first and second order derivatives
    for (std::size_t i = 0; i < cdsSpreadVals.size(); ++i) {
        EXPECT_NEAR_ABS(first_derivs_adhoc[i], first_derivs_res[i], 1e-10);
        for (std::size_t j = 0; j < cdsSpreadVals.size(); ++j) {
            EXPECT_NEAR_ABS(second_derivs_adhoc(i, j), second_derivs_res(i, j), 2e-7);
        }
    }
}

void
compressed_test()
{
    // cds?
    double x = 0.0050000050000000001;
    double cds_quote = 0.0025571700000000001;
    auto sc = credit_curve(pricingDate);
    sc.addLastPoint(2016y / 12 / 22, 1.0);

    const auto yc = get_yc(discountCurveParams);
    std::vector<CDS> cdsInstruments = get_cds();
    auto& cds = cdsInstruments[0];
    double df1 = 0;
    double f1 = 0;

    double df2 = 0;
    double f2 = 0;

    using adhoc_mode = adhoc::opcode<double>;
    using D = adhoc_mode::type;

    {
        D inputD = x;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tapeptr;
        tapeptr->configure(adhoc::Method::BwdBuffer, 1, 1);
        auto& tape = *tapeptr;
        tape.register_variable(inputD);
        D outputD = 0.0;

        D cds_quote2 = cds_quote;
        credit_curve<D> sc2(sc);

        set_rate_and_price(sc2, yc, cds, cds_quote2, pricingDate, inputD, outputD);

        tape.register_output_variable(outputD);

        tape.set_derivative(outputD, 1.0);
        tape.backpropagate();
        df1 = tape.get_derivative(inputD);
        f1 = adhoc::passive_value(outputD);
    }

    {
        D inputD = x;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tapeptr;
        tapeptr->configure(adhoc::Method::BwdBufferCompressed, 1, 1);
        auto& tape = *tapeptr;
        tape.register_variable(inputD);
        D outputD = 0.0;

        D cds_quote2 = cds_quote;
        credit_curve<D> sc2(sc);

        set_rate_and_price(sc2, yc, cds, cds_quote2, pricingDate, inputD, outputD);

        tape.register_output_variable(outputD);

        tape.set_derivative(outputD, 1.0);
        tape.backpropagate();
        df2 = tape.get_derivative(inputD);
        f2 = adhoc::passive_value(outputD);
    }

    EXPECT_NEAR_ABS(df1, df2, 2e-14);
    EXPECT_NEAR_ABS(f1, f2, 2e-15);
}

void
lossy_test()
{
    // cds?
    double x1 = 0.0050000050000000001;
    double x2 = 0.0043152164958849123;
    double cds_quote = 0.0025571700000000001;
    auto sc = credit_curve(pricingDate);
    sc.addLastPoint(2016y / 12 / 22, 1.0);

    const auto yc = get_yc(discountCurveParams);
    std::vector<CDS> cdsInstruments = get_cds();
    auto& cds = cdsInstruments[0];
    double df11 = 0;
    double f11 = 0;
    double df12 = 0;
    double f12 = 0;

    double df21 = 0;
    double f21 = 0;
    double df22 = 0;
    double f22 = 0;

    using adhoc_mode = adhoc::opcode<double>;
    using D = adhoc_mode::type;

    {
        D inputD = x1;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tapeptr;
        tapeptr->configure(adhoc::Method::Bwd, 1, 1);
        auto& tape = *tapeptr;
        tape.register_variable(inputD);
        D outputD = 0.0;

        D cds_quote2 = cds_quote;
        credit_curve<D> sc2(sc);

        set_rate_and_price(sc2, yc, cds, cds_quote2, pricingDate, inputD, outputD);

        tape.register_output_variable(outputD);

        tape.set_derivative(outputD, 1.0);
        tape.backpropagate();
        df11 = tape.get_derivative(inputD);
        f11 = adhoc::passive_value(outputD);

        D inputD2 = x2;
        tape.register_variable(inputD2);
        D outputD2 = 0.0;

        set_rate_and_price(sc2, yc, cds, cds_quote2, pricingDate, inputD2, outputD2);

        tape.register_output_variable(outputD2);

        tape.set_derivative(outputD2, 1.0);
        tape.backpropagate();
        df12 = tape.get_derivative(inputD2);
        f12 = adhoc::passive_value(outputD2);
    }

    {
        D inputD = x1;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tapeptr;
        tapeptr->configure(adhoc::Method::BwdBuffer, 1, 1);
        auto& tape = *tapeptr;
        tape.register_variable(inputD);
        D outputD = 0.0;

        D cds_quote2 = cds_quote;
        credit_curve<D> sc2(sc);

        set_rate_and_price(sc2, yc, cds, cds_quote2, pricingDate, inputD, outputD);

        tape.register_output_variable(outputD);

        tape.set_derivative(outputD, 1.0);
        tape.backpropagate();
        df21 = tape.get_derivative(inputD);
        f21 = adhoc::passive_value(outputD);

        D inputD2 = x2;
        tape.register_variable(inputD2);
        D outputD2 = 0.0;

        set_rate_and_price(sc2, yc, cds, cds_quote2, pricingDate, inputD2, outputD2);

        tape.register_output_variable(outputD2);

        tape.set_derivative(outputD2, 1.0);
        tape.backpropagate();
        df22 = tape.get_derivative(inputD2);
        f22 = adhoc::passive_value(outputD2);
    }

    EXPECT_NEAR_ABS(df11, df21, 2e-14);
    EXPECT_NEAR_ABS(f11, f21, 2e-15);
    EXPECT_NEAR_ABS(df12, df22, 2e-14);
    EXPECT_NEAR_ABS(f12, f22, 2e-15);
}

void
lossy_reuse_test()
{
    // cds?
    double x = 0.0050000050000000001;
    double cds_quote = 0.0025571700000000001;
    auto sc = credit_curve(pricingDate);
    sc.addLastPoint(2016y / 12 / 22, 1.0);

    const auto yc = get_yc(discountCurveParams);
    std::vector<CDS> cdsInstruments = get_cds();
    auto& cds = cdsInstruments[0];
    double df1 = 0;
    double f1 = 0;

    double df2 = 0;
    double f2 = 0;

    using adhoc_mode = adhoc::opcode<double>;
    using D = adhoc_mode::type;

    {
        D inputD = x;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tapeptr;
        tapeptr->configure(adhoc::Method::BwdBuffer, 1, 1);
        auto& tape = *tapeptr;
        tape.register_variable(inputD);
        D outputD = 0.0;

        D cds_quote2 = cds_quote;
        credit_curve<D> sc2(sc);

        set_rate_and_price(sc2, yc, cds, cds_quote2, pricingDate, inputD, outputD);

        tape.register_output_variable(outputD);

        tape.set_derivative(outputD, 1.0);
        tape.backpropagate();
        df1 = tape.get_derivative(inputD);
        f1 = adhoc::passive_value(outputD);
    }

    {
        D inputD = x;
        adhoc::smart_tape_ptr_t<adhoc::opcode<double> > tapeptr;
        tapeptr->configure(adhoc::Method::BwdBufferPathReuse, 1, 1);
        auto& tape = *tapeptr;
        tape.register_variable(inputD);
        D outputD = 0.0;

        D cds_quote2 = cds_quote;
        credit_curve<D> sc2(sc);

        set_rate_and_price(sc2, yc, cds, cds_quote2, pricingDate, inputD, outputD);

        tape.register_output_variable(outputD);

        tape.set_derivative(outputD, 1.0);
        tape.backpropagate();
        df2 = tape.get_derivative(inputD);
        f2 = adhoc::passive_value(outputD);
    }

    EXPECT_NEAR_ABS(df1, df2, 2e-14);
    EXPECT_NEAR_ABS(f1, f2, 2e-15);
}

auto
main() -> int
{
    // init vectors
    cdsSpreadVals = loadCSVToVector("cds_spread_vals.csv");
    survivalProbVals = loadCSVToVector("survival_prob_vals.csv");
    expectedRates = loadCSVToVector("expected_rates.csv");
    {
        auto [vals, span] = loadCSVToTensor("expected_first_derivatives_vals.csv", 8, 8);
        expectedFirstDerivativesVals = std::move(vals);
        expectedFirstDerivatives = span;
    }

    {
        auto [vals, span] = loadCSVToTensor("expected_second_derivatives_vals.csv", 8, 8, 8);
        expectedSecondDerivativesVals = std::move(vals);
        expectedSecondDerivatives = span;
    }

    {
        auto [vals, span] = loadCSVToTensor("expected_second_derivatives_ADHOC_vals.csv", 8, 8, 8);
        expectedSecondDerivativesADHOCVals = std::move(vals);
        expectedSecondDerivativesADHOC = span;
    }

    {
        auto [vals, span] = loadCSVToTensor("expected_first_derivatives_dfs_vals.csv", 8, 8);
        expectedFirstDerivativesDFsVals = std::move(vals);
        expectedFirstDerivativesDFs = span;
    }

    pricingDate = loadCSVToVectorDate("pricing_date.csv").front();
    settleDate = loadCSVToVectorDate("settle_date.csv").front();
    discountCurveDates = loadCSVToVectorDate("discount_curve_dates.csv");
    discountCurveParams = loadCSVToVector("discount_curve_params.csv");

    bwd_adhoc_df();

    lossy_reuse_test();

    lossy_test();
    compressed_test();
    normal_calib();
    fd1();
    create_curve_adhoc();
    bwd_adhoc();
    bwd_adhoc2();
    // bwd_adhoc2_payoff();

    TEST_END;
}
