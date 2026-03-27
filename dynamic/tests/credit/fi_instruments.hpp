/*
 * AD-HOC, Automatic Differentiation for High Order Calculations
 *
 * This file is part of the AD-HOC distribution
 * (https://github.com/juanlucasrey/AD-HOC).
 * Copyright (c) 2026 Juan Lucas Rey
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef BMC_FI_INSTRUMENTS
#define BMC_FI_INSTRUMENTS

#include "dates.hpp"

#include <chrono>
#include <vector>

namespace adhoc {

enum class DayCountConvention {
    e_ACTUAL_360,
    e_ACTUAL_365_FIXED,
    e_ISDA_ACTUAL_ACTUAL,
    e_ISMA_30_360,
};

auto
GetCoverage(DayCountConvention const& dcc,
            std::chrono::year_month_day const& start_date,
            std::chrono::year_month_day const& end_date) -> double
{
    double coverage = 0.0;

    switch (dcc) {
        case DayCountConvention::e_ACTUAL_360:
            coverage = static_cast<double>(end_date - start_date) / 360.0;
            break;
        case DayCountConvention::e_ACTUAL_365_FIXED:
            coverage = static_cast<double>(end_date - start_date) / 365.0;
            break;
        default:
            coverage = static_cast<double>(end_date - start_date) / 360.0;
            break;
    }

    return coverage;
}

template<class Float = double>
struct RFRRate {
    std::chrono::year_month_day start_date;
    std::chrono::year_month_day end_date;
    DayCountConvention dcc = DayCountConvention::e_ACTUAL_360;
    Float spread;
    auto GetCoverage() const -> double
    {
        double coverage = 0.0;

        switch (dcc) {
            case DayCountConvention::e_ACTUAL_360:
                coverage = static_cast<double>(end_date - start_date) / 360.0;
                break;
            default:
                coverage = static_cast<double>(end_date - start_date) / 360.0;
                break;
        }

        return coverage;
    }
};

template<class Float = double>
struct CMSIndex {
    std::chrono::year_month_day fix_date;
    std::vector<std::chrono::year_month_day> oFixSch;
    std::vector<double> pFixCvg;
    std::vector<std::chrono::year_month_day> oFltSch;
    std::vector<double> pFltCvg;
};

// eventually these should belong to "formula" type
template<class Float = double>
struct CashflowSimple {
    std::chrono::year_month_day pay_date;
    Float amount;
};

template<class Float = double>
struct CashflowSimple2 {
    std::chrono::year_month_day pay_date;
    Float dcf;
};

template<class Float = double>
struct CashflowLinear {
    std::chrono::year_month_day pay_date;
    double amount;
    RFRRate<double> rfr;
};

template<class Float = double>
struct CashflowLinear2 {
    std::chrono::year_month_day pay_date;
    RFRRate<double> rfr;
};

template<class Float = double>
struct Swap {
    std::vector<CashflowSimple<Float> > fixed_leg;
    std::vector<CashflowLinear<Float> > float_leg;
};

template<class Float = double>
struct Swap2 {
    std::vector<CashflowSimple2<Float> > fixed_leg;
    std::vector<CashflowLinear2<Float> > float_leg;
};

template<class Float = double>
struct Option {
    std::chrono::year_month_day fix_date;
    std::optional<CMSIndex<Float> > cms;
    std::optional<RFRRate<Float> > rfr;
};

} // namespace adhoc

#endif // BMC_FI_INSTRUMENTS
