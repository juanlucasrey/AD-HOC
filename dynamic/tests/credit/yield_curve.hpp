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

// 1-piecewiselinear simple-comp
// interpolationMethod: LINEAR_SIMPLE
// ConstantSpotDFFunction, LinearSpotDFFunction, ConstantSpotDFFunction

// 2-smooth forward:
// interpolationMethod: SMOOTH_FORWARD
// Piecewise Cubic Function

// 3-step-function forward:
// interpolationMethod: STEP_FORWARD
// AnchoredStepForwardDFFunction, StepForwardDFFunction, ExponentialFunction

// 4-piecewiselinear continuous-comp
// interpolationMethod: LINEAR_CONTINUOUS
// ConstantCCSpotDFFunction, LinearCCSpotDFFunction, ConstantCCSpotDFFunction

#ifndef BMC_YIELD_CURVE
#define BMC_YIELD_CURVE

#include "dates.hpp"

#include <chrono>
#include <vector>

namespace adhoc {

template<class Float = double>
class yield_curve final {
  public:
    explicit yield_curve(std::chrono::year_month_day const& asofdate_in,
                         std::vector<Float> const& params_in,
                         std::vector<std::chrono::year_month_day> const& segments_in)
      : asofdate(asofdate_in)
    {

        dates.reserve(segments_in.size() + 1);
        dates.push_back(asofdate_in);
        dates.insert(end(dates), begin(segments_in), end(segments_in));

        params.reserve(params_in.size() + 1);
        params.push_back(1.);
        params.insert(end(params), begin(params_in), end(params_in));
    }

    auto getDiscountFactor(std::chrono::year_month_day const& date) const -> Float
    {
        auto const it = std::min(std::upper_bound(this->dates.begin(), this->dates.end(), date), this->dates.end() - 1);

        auto const d = std::distance(this->dates.begin(), it);
        auto const& df1 = this->params[static_cast<std::size_t>(d - 1)];
        auto const& df2 = this->params[static_cast<std::size_t>(d)];
        auto const df = df2 / df1;

        auto const proportion =
          static_cast<Float>(date - this->dates[static_cast<std::size_t>(d - 1)]) /
          static_cast<Float>(this->dates[static_cast<std::size_t>(d)] - this->dates[static_cast<std::size_t>(d - 1)]);

        using std::exp;
        using std::log;
        return df1 * exp(log(df) * proportion);
    }

  private:
    std::chrono::year_month_day asofdate;
    std::vector<Float> params;
    std::vector<std::chrono::year_month_day> dates;
};

template<class Float = double>
class yield_curve2 final {
  public:
    explicit yield_curve2(std::chrono::year_month_day const& asofdate_in,
                          std::vector<Float> const& params_in,
                          std::vector<std::chrono::year_month_day> const& segments_in)
      : asofdate(asofdate_in)
    {

        dates.reserve(segments_in.size() + 1);
        dates.push_back(asofdate_in);
        dates.insert(end(dates), begin(segments_in), end(segments_in));

        params.reserve(params_in.size() + 1);
        params.push_back(1.);
        params.insert(end(params), begin(params_in), end(params_in));
    }

    auto getDiscountFactor(std::chrono::year_month_day const& date) const -> Float
    {
        auto const it = std::min(std::upper_bound(this->dates.begin(), this->dates.end(), date), this->dates.end() - 1);

        auto const d = std::distance(this->dates.begin(), it);
        auto const& df1 = this->params[static_cast<std::size_t>(d - 1)];
        auto const& df2 = this->params[static_cast<std::size_t>(d)];
        auto const df = df2 * inv(df1);

        auto const proportion =
          static_cast<Float>(date - this->dates[static_cast<std::size_t>(d - 1)]) /
          static_cast<Float>(this->dates[static_cast<std::size_t>(d)] - this->dates[static_cast<std::size_t>(d - 1)]);

        using std::exp;
        using std::log;
        return df1 * exp(log(df) * proportion);
    }

  private:
    std::chrono::year_month_day asofdate;
    std::vector<Float> params;
    std::vector<std::chrono::year_month_day> dates;
};

} // namespace adhoc

#endif // BMC_YIELD_CURVE
