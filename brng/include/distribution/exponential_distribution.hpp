
/*
 * AD-HOC, Automatic Differentiation for High Order Calculations
 *
 * This file is part of the AD-HOC distribution
 * (https://github.com/juanlucasrey/AD-HOC).
 * Copyright (c) 2024 Juan Lucas Rey
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

#ifndef BRNG_DISTRIBUTION_EXPONENTIAL_DISTRIBUTION
#define BRNG_DISTRIBUTION_EXPONENTIAL_DISTRIBUTION

#include "centered_canonical_distribution.hpp"

#include <cmath>
#include <limits>
#include <type_traits>

namespace adhoc {

template <class RealType = double> class exponential_distribution {

    static_assert(std::is_floating_point_v<RealType>);

  public:
    using result_type = RealType;

  private:
    result_type _lambda;
    result_type _lambda_inv;

  public:
    exponential_distribution() : exponential_distribution(1) {}

    explicit exponential_distribution(result_type lambda_in)
        : _lambda(lambda_in), _lambda_inv(1.0 / lambda_in) {}

    void reset() {}

    template <class RNG> auto operator()(RNG &g) -> result_type {
        adhoc::centered_canonical_distribution<result_type> u;
        auto const val = u(g);

        if (std::copysign(1.0, val) > 0) {
            return -std::log(result_type{1} - val) * _lambda_inv;
        }

        return -std::log(-val) * _lambda_inv;
    }

    result_type lambda() const { return _lambda; }

    result_type min() const { return 0; }
    result_type max() const {
        return std::numeric_limits<result_type>::infinity();
    }

    auto operator==(const exponential_distribution &rhs) const -> bool {
        return (this->_lambda == rhs._lambda) &&
               (this->_lambda_inv == rhs._lambda_inv);
    }

    auto operator!=(const exponential_distribution &rhs) const -> bool {
        return !(this->operator==(rhs));
    }
};

} // namespace adhoc

#endif // BRNG_DISTRIBUTION_EXPONENTIAL_DISTRIBUTION
