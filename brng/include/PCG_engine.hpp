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

#ifndef BRNG_PCG_ENGINE
#define BRNG_PCG_ENGINE

#include "tools/modular_arithmetic.hpp"

#include <cstdint>

namespace adhoc {

template <class UIntType = std::uint32_t> class PCG_engine final {

  public:
    using result_type = UIntType;

    static constexpr std::uint64_t default_seed1 = 0x853c49e6748fea9bULL;
    static constexpr std::uint64_t default_seed2 = 0xda3e39cb94b95bdbULL;

    PCG_engine() : PCG_engine(default_seed1, default_seed2) {}

    explicit PCG_engine(std::uint64_t initstate, std::uint64_t initseq)
        : state(initstate), inc(initseq) {}

    template <bool FwdDirection = true>
    inline auto operator()() -> result_type {
        if constexpr (FwdDirection) {
            std::uint64_t oldstate = this->state;
            this->state = oldstate * 6364136223846793005ULL + this->inc;
            UIntType xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
            UIntType rot = oldstate >> 59u;
            return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
        } else {
            uint64_t inc_inv = -this->inc;
            static constexpr std::uint64_t multiplier_inverse =
                modular_multiplicative_inverse_max_plus_one(
                    6364136223846793005ULL);
            this->state = multiplier_inverse * (this->state + inc_inv);
            UIntType xorshifted = ((this->state >> 18u) ^ this->state) >> 27u;
            UIntType rot = this->state >> 59u;
            return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
        }
    }

    auto operator==(const PCG_engine &rhs) const -> bool {
        return (this->state == rhs.state) && (this->inc == rhs.inc);
    }

  private:
    std::uint64_t state;
    std::uint64_t inc;
};

} // namespace adhoc

#endif // BRNG_PCG_ENGINE
