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

#ifndef BRNG_TOOLS_MODULAR_ARITHMETIC
#define BRNG_TOOLS_MODULAR_ARITHMETIC

#include "mask.hpp"

#include <type_traits>

namespace adhoc {

template <class UIntType>
constexpr auto gcd(UIntType a, UIntType b) -> UIntType {
    while (b != 0) {
        auto tmp = b;
        b = a % b;
        a = tmp;
    }
    return a;
};

template <class UIntType>
constexpr auto modular_multiplicative_inverse(UIntType n, UIntType b)
    -> UIntType {
    static_assert(std::is_unsigned_v<UIntType>);
    UIntType n_original = n;
    UIntType t1 = 0;
    bool t1_b = true; // true means positive, false negative
    UIntType t2 = 1;
    bool t2_b = true; // true means positive, false negative
    UIntType t3 = 0;
    bool t3_b = true; // true means positive, false negative
    UIntType r = 1;
    while (r != 0) {
        UIntType q = n / b;
        r = n % b;

        // t3 = t1 - q × t2
        if (t1_b && !t2_b) {
            t3_b = true;
            t3 = t1 + q * t2;
        } else if (!t1_b && t2_b) {
            t3_b = false;
            t3 = t1 + q * t2;
        } else if (t1_b && t2_b) {
            t3_b = t1 >= q * t2;
            if (t3_b) {
                t3 = t1 - q * t2;
            } else {
                t3 = q * t2 - t1;
            }
        } else {
            t3_b = q * t2 >= t1;
            if (t3_b) {
                t3 = q * t2 - t1;
            } else {
                t3 = t1 - q * t2;
            }
        }

        n = b;
        b = r;
        t1 = t2;
        t1_b = t2_b;
        t2 = t3;
        t2_b = t3_b;
    }

    if (t1_b) {
        UIntType q = t1 / n_original;
        r = t1 % n_original;
        return r;
    } else {
        UIntType q = t1 / n_original;
        r = t1 % n_original;
        return n_original - r;
    }
}

template <class UIntType, bool AltVersion = false>
constexpr auto modular_multiplicative_inverse_pow2(std::size_t pow2,
                                                   UIntType b) {
    static_assert(std::is_unsigned_v<UIntType>);
    UIntType power = mask<UIntType>(pow2);

    if constexpr (AltVersion) {
        UIntType result = UIntType(1);
        UIntType multiplier = b;
        while (power != UIntType(0)) {
            UIntType thismult = power & UIntType(1) ? multiplier : UIntType(1);
            result *= thismult;
            power >>= 1U;
            multiplier *= multiplier;
        }
        return result;
    } else {
        UIntType t1 = 0;
        bool t1_b = true; // true means positive, false negative
        UIntType t2 = 1;
        bool t2_b = true; // true means positive, false negative
        UIntType t3 = 0;
        bool t3_b = true; // true means positive, false negative

        // awkward first step to avoid overflow
        UIntType r = power - (b - 1);
        UIntType q = r / b + 1;
        r = r % b;

        // t3 = t1 - q × t2
        if (t1_b && !t2_b) {
            t3_b = true;
            t3 = t1 + q * t2;
        } else if (!t1_b && t2_b) {
            t3_b = false;
            t3 = t1 + q * t2;
        } else if (t1_b && t2_b) {
            t3_b = t1 >= q * t2;
            if (t3_b) {
                t3 = t1 - q * t2;
            } else {
                t3 = q * t2 - t1;
            }
        } else {
            t3_b = q * t2 >= t1;
            if (t3_b) {
                t3 = q * t2 - t1;
            } else {
                t3 = t1 - q * t2;
            }
        }

        UIntType n = b;
        b = r;
        t1 = t2;
        t1_b = t2_b;
        t2 = t3;
        t2_b = t3_b;

        while (r != 0) {
            UIntType q = n / b;
            r = n % b;

            // t3 = t1 - q × t2
            if (t1_b && !t2_b) {
                t3_b = true;
                t3 = t1 + q * t2;
            } else if (!t1_b && t2_b) {
                t3_b = false;
                t3 = t1 + q * t2;
            } else if (t1_b && t2_b) {
                t3_b = t1 >= q * t2;
                if (t3_b) {
                    t3 = t1 - q * t2;
                } else {
                    t3 = q * t2 - t1;
                }
            } else {
                t3_b = q * t2 >= t1;
                if (t3_b) {
                    t3 = q * t2 - t1;
                } else {
                    t3 = t1 - q * t2;
                }
            }

            n = b;
            b = r;
            t1 = t2;
            t1_b = t2_b;
            t2 = t3;
            t2_b = t3_b;
        }

        if (t1_b) {
            return t1;
        } else {
            UIntType q = power - (t1 - 1);
            return q;
        }
    }
}

} // namespace adhoc

#endif // BRNG_TOOLS_MODULAR_ARITHMETIC
