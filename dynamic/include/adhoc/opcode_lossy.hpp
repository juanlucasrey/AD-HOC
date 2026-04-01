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

#ifndef ADHOC_OPCODE_LOSSY_HPP
#define ADHOC_OPCODE_LOSSY_HPP

#include <cstdint>

namespace adhoc {

// Simplified opcodes for lossy tape - only additions and multiplications
enum class LossyOpCode : std::uint8_t {
    ADD_FROM_SEED,
    COPY,          // result = source
    COPY_MINUS,    // result = -source
    ADD,           // result += source
    SUB,           // result -= source
    MINUS_INPLACE, // result = -result
    MUL_INPLACE,   // result *= factor
    MUL_ADD,       // result += factor * source (multiply and accumulate)
    MUL_SET,       // result = factor * source (multiply and set)
};

} // namespace adhoc

#endif // ADHOC_OPCODE_LOSSY_HPP
