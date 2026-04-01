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

#ifndef ADHOC_TAPE_DATA_HPP
#define ADHOC_TAPE_DATA_HPP

#include <cstdint>
#include <vector>

namespace adhoc {

enum class OpCode : std::uint8_t {
    ADD,
    SUB,
    MUL,
    ADD_C,
    SUB_C,
    MUL_C,
    NORM, // a.k.a. SQUARE
    ABS,
    INV,
    SQRT,
    POW_C,
    EXP,
    LOG,
    ERF,
    ERFC,
    COS,
    REG_INPUT,
    REG_OUTPUT,
    // IFT
    // IFT_START,
    // IFT_END,
    // REG_INPUT_IFT,
    // REG_OUTPUT_IFT,
};

struct TapeData {
    std::vector<OpCode> ops;
    std::vector<double> vals;
    std::vector<std::size_t> ids;
    size_t next_id{ 0 };
};

} // namespace adhoc

#endif // ADHOC_TAPE_DATA_HPP
