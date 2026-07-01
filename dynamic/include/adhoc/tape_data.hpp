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

#include "vector_enum.hpp"
#include "vector_idx.hpp"

#include <cstdint>
#include <vector>

namespace adhoc {

// Every time a new type of operation is added, the NumValuesOpcode must be updated to reflect the new number of
// operations. This is used for the vector_enum2 type to determine how many bits are needed to store the operation codes
// efficiently.
static constexpr std::size_t NumValuesOpcode = 18;
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

enum class EnumVectorType { Simple, BitCompression, Valuecompression };
enum class IdxVectorType { Simple, BitCompression };

template<EnumVectorType enumvectype, IdxVectorType idxvectype>
class TapeData {
  private:
    template<class T, std::size_t NumValues>
    using enumvector_t = std::conditional_t<enumvectype == EnumVectorType::Simple,
                                            std::vector<T>,
                                            std::conditional_t<enumvectype == EnumVectorType::BitCompression,
                                                               vector_enum<T, NumValues>,
                                                               vector_enum2<T, NumValues> > >;

    template<class T>
    using idxvector_t = std::conditional_t<idxvectype == IdxVectorType::Simple, std::vector<T>, vector_idx<T> >;

  public:
    static constexpr EnumVectorType tape_enumvector_t = enumvectype;
    static constexpr IdxVectorType tape_idxvector_t = idxvectype;

    enumvector_t<OpCode, NumValuesOpcode> ops;
    idxvector_t<std::size_t> ids;
    std::vector<double> vals;
    std::size_t next_id{ 0 };

    void record_unary(OpCode op, std::size_t arg_id, std::size_t res_id)
    {
        this->ops.push_back(op);
        this->ids.push_back(arg_id);
        this->ids.push_back(res_id);
    }

    void record_binary(OpCode op, std::size_t lhs_id, std::size_t rhs_id, std::size_t res_id)
    {
        this->ops.push_back(op);
        this->ids.push_back(lhs_id);
        this->ids.push_back(rhs_id);
        this->ids.push_back(res_id);
    }

    auto generate_id() -> std::size_t { return this->next_id++; }

    void record_value(double const value) { this->vals.push_back(value); }
};

} // namespace adhoc

#endif // ADHOC_TAPE_DATA_HPP
