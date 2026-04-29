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

#ifndef ADHOC_BACKPROPAGATOR1LOSSYCOMPRESSED_HPP
#define ADHOC_BACKPROPAGATOR1LOSSYCOMPRESSED_HPP

#include "passive_id.hpp"
#include "tape_data.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <vector>

namespace adhoc {

template<class Float, bool Vectorised = false>
class BackPropagatorLossyCompressed {
  private:
    enum class LossyOpCode : std::uint8_t {
        COPY,          // result = source
        COPY_MINUS,    // result = -source
        ADD,           // result += source
        SUB,           // result -= source
        MINUS_INPLACE, // result = -result
        MUL_INPLACE,   // result *= factor
        MUL_ADD,       // result += factor * source (multiply and accumulate)
        MUL_SET,       // result = factor * source (multiply and set)
    };

    std::size_t m_num_lanes{ 1 };

    std::vector<std::size_t> node_location_on_buffer;

    // lossy tape
    std::vector<std::uint8_t> on_which_buffer;
    std::vector<std::size_t> pos;
    std::vector<LossyOpCode> lossy_op;
    std::vector<double> values;

    std::size_t from_prev{ 0 };
    std::size_t to_prev{ 0 };

    struct buffer_t {
        std::vector<double> values;
        std::size_t size{ 0 };
        std::vector<std::size_t> free_positions;
    };

    std::vector<std::size_t> checkpoints{ 0 };
    std::vector<buffer_t> buffers{ buffer_t{} };

  public:
    explicit BackPropagatorLossyCompressed() = default;

    void set_checkpoint(std::size_t ops_size)
    {
        if (this->checkpoints.back() > ops_size) {
            // should not happen
            throw;
        }

        if (this->checkpoints.back() != ops_size) {
            this->checkpoints.push_back(ops_size);
            this->buffers.push_back(buffer_t{});
        }
    }
    void set_lanes(std::size_t num_lanes)
    {
        if constexpr (!Vectorised) {
            if (this->m_num_lanes != 1) {
                // This backpropagator is not designed for multiple lanes
                throw;
            }
        }
        this->m_num_lanes = num_lanes;
    }
    auto get_lanes() const -> std::size_t { return this->m_num_lanes; }
    void reserve_input(std::size_t count_registered)
    {
        this->buffers.back().values.reserve(this->buffers.back().values.size() +
                                            (this->m_num_lanes * count_registered));
    }
    void reserve_output(std::size_t count_registered)
    {
        this->buffers.back().values.reserve(this->buffers.back().values.size() +
                                            (this->m_num_lanes * count_registered));
    }
    void register_variable(std::size_t var_id)
    {
        this->node_location_on_buffer.resize(std::max(this->node_location_on_buffer.size(), var_id + 1),
                                             passive_id<std::size_t>);
        std::size_t& var_pos = this->node_location_on_buffer[var_id];
        if (var_pos == passive_id<std::size_t>) {
            auto it = std::upper_bound(checkpoints.cbegin(), checkpoints.cend(), var_id);
            auto& var_buffer = buffers[std::distance(checkpoints.cbegin(), it) - 1];
            if (var_buffer.free_positions.empty()) {
                var_pos = var_buffer.size;
                ++var_buffer.size;
                var_buffer.values.resize(var_buffer.size * this->m_num_lanes);
            }
            else {
                var_pos = var_buffer.free_positions.back();
                var_buffer.free_positions.pop_back();
                if constexpr (Vectorised) {
                    double* dest = &var_buffer.values[var_pos * this->m_num_lanes];
#pragma omp simd
                    for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                        dest[i] = 0.;
                    }
                }
                else {
                    var_buffer.values[var_pos] = 0.;
                }
            }
        }
    }

    void register_output_variable(std::size_t var_id, std::size_t /* ops_size */)
    {
        this->node_location_on_buffer.resize(std::max(this->node_location_on_buffer.size(), var_id + 1),
                                             passive_id<std::size_t>);
        std::size_t& var_pos = this->node_location_on_buffer[var_id];
        if (var_pos == passive_id<std::size_t>) {
            auto it = std::upper_bound(checkpoints.cbegin(), checkpoints.cend(), var_id);
            auto& var_buffer = buffers[std::distance(checkpoints.cbegin(), it) - 1];
            if (var_buffer.free_positions.empty()) {
                var_pos = var_buffer.size;
                ++var_buffer.size;
                var_buffer.values.resize(var_buffer.size * this->m_num_lanes);
            }
            else {
                var_pos = var_buffer.free_positions.back();
                var_buffer.free_positions.pop_back();
                if constexpr (Vectorised) {
                    double* dest = &var_buffer.values[var_pos * this->m_num_lanes];
#pragma omp simd
                    for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                        dest[i] = 0.;
                    }
                }
                else {
                    var_buffer.values[var_pos] = 0.;
                }
            }
        }
    }

    void set_derivative(std::size_t var_id, double deriv, std::size_t /* ops_size */, std::size_t lane = 0)
    {
        if (lane < this->m_num_lanes) {
            std::size_t const var_pos = this->node_location_on_buffer[var_id];
            if (var_pos == passive_id<std::size_t>) {
                // this derivative is not on buffer.
                // this is probably nor an input nor an output.
                throw;
            }

            auto it = std::upper_bound(checkpoints.cbegin(), checkpoints.cend(), var_id);
            auto buffed_id = static_cast<std::uint8_t>(std::distance(checkpoints.cbegin(), it) - 1);
            buffers[buffed_id].values[(var_pos * this->m_num_lanes) + lane] = deriv;
        }
        else {
            throw;
        }
    }

    void set_derivative(std::size_t /* var_id1 */, std::size_t /* var_id2 */, double /* deriv */)
    {
        // not supported for lossy tape (first order only)
        throw;
    }

    auto get_derivative(std::size_t var_id, std::size_t lane) const -> double
    {
        if (lane < this->m_num_lanes) {
            std::size_t const var_pos = this->node_location_on_buffer[var_id];
            if (var_pos == passive_id<std::size_t>) {
                // this derivative is not on buffer.
                // this is probably nor an input nor an output.
                throw;
            }

            auto it = std::upper_bound(checkpoints.cbegin(), checkpoints.cend(), var_id);
            return buffers[std::distance(checkpoints.cbegin(), it) - 1].values[(var_pos * this->m_num_lanes) + lane];
        }
        else {
            throw;
        }
        return 0.;
    }

    auto get_derivative(std::size_t /* var_id1 */, std::size_t /* var_id2 */, std::size_t /* lane */) const -> double
    {
        // order is too high for lossy tape
        return 0.;
    }

    void clear() {}

    void zero_adjoints()
    {
        for (auto& b : buffers) {
            std::fill(b.values.begin(), b.values.end(), 0.0);
        }
    }

    auto size_of(bool capacity = false) const -> std::size_t
    {
        std::size_t size = 0;
        size += 3 * sizeof(std::size_t); // m_num_lanes, from_prev, to_prev
        size += sizeof(std::size_t) * (capacity ? node_location_on_buffer.capacity() : node_location_on_buffer.size());
        size += sizeof(std::uint8_t) * (capacity ? on_which_buffer.capacity() : on_which_buffer.size());
        size += sizeof(std::size_t) * (capacity ? pos.capacity() : pos.size());
        size += sizeof(LossyOpCode) * (capacity ? lossy_op.capacity() : lossy_op.size());
        size += sizeof(double) * (capacity ? values.capacity() : values.size());
        size += sizeof(std::size_t) * (capacity ? checkpoints.capacity() : checkpoints.size());
        for (const auto& buffer : buffers) {
            size += sizeof(double) * (capacity ? buffer.values.capacity() : buffer.values.size());
            size += sizeof(std::size_t) * (capacity ? buffer.free_positions.capacity() : buffer.free_positions.size());
        }
        return size;
    }

    template<bool Reset = false, bool ResetInPlace = false, bool Log = false>
    void backpropagate_to(std::size_t to, TapeData& data);
};

template<class Float, bool Vectorised>
template<bool Reset, bool ResetInPlace, bool Log>
void
BackPropagatorLossyCompressed<Float, Vectorised>::backpropagate_to(std::size_t to, TapeData& data)
{
    std::size_t from = data.next_id;
    if (from == checkpoints.back()) {
        checkpoints.pop_back();
        buffers.pop_back();
    }

    const auto& ops = data.ops;
    const auto& vals = data.vals;
    const auto& ids = data.ids;

    std::size_t val_idx = vals.size();
    std::size_t id_idx = ids.size();
    if (this->lossy_op.empty() || from_prev != from || to_prev != to) {
        from_prev = from;
        to_prev = to;

        this->on_which_buffer.clear();
        this->pos.clear();
        this->lossy_op.clear();
        this->values.clear();

        this->node_location_on_buffer.resize(ops.size(), passive_id<std::size_t>);

        auto& buffer_free_positions = buffers.back().free_positions;

        std::vector<std::size_t> number_dependents(from - to);

        val_idx = vals.size();
        id_idx = ids.size();
        for (std::size_t op_idx = from; op_idx-- > to;) {
            OpCode const& op = ops[op_idx];
            bool const use_this_op =
              this->node_location_on_buffer[op_idx] != passive_id<std::size_t> || (number_dependents[op_idx - to] > 0);

            switch (op) {
                case OpCode::REG_INPUT:
                case OpCode::REG_OUTPUT:
                case OpCode::ADD:
                case OpCode::SUB:
                case OpCode::ADD_C:
                case OpCode::SUB_C: {
                    break;
                }
                case OpCode::MUL_C:
                case OpCode::NORM:
                case OpCode::INV:
                case OpCode::ABS:
                case OpCode::EXP:
                case OpCode::LOG:
                case OpCode::ERF:
                case OpCode::ERFC: {
                    val_idx -= 1;
                    break;
                }
                case OpCode::MUL:
                case OpCode::COS:
                case OpCode::SQRT: {
                    val_idx -= 2;
                    break;
                }
                case OpCode::POW_C: {
                    val_idx -= 3;
                    break;
                }
            }

            switch (op) {
                case OpCode::REG_INPUT: {
                    id_idx -= 1;
                    break;
                }
                case OpCode::REG_OUTPUT:
                case OpCode::ADD_C:
                case OpCode::SUB_C:
                case OpCode::MUL_C:
                case OpCode::NORM:
                case OpCode::INV:
                case OpCode::ABS:
                case OpCode::EXP:
                case OpCode::LOG:
                case OpCode::ERF:
                case OpCode::ERFC:
                case OpCode::COS:
                case OpCode::SQRT:
                case OpCode::POW_C: {
                    id_idx -= 2;
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        if (arg_id >= to) {
                            number_dependents[arg_id - to] += 1;
                        }
                    }
                    break;
                }
                case OpCode::ADD:
                case OpCode::SUB:
                case OpCode::MUL: {
                    id_idx -= 3;
                    if (use_this_op) {
                        std::size_t const lhs_id = ids[id_idx];
                        std::size_t const rhs_id = ids[id_idx + 1];
                        if (lhs_id >= to) {
                            number_dependents[lhs_id - to] += 1;
                        }
                        if (rhs_id >= to) {
                            number_dependents[rhs_id - to] += 1;
                        }
                    }
                    break;
                }
            }
        }
        std::size_t id_idx_start = id_idx;
        std::size_t val_idx_start = val_idx;

        // loop in foward order to detect multiplication chains and bivariate operators with same agument
        buffer_t buffer_multipliers;

        std::vector<std::size_t> pos_multiplier;
        std::vector<bool> op_multiply;
        std::vector<std::size_t> multiplier_loc_from;
        std::vector<std::size_t> multiplier_loc_to;
        std::vector<bool> bivariate_consolidate;
        std::vector<bool> univariate_consolidate;
        std::vector<std::size_t> multiplier_origin(from - to, passive_id<std::size_t>);
        std::vector<std::size_t> multiplier_origin_bi((from - to) * 2, passive_id<std::size_t>);
        std::vector<bool> is_output(from - to, false);

        id_idx = id_idx_start;
        for (std::size_t op_idx = to; op_idx < from; ++op_idx) {
            OpCode const& op = ops[op_idx];
            bool const use_this_op =
              this->node_location_on_buffer[op_idx] != passive_id<std::size_t> || (number_dependents[op_idx - to] > 0);

            switch (op) {
                case OpCode::REG_INPUT: {
                    id_idx += 1;
                    break;
                }
                case OpCode::REG_OUTPUT: {
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];

                        is_output[res_id - to] = true;
                        bool arg_is_induced_path = (arg_id >= to) && (number_dependents[arg_id - to] == 1) &&
                                                   (multiplier_origin[arg_id - to] != passive_id<std::size_t>);
                        if (arg_is_induced_path) {
                            auto& mul_origin = multiplier_origin[arg_id - to];
                            multiplier_loc_to[mul_origin] = res_id;
                            multiplier_origin[res_id - to] = mul_origin;
                            mul_origin = passive_id<std::size_t>;
                        }
                        else {
                            op_multiply.push_back(false);

                            std::size_t mult_loc_res = 0;
                            if (buffer_multipliers.free_positions.empty()) {
                                mult_loc_res = buffer_multipliers.size;
                                ++buffer_multipliers.size;
                                multiplier_loc_from.resize(buffer_multipliers.size);
                                multiplier_loc_to.resize(buffer_multipliers.size);
                            }
                            else {
                                mult_loc_res = buffer_multipliers.free_positions.back();
                                buffer_multipliers.free_positions.pop_back();
                            }

                            pos_multiplier.push_back(mult_loc_res);
                            multiplier_loc_from[mult_loc_res] = arg_id;
                            multiplier_loc_to[mult_loc_res] = res_id;
                            multiplier_origin[res_id - to] = mult_loc_res;
                        }
                    }
                    id_idx += 2;
                    break;
                }
                // case OpCode::SUB: {
                //     if (use_this_op) {
                //         // todo
                //         throw;
                //         std::size_t const lhs_id = ids[id_idx];
                //         std::size_t const rhs_id = ids[id_idx + 1];
                //         std::size_t const res_id = ids[id_idx + 2];
                //     }
                //     id_idx += 3;
                //     break;
                // }
                case OpCode::SUB:
                case OpCode::ADD:
                case OpCode::MUL: {
                    if (use_this_op) {
                        std::size_t const lhs_id = ids[id_idx];
                        std::size_t const rhs_id = ids[id_idx + 1];
                        std::size_t const res_id = ids[id_idx + 2];

                        bool lhs_is_induced_path = (lhs_id >= to) && (number_dependents[lhs_id - to] == 1) &&
                                                   (multiplier_origin[lhs_id - to] != passive_id<std::size_t>);
                        bool rhs_is_induced_path = (rhs_id >= to) && (number_dependents[rhs_id - to] == 1) &&
                                                   (multiplier_origin[rhs_id - to] != passive_id<std::size_t>);

                        std::size_t multiplier_loc_lhs = 0;
                        if (lhs_is_induced_path) {
                            op_multiply.push_back(true);
                            auto& coming_from = multiplier_origin[lhs_id - to];
                            pos_multiplier.push_back(coming_from);
                            multiplier_loc_to[coming_from] = res_id;
                            multiplier_loc_lhs = coming_from;
                            coming_from = passive_id<std::size_t>;
                        }
                        else {
                            op_multiply.push_back(false);

                            std::size_t mult_loc_res = 0;
                            if (buffer_multipliers.free_positions.empty()) {
                                mult_loc_res = buffer_multipliers.size;
                                ++buffer_multipliers.size;
                                multiplier_loc_from.resize(buffer_multipliers.size);
                                multiplier_loc_to.resize(buffer_multipliers.size);
                            }
                            else {
                                mult_loc_res = buffer_multipliers.free_positions.back();
                                buffer_multipliers.free_positions.pop_back();
                            }

                            pos_multiplier.push_back(mult_loc_res);
                            multiplier_loc_from[mult_loc_res] = lhs_id;
                            multiplier_loc_to[mult_loc_res] = res_id;

                            multiplier_loc_lhs = mult_loc_res;
                        }

                        std::size_t multiplier_loc_rhs = 0;
                        if (rhs_is_induced_path) {
                            op_multiply.push_back(true);
                            auto& coming_from = multiplier_origin[rhs_id - to];
                            pos_multiplier.push_back(coming_from);
                            multiplier_loc_to[coming_from] = res_id;
                            multiplier_loc_rhs = coming_from;
                            coming_from = passive_id<std::size_t>;
                        }
                        else {
                            op_multiply.push_back(false);

                            std::size_t mult_loc_res = 0;
                            if (buffer_multipliers.free_positions.empty()) {
                                mult_loc_res = buffer_multipliers.size;
                                ++buffer_multipliers.size;
                                multiplier_loc_from.resize(buffer_multipliers.size);
                                multiplier_loc_to.resize(buffer_multipliers.size);
                            }
                            else {
                                mult_loc_res = buffer_multipliers.free_positions.back();
                                buffer_multipliers.free_positions.pop_back();
                            }

                            pos_multiplier.push_back(mult_loc_res);
                            multiplier_loc_from[mult_loc_res] = rhs_id;
                            multiplier_loc_to[mult_loc_res] = res_id;

                            multiplier_loc_rhs = mult_loc_res;
                        }

                        bool has_induced_path = lhs_is_induced_path || rhs_is_induced_path;
                        // there is potential for a bivariate operator on the same argument, we need to check if
                        // this is the case and update the multiplier if so
                        bool bivariate_consolidate_this =
                          has_induced_path &&
                          (multiplier_loc_from[multiplier_loc_lhs] == multiplier_loc_from[multiplier_loc_rhs]);
                        bivariate_consolidate.push_back(bivariate_consolidate_this);

                        if (bivariate_consolidate_this) {
                            std::size_t origin_id = multiplier_loc_from[multiplier_loc_lhs];
                            pos_multiplier.push_back(multiplier_loc_lhs);
                            pos_multiplier.push_back(multiplier_loc_rhs);
                            buffer_multipliers.free_positions.push_back(multiplier_loc_rhs);
                            multiplier_loc_from[multiplier_loc_rhs] = passive_id<std::size_t>;
                            multiplier_loc_to[multiplier_loc_rhs] = passive_id<std::size_t>;
                            multiplier_origin[res_id - to] = multiplier_loc_lhs;

                            auto& number_dependents_to_update = number_dependents[origin_id - to];
                            --number_dependents_to_update;

                            bool univariate_consolidate_this =
                              (number_dependents_to_update == 1) &&
                              (multiplier_origin[origin_id - to] != passive_id<std::size_t>);
                            univariate_consolidate.push_back(univariate_consolidate_this);
                            if (univariate_consolidate_this) {
                                // this node now has only one dependent, we can reintroduce a
                                // multiplication chain
                                auto& coming_from = multiplier_origin[origin_id - to];
                                pos_multiplier.push_back(coming_from);
                                multiplier_loc_to[coming_from] = res_id;

                                buffer_multipliers.free_positions.push_back(multiplier_loc_lhs);
                                multiplier_loc_from[multiplier_loc_lhs] = passive_id<std::size_t>;
                                multiplier_loc_to[multiplier_loc_lhs] = passive_id<std::size_t>;

                                multiplier_origin[res_id - to] = coming_from;
                                coming_from = passive_id<std::size_t>;
                            }
                        }
                        else {
                            multiplier_origin_bi[(res_id - to) * 2] = multiplier_loc_lhs;
                            multiplier_origin_bi[((res_id - to) * 2) + 1] = multiplier_loc_rhs;
                        }
                    }
                    id_idx += 3;
                    break;
                }
                // case OpCode::ADD_C: {
                //     if (use_this_op) {
                //         // todo
                //         throw;
                //         std::size_t const arg_id = ids[id_idx];
                //         std::size_t const res_id = ids[id_idx + 1];
                //     }
                //     id_idx += 2;
                //     break;
                // }
                // case OpCode::SUB_C: {
                //     if (use_this_op) {
                //         // todo
                //         throw;
                //         std::size_t const arg_id = ids[id_idx];
                //         std::size_t const res_id = ids[id_idx + 1];
                //     }
                //     id_idx += 2;
                //     break;
                // }
                case OpCode::SUB_C:
                case OpCode::ADD_C:
                case OpCode::MUL_C:
                case OpCode::NORM:
                case OpCode::INV:
                case OpCode::ABS:
                case OpCode::EXP:
                case OpCode::LOG:
                case OpCode::ERF:
                case OpCode::ERFC:
                case OpCode::COS:
                case OpCode::SQRT:
                case OpCode::POW_C: {
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];

                        bool arg_is_induced_path = (arg_id >= to) && (number_dependents[arg_id - to] == 1) &&
                                                   (multiplier_origin[arg_id - to] != passive_id<std::size_t>);

                        if (arg_is_induced_path) {
                            op_multiply.push_back(true);
                            auto& coming_from = multiplier_origin[arg_id - to];
                            pos_multiplier.push_back(coming_from);
                            multiplier_loc_to[coming_from] = res_id;
                            multiplier_origin[res_id - to] = coming_from;
                            coming_from = passive_id<std::size_t>;
                        }
                        else {
                            auto& mult_loc_res = multiplier_origin[res_id - to];
                            if (buffer_multipliers.free_positions.empty()) {
                                mult_loc_res = buffer_multipliers.size;
                                ++buffer_multipliers.size;
                                multiplier_loc_from.resize(buffer_multipliers.size);
                                multiplier_loc_to.resize(buffer_multipliers.size);
                            }
                            else {
                                mult_loc_res = buffer_multipliers.free_positions.back();
                                buffer_multipliers.free_positions.pop_back();
                            }

                            op_multiply.push_back(false);
                            pos_multiplier.push_back(mult_loc_res);
                            multiplier_loc_from[mult_loc_res] = arg_id;
                            multiplier_loc_to[mult_loc_res] = res_id;
                        }
                    }
                    id_idx += 2;
                    break;
                }
            }
        }

        buffer_multipliers.values.resize(buffer_multipliers.size * this->m_num_lanes, 1.);

        id_idx = id_idx_start;
        val_idx = val_idx_start;

        std::size_t mult_op_idx = 0;
        std::size_t mult_pos_idx = 0;
        std::size_t mult_bc_idx = 0;
        std::size_t mult_uc_idx = 0;
        for (std::size_t op_idx = to; op_idx < from; ++op_idx) {
            OpCode const& op = ops[op_idx];
            bool const use_this_op =
              this->node_location_on_buffer[op_idx] != passive_id<std::size_t> || (number_dependents[op_idx - to] > 0);

            switch (op) {
                case OpCode::REG_INPUT: {
                    id_idx += 1;
                    break;
                }
                case OpCode::REG_OUTPUT: {
                    id_idx += 2;
                    break;
                }
                case OpCode::ADD: {
                    if (use_this_op) {
                        double const lhs_val = 1.0;
                        double const rhs_val = 1.0;

                        // lhs multiplier
                        {
                            double const multiplier = rhs_val;
                            bool const op_m = op_multiply[mult_op_idx++];
                            std::size_t const pos = pos_multiplier[mult_pos_idx++];
                            if (op_m) {
                                buffer_multipliers.values[pos] *= multiplier;
                            }
                            else {
                                buffer_multipliers.values[pos] = multiplier;
                            }
                        }

                        // rhs multiplier
                        {
                            double const multiplier = lhs_val;
                            bool const op_m = op_multiply[mult_op_idx++];
                            std::size_t const pos = pos_multiplier[mult_pos_idx++];
                            if (op_m) {
                                buffer_multipliers.values[pos] *= multiplier;
                            }
                            else {
                                buffer_multipliers.values[pos] = multiplier;
                            }
                        }

                        if (bivariate_consolidate[mult_bc_idx++]) {
                            std::size_t const pos_lhs = pos_multiplier[mult_pos_idx++];
                            std::size_t const pos_rhs = pos_multiplier[mult_pos_idx++];
                            buffer_multipliers.values[pos_lhs] += buffer_multipliers.values[pos_rhs];
                            if (univariate_consolidate[mult_uc_idx++]) {
                                std::size_t const pos = pos_multiplier[mult_pos_idx++];
                                buffer_multipliers.values[pos] *= buffer_multipliers.values[pos_lhs];
                            }
                        }
                    }
                    id_idx += 3;
                    break;
                }
                case OpCode::SUB: {
                    if (use_this_op) {
                        double const lhs_val = 1.0;
                        double const rhs_val = -1.0;

                        // lhs multiplier
                        {
                            double const multiplier = rhs_val;
                            bool const op_m = op_multiply[mult_op_idx++];
                            std::size_t const pos = pos_multiplier[mult_pos_idx++];
                            if (op_m) {
                                buffer_multipliers.values[pos] *= multiplier;
                            }
                            else {
                                buffer_multipliers.values[pos] = multiplier;
                            }
                        }

                        // rhs multiplier
                        {
                            double const multiplier = lhs_val;
                            bool const op_m = op_multiply[mult_op_idx++];
                            std::size_t const pos = pos_multiplier[mult_pos_idx++];
                            if (op_m) {
                                buffer_multipliers.values[pos] *= multiplier;
                            }
                            else {
                                buffer_multipliers.values[pos] = multiplier;
                            }
                        }

                        if (bivariate_consolidate[mult_bc_idx++]) {
                            std::size_t const pos_lhs = pos_multiplier[mult_pos_idx++];
                            std::size_t const pos_rhs = pos_multiplier[mult_pos_idx++];
                            buffer_multipliers.values[pos_lhs] += buffer_multipliers.values[pos_rhs];
                            if (univariate_consolidate[mult_uc_idx++]) {
                                std::size_t const pos = pos_multiplier[mult_pos_idx++];
                                buffer_multipliers.values[pos] *= buffer_multipliers.values[pos_lhs];
                            }
                        }
                    }
                    id_idx += 3;
                    break;
                }
                case OpCode::MUL: {
                    if (use_this_op) {

                        double const lhs_val = vals[val_idx];
                        double const rhs_val = vals[val_idx + 1];

                        // lhs multiplier
                        {
                            double const multiplier = rhs_val;
                            bool const op_m = op_multiply[mult_op_idx++];
                            std::size_t const pos = pos_multiplier[mult_pos_idx++];
                            if (op_m) {
                                buffer_multipliers.values[pos] *= multiplier;
                            }
                            else {
                                buffer_multipliers.values[pos] = multiplier;
                            }
                        }

                        // rhs multiplier
                        {
                            double const multiplier = lhs_val;
                            bool const op_m = op_multiply[mult_op_idx++];
                            std::size_t const pos = pos_multiplier[mult_pos_idx++];
                            if (op_m) {
                                buffer_multipliers.values[pos] *= multiplier;
                            }
                            else {
                                buffer_multipliers.values[pos] = multiplier;
                            }
                        }

                        if (bivariate_consolidate[mult_bc_idx++]) {
                            std::size_t const pos_lhs = pos_multiplier[mult_pos_idx++];
                            std::size_t const pos_rhs = pos_multiplier[mult_pos_idx++];
                            buffer_multipliers.values[pos_lhs] += buffer_multipliers.values[pos_rhs];
                            if (univariate_consolidate[mult_uc_idx++]) {
                                std::size_t const pos = pos_multiplier[mult_pos_idx++];
                                buffer_multipliers.values[pos] *= buffer_multipliers.values[pos_lhs];
                            }
                        }
                    }
                    val_idx += 2;
                    id_idx += 3;
                    break;
                }
                case OpCode::ADD_C: {
                    if (use_this_op) {
                        double const multiplier = 1.;
                        bool const op_m = op_multiply[mult_op_idx++];
                        std::size_t const pos = pos_multiplier[mult_pos_idx++];
                        if (op_m) {
                            buffer_multipliers.values[pos] *= multiplier;
                        }
                        else {
                            buffer_multipliers.values[pos] = multiplier;
                        }
                    }
                    id_idx += 2;
                    break;
                }
                case OpCode::SUB_C: {
                    if (use_this_op) {
                        double const multiplier = -1.;
                        bool const op_m = op_multiply[mult_op_idx++];
                        std::size_t const pos = pos_multiplier[mult_pos_idx++];
                        if (op_m) {
                            buffer_multipliers.values[pos] *= multiplier;
                        }
                        else {
                            buffer_multipliers.values[pos] = multiplier;
                        }
                    }
                    id_idx += 2;
                    break;
                }
                case OpCode::MUL_C: {
                    if (use_this_op) {
                        double const multiplier = vals[val_idx];
                        bool const op_m = op_multiply[mult_op_idx++];
                        std::size_t const pos = pos_multiplier[mult_pos_idx++];
                        if (op_m) {
                            buffer_multipliers.values[pos] *= multiplier;
                        }
                        else {
                            buffer_multipliers.values[pos] = multiplier;
                        }
                    }
                    val_idx += 1;
                    id_idx += 2;
                    break;
                }
                case OpCode::NORM: {
                    if (use_this_op) {
                        double const multiplier = 2.0 * vals[val_idx];
                        bool const op_m = op_multiply[mult_op_idx++];
                        std::size_t const pos = pos_multiplier[mult_pos_idx++];
                        if (op_m) {
                            buffer_multipliers.values[pos] *= multiplier;
                        }
                        else {
                            buffer_multipliers.values[pos] = multiplier;
                        }
                    }
                    val_idx += 1;
                    id_idx += 2;
                    break;
                }
                case OpCode::INV: {
                    if (use_this_op) {
                        double const multiplier = -vals[val_idx] * vals[val_idx];
                        bool const op_m = op_multiply[mult_op_idx++];
                        std::size_t const pos = pos_multiplier[mult_pos_idx++];
                        if (op_m) {
                            buffer_multipliers.values[pos] *= multiplier;
                        }
                        else {
                            buffer_multipliers.values[pos] = multiplier;
                        }
                    }
                    val_idx += 1;
                    id_idx += 2;
                    break;
                }
                case OpCode::EXP: {
                    if (use_this_op) {
                        double const multiplier = vals[val_idx];
                        bool const op_m = op_multiply[mult_op_idx++];
                        std::size_t const pos = pos_multiplier[mult_pos_idx++];
                        if (op_m) {
                            buffer_multipliers.values[pos] *= multiplier;
                        }
                        else {
                            buffer_multipliers.values[pos] = multiplier;
                        }
                    }
                    val_idx += 1;
                    id_idx += 2;
                    break;
                }
                case OpCode::LOG: {
                    if (use_this_op) {
                        double const multiplier = 1.0 / vals[val_idx];
                        bool const op_m = op_multiply[mult_op_idx++];
                        std::size_t const pos = pos_multiplier[mult_pos_idx++];
                        if (op_m) {
                            buffer_multipliers.values[pos] *= multiplier;
                        }
                        else {
                            buffer_multipliers.values[pos] = multiplier;
                        }
                    }
                    val_idx += 1;
                    id_idx += 2;
                    break;
                }
                case OpCode::ERF: {
                    if (use_this_op) {
                        constexpr double two_over_root_pi = 2. * std::numbers::inv_sqrtpi_v<double>;
                        double const multiplier = std::exp(-vals[val_idx] * vals[val_idx]) * two_over_root_pi;
                        bool const op_m = op_multiply[mult_op_idx++];
                        std::size_t const pos = pos_multiplier[mult_pos_idx++];
                        if (op_m) {
                            buffer_multipliers.values[pos] *= multiplier;
                        }
                        else {
                            buffer_multipliers.values[pos] = multiplier;
                        }
                    }
                    val_idx += 1;
                    id_idx += 2;
                    break;
                }
                case OpCode::ERFC: {
                    if (use_this_op) {
                        constexpr double minus_two_over_root_pi = -2. * std::numbers::inv_sqrtpi_v<double>;
                        double const multiplier = std::exp(-vals[val_idx] * vals[val_idx]) * minus_two_over_root_pi;
                        bool const op_m = op_multiply[mult_op_idx++];
                        std::size_t const pos = pos_multiplier[mult_pos_idx++];
                        if (op_m) {
                            buffer_multipliers.values[pos] *= multiplier;
                        }
                        else {
                            buffer_multipliers.values[pos] = multiplier;
                        }
                    }
                    val_idx += 1;
                    id_idx += 2;
                    break;
                }
                case OpCode::COS: {
                    if (use_this_op) {
                        double const multiplier = -std::sin(vals[val_idx]);
                        bool const op_m = op_multiply[mult_op_idx++];
                        std::size_t const pos = pos_multiplier[mult_pos_idx++];
                        if (op_m) {
                            buffer_multipliers.values[pos] *= multiplier;
                        }
                        else {
                            buffer_multipliers.values[pos] = multiplier;
                        }
                    }
                    val_idx += 2;
                    id_idx += 2;
                    break;
                }
                default: {
                    throw;
                    break;
                }
            }
        }

        auto const& checkpoints_c = this->checkpoints;
        auto get_loc = [checkpoints_c](std::size_t id) -> std::tuple<bool, std::uint8_t> {
            auto it = std::upper_bound(checkpoints_c.begin(), checkpoints_c.end(), id);
            auto buffer_id = static_cast<std::uint8_t>(std::distance(checkpoints_c.cbegin(), it) - 1);
            return { it == checkpoints_c.end(), buffer_id };
        };

        auto get_loc_buffer = [this](std::size_t& arg_pos, std::uint8_t buffer_id) {
            if (arg_pos == passive_id<std::size_t>) {
                auto& arg_buffer = this->buffers[buffer_id];
                if (arg_buffer.free_positions.empty()) {
                    arg_pos = arg_buffer.size;
                    ++arg_buffer.size;
                }
                else {
                    arg_pos = arg_buffer.free_positions.back();
                    arg_buffer.free_positions.pop_back();
                }
            }
        };

        auto update_loc = [this, get_loc_buffer](std::size_t res_pos, std::size_t& arg_pos, std::uint8_t buffer_id) {
            get_loc_buffer(arg_pos, buffer_id);
            this->pos.push_back(res_pos);
            this->pos.push_back(arg_pos);
            this->on_which_buffer.push_back(buffer_id);
        };

        auto copy_mul =
          [this, update_loc](std::size_t res_pos, std::size_t& arg_pos, std::uint8_t buffer_id, double multiplier) {
              bool arg_is_new = (arg_pos == passive_id<std::size_t>);
              update_loc(res_pos, arg_pos, buffer_id);

              if (arg_is_new) {
                  // this is a new value, we NEED to override
                  this->lossy_op.push_back(LossyOpCode::MUL_SET);
              }
              else {
                  this->lossy_op.push_back(LossyOpCode::MUL_ADD);
              }
              this->values.push_back(multiplier);
          };

        for (std::size_t op_idx = from; op_idx-- > to;) {
            std::size_t const this_multiplier_origin = multiplier_origin[op_idx - to];
            if (this_multiplier_origin != passive_id<std::size_t>) {
                std::size_t const arg_id = multiplier_loc_from[this_multiplier_origin];
                std::size_t const res_id = multiplier_loc_to[this_multiplier_origin];

                bool const this_is_output = is_output[res_id - to];
                std::size_t const res_pos = this->node_location_on_buffer[res_id];
                if (!this_is_output) {
                    node_location_on_buffer[res_id] = passive_id<std::size_t>;
                }
                std::size_t& arg_pos = node_location_on_buffer[arg_id];

                double const multiplier = buffer_multipliers.values[this_multiplier_origin];

                auto const arg_pos_data = get_loc(arg_id);
                bool arg_is_new = (arg_pos == passive_id<std::size_t>);
                bool arg_inplace = arg_is_new && std::get<0>(arg_pos_data);

                if (arg_inplace && !this_is_output) {
                    // res id should now be arg id, avoiding a copy and a potential buffer increase
                    arg_pos = res_pos;
                    this->pos.push_back(res_pos);
                    this->values.push_back(multiplier);
                    this->lossy_op.push_back(LossyOpCode::MUL_INPLACE);
                }
                else {
                    // std::uint8_t const buffer_id = std::get<1>(arg_pos_data);

                    // if (arg_is_new) {
                    //     auto& arg_buffer = buffers[buffer_id];
                    //     if (arg_buffer.free_positions.empty()) {
                    //         arg_pos = arg_buffer.size;
                    //         ++arg_buffer.size;
                    //     }
                    //     else {
                    //         arg_pos = arg_buffer.free_positions.back();
                    //         arg_buffer.free_positions.pop_back();
                    //     }
                    //     lossy_op.push_back(LossyOpCode::MUL_SET);
                    // }
                    // else {
                    //     lossy_op.push_back(LossyOpCode::MUL_ADD);
                    // }

                    // pos.push_back(res_pos);
                    // pos.push_back(arg_pos);
                    // on_which_buffer.push_back(buffer_id);
                    // values.push_back(multiplier);

                    copy_mul(res_pos, arg_pos, std::get<1>(arg_pos_data), multiplier);

                    if (!this_is_output) {
                        buffer_free_positions.push_back(res_pos);
                    }
                }
            }
            else {
                std::size_t const this_multiplier_origin_bi1 = multiplier_origin_bi[(op_idx - to) * 2];
                std::size_t const this_multiplier_origin_bi2 = multiplier_origin_bi[((op_idx - to) * 2) + 1];
                if (this_multiplier_origin_bi1 != passive_id<std::size_t>) {
                    double const mul_1 = buffer_multipliers.values[this_multiplier_origin_bi1];
                    double const mul_2 = buffer_multipliers.values[this_multiplier_origin_bi2];

                    std::size_t const lhs_id = multiplier_loc_from[this_multiplier_origin_bi1];
                    std::size_t const rhs_id = multiplier_loc_from[this_multiplier_origin_bi2];
                    std::size_t const res_id = multiplier_loc_to[this_multiplier_origin_bi1];

                    auto const res_pos = node_location_on_buffer[res_id];
                    node_location_on_buffer[res_id] = passive_id<std::size_t>;
                    std::size_t& lhs_pos = node_location_on_buffer[lhs_id];
                    std::size_t& rhs_pos = node_location_on_buffer[rhs_id];

                    auto const lhs_pos_data = get_loc(lhs_id);
                    auto const rhs_pos_data = get_loc(rhs_id);

                    bool lhs_is_new = (lhs_pos == passive_id<std::size_t>);
                    bool rhs_is_new = (rhs_pos == passive_id<std::size_t>);
                    bool lhs_inplace = lhs_is_new && std::get<0>(lhs_pos_data);
                    bool rhs_inplace = !lhs_inplace && rhs_is_new && std::get<0>(rhs_pos_data);

                    if (!lhs_inplace) {
                        copy_mul(res_pos, lhs_pos, std::get<1>(lhs_pos_data), mul_1);
                    }

                    if (!rhs_inplace) {
                        copy_mul(res_pos, rhs_pos, std::get<1>(rhs_pos_data), mul_2);
                    }

                    if (lhs_inplace) {
                        // res id should now be lhs id, avoiding a copy and a potential buffer increase
                        lhs_pos = res_pos;
                        this->pos.push_back(lhs_pos);
                        this->values.push_back(mul_1);
                        this->lossy_op.push_back(LossyOpCode::MUL_INPLACE);
                    }
                    else if (rhs_inplace) {
                        // res id should now be rhs id, avoiding a copy and a potential buffer increase
                        rhs_pos = res_pos;
                        this->pos.push_back(rhs_pos);
                        this->values.push_back(mul_2);
                        this->lossy_op.push_back(LossyOpCode::MUL_INPLACE);
                    }
                    else {
                        // don't forget to free res_id from the buffer!
                        buffer_free_positions.push_back(res_pos);
                    }
                }
            }
        }

        if constexpr (Reset) {
            data.ops.resize(to);
            data.vals.resize(val_idx_start);
            data.ids.resize(id_idx_start);
            data.next_id = to;

            this->node_location_on_buffer.resize(to);

            from_prev = 0;
            to_prev = 0;
        }

        for (auto& b : buffers) {
            b.values.resize(b.size * this->m_num_lanes, 0.);
        }
        auto& buffer = buffers.back().values;

        std::size_t val_idx_l = 0;
        std::size_t id_idx_l = 0;
        std::size_t on_which_buffer_idx = 0;

        for (std::size_t lossy_op_idx = 0; lossy_op_idx < this->lossy_op.size(); ++lossy_op_idx) {
            auto const& op = this->lossy_op[lossy_op_idx];

            switch (op) {
                case LossyOpCode::COPY: {
                    std::uint8_t const which = this->on_which_buffer[on_which_buffer_idx++];
                    std::size_t const out_pos = this->pos[id_idx_l++];
                    std::size_t const in_pos = this->pos[id_idx_l++];
                    if constexpr (Vectorised) {
                        const double* src = &buffer[out_pos * this->m_num_lanes];
                        double* dest = &buffers[which].values[in_pos * this->m_num_lanes];
#pragma omp simd
                        for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                            dest[i] = src[i];
                        }
                    }
                    else {
                        buffers[which].values[in_pos] = buffer[out_pos];
                    }
                    break;
                }
                case LossyOpCode::COPY_MINUS: {
                    std::uint8_t const which = this->on_which_buffer[on_which_buffer_idx++];
                    std::size_t const out_pos = this->pos[id_idx_l++];
                    std::size_t const in_pos = this->pos[id_idx_l++];
                    if constexpr (Vectorised) {
                        const double* src = &buffer[out_pos * this->m_num_lanes];
                        double* dest = &buffers[which].values[in_pos * this->m_num_lanes];
#pragma omp simd
                        for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                            dest[i] = -src[i];
                        }
                    }
                    else {
                        buffers[which].values[in_pos] = -buffer[out_pos];
                    }
                    break;
                }
                case LossyOpCode::ADD: {
                    std::uint8_t const which = this->on_which_buffer[on_which_buffer_idx++];
                    std::size_t const out_pos = this->pos[id_idx_l++];
                    std::size_t const in_pos = this->pos[id_idx_l++];
                    if constexpr (Vectorised) {
                        const double* src = &buffer[out_pos * this->m_num_lanes];
                        double* dest = &buffers[which].values[in_pos * this->m_num_lanes];
#pragma omp simd
                        for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                            dest[i] += src[i];
                        }
                    }
                    else {
                        buffers[which].values[in_pos] += buffer[out_pos];
                    }
                    break;
                }
                case LossyOpCode::SUB: {
                    std::uint8_t const which = this->on_which_buffer[on_which_buffer_idx++];
                    std::size_t const out_pos = this->pos[id_idx_l++];
                    std::size_t const in_pos = this->pos[id_idx_l++];
                    if constexpr (Vectorised) {
                        const double* src = &buffer[out_pos * this->m_num_lanes];
                        double* dest = &buffers[which].values[in_pos * this->m_num_lanes];
#pragma omp simd
                        for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                            dest[i] -= src[i];
                        }
                    }
                    else {
                        buffers[which].values[in_pos] -= buffer[out_pos];
                    }
                    break;
                }
                case LossyOpCode::MINUS_INPLACE: {
                    std::size_t const pos = this->pos[id_idx_l++];
                    if constexpr (Vectorised) {
                        double* dest = &buffer[pos * this->m_num_lanes];
#pragma omp simd
                        for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                            dest[i] = -dest[i];
                        }
                    }
                    else {
                        buffer[pos] = -buffer[pos];
                    }
                    break;
                }
                case LossyOpCode::MUL_INPLACE: {
                    std::size_t const pos = this->pos[id_idx_l++];
                    double const multiplier = this->values[val_idx_l++];
                    if constexpr (Vectorised) {
                        double* dest = &buffer[pos * this->m_num_lanes];
#pragma omp simd
                        for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                            dest[i] *= multiplier;
                        }
                    }
                    else {
                        buffer[pos] *= multiplier;
                    }
                    break;
                }
                case LossyOpCode::MUL_ADD: {
                    std::uint8_t const which = this->on_which_buffer[on_which_buffer_idx++];
                    std::size_t const out_pos = this->pos[id_idx_l++];
                    std::size_t const in_pos = this->pos[id_idx_l++];
                    double const multiplier = this->values[val_idx_l++];
                    if constexpr (Vectorised) {
                        const double* src = &buffer[out_pos * this->m_num_lanes];
                        double* dest = &buffers[which].values[in_pos * this->m_num_lanes];
#pragma omp simd
                        for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                            dest[i] += src[i] * multiplier;
                        }
                    }
                    else {
                        buffers[which].values[in_pos] += buffer[out_pos] * multiplier;
                    }
                    break;
                }
                case LossyOpCode::MUL_SET: {
                    std::uint8_t const which = this->on_which_buffer[on_which_buffer_idx++];
                    std::size_t const out_pos = this->pos[id_idx_l++];
                    std::size_t const in_pos = this->pos[id_idx_l++];
                    double const multiplier = this->values[val_idx_l++];
                    if constexpr (Vectorised) {
                        const double* src = &buffer[out_pos * this->m_num_lanes];
                        double* dest = &buffers[which].values[in_pos * this->m_num_lanes];
#pragma omp simd
                        for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                            dest[i] = src[i] * multiplier;
                        }
                    }
                    else {
                        buffers[which].values[in_pos] = buffer[out_pos] * multiplier;
                    }
                    break;
                }
            }
        }

        return;

        std::vector<std::size_t> multiplier_loc(from - to, passive_id<std::size_t>);

        // output is: is_new, is_current, buffer_id, position
        // auto const& checkpoints_c = this->checkpoints;
        // auto get_loc = [checkpoints_c](std::size_t id) -> std::tuple<bool, std::uint8_t> {
        //     auto it = std::upper_bound(checkpoints_c.begin(), checkpoints_c.end(), id);
        //     auto buffed_id = static_cast<std::uint8_t>(std::distance(checkpoints_c.cbegin(), it) - 1);
        //     return { it == checkpoints_c.end(), buffed_id };
        // };

        auto update_loc2 = [&buffers = this->buffers](std::size_t& arg_pos, std::uint8_t buffer_id) {
            if (arg_pos == passive_id<std::size_t>) {
                auto& arg_buffer = buffers[buffer_id];
                if (arg_buffer.free_positions.empty()) {
                    arg_pos = arg_buffer.size;
                    ++arg_buffer.size;
                }
                else {
                    arg_pos = arg_buffer.free_positions.back();
                    arg_buffer.free_positions.pop_back();
                }
            }
        };

        auto copy_add = [this, update_loc](std::size_t res_pos, std::size_t& arg_pos, std::uint8_t buffer_id) {
            bool arg_is_new = (arg_pos == passive_id<std::size_t>);
            update_loc(res_pos, arg_pos, buffer_id);

            if (arg_is_new) {
                // this is a new value, we NEED to override
                this->lossy_op.push_back(LossyOpCode::COPY);
            }
            else {
                this->lossy_op.push_back(LossyOpCode::ADD);
            }
            this->values.push_back(1.0);
        };

        auto copy_sub = [this, update_loc](std::size_t res_pos, std::size_t& arg_pos, std::uint8_t buffer_id) {
            bool arg_is_new = (arg_pos == passive_id<std::size_t>);
            update_loc(res_pos, arg_pos, buffer_id);

            if (arg_is_new) {
                // this is a new value, we NEED to override
                this->lossy_op.push_back(LossyOpCode::COPY_MINUS);
            }
            else {
                this->lossy_op.push_back(LossyOpCode::SUB);
            }
            this->values.push_back(-1.0);
        };

        auto update_univariate = [this, &multiplier_loc, to, get_loc_buffer, copy_mul, get_loc, &buffer_free_positions](
                                   std::size_t const arg_id, std::size_t const res_id, double const der_local_1) {
            auto const res_pos = this->node_location_on_buffer[res_id];
            this->node_location_on_buffer[res_id] = passive_id<std::size_t>;
            std::size_t& arg_pos = this->node_location_on_buffer[arg_id];

            auto const arg_pos_data = get_loc(arg_id);
            bool arg_is_new = (arg_pos == passive_id<std::size_t>);
            bool arg_inplace = arg_is_new && std::get<0>(arg_pos_data);

            auto const multiplier_location = multiplier_loc[res_id - to];
            if (multiplier_location != passive_id<std::size_t>) {
                if (arg_inplace) {
                    arg_pos = res_pos;
                }
                else {
                    get_loc_buffer(arg_pos, std::get<1>(arg_pos_data));
                }
                this->pos[(multiplier_location * 2) + 1] = arg_pos;
                this->values[multiplier_location] *= der_local_1;
                this->on_which_buffer[multiplier_location] = std::get<1>(arg_pos_data);
                auto const original_op = this->lossy_op[multiplier_location];
                if (original_op == LossyOpCode::COPY) {
                    if (arg_is_new) {
                        this->lossy_op[multiplier_location] = LossyOpCode::MUL_SET;
                    }
                    else {
                        this->lossy_op[multiplier_location] = LossyOpCode::MUL_ADD;
                    }
                }
                else if (original_op == LossyOpCode::MUL_SET) {
                    if (arg_is_new) {
                        this->lossy_op[multiplier_location] = LossyOpCode::MUL_SET;
                    }
                    else {
                        this->lossy_op[multiplier_location] = LossyOpCode::MUL_ADD;
                    }
                }
                else if (original_op == LossyOpCode::MUL_INPLACE) {
                    if (arg_is_new) {
                        if (arg_inplace) {
                            this->lossy_op[multiplier_location] = LossyOpCode::MUL_INPLACE;
                        }
                        else {
                            this->lossy_op[multiplier_location] = LossyOpCode::MUL_SET;
                        }
                    }
                    else {
                        this->lossy_op[multiplier_location] = LossyOpCode::MUL_ADD;
                    }
                }

                if (arg_id >= to) {
                    if (multiplier_loc[arg_id - to] != passive_id<std::size_t>) {
                        // we CUT the multiplication chain, because this node has multiple outputs
                        multiplier_loc[arg_id - to] = passive_id<std::size_t>;
                    }
                    else {
                        multiplier_loc[arg_id - to] = multiplier_loc[res_id - to];
                    }
                }
            }
            else {
                if (arg_inplace) {
                    // res id should now be arg id, avoiding a copy and a potential buffer increase
                    arg_pos = res_pos;
                    this->pos.push_back(res_pos);
                    this->pos.push_back(arg_pos);
                    this->values.push_back(der_local_1);
                    this->on_which_buffer.push_back(std::get<1>(arg_pos_data));
                    this->lossy_op.push_back(LossyOpCode::MUL_INPLACE);

                    if (arg_id >= to) {
                        multiplier_loc[arg_id - to] = this->values.size() - 1;
                    }
                }
                else {
                    copy_mul(res_pos, arg_pos, std::get<1>(arg_pos_data), der_local_1);
                    buffer_free_positions.push_back(res_pos);
                }
            }
        };

        enum class mul_type : std::uint8_t {
            ANY,
            ONE,
            MINUS_ONE,
        };

        std::vector<bool> use_op(from - to);
        std::vector<mul_type> values_type;

        auto add_pos_val = [this, &multiplier_loc, &use_op, &values_type, to](
                             std::size_t res_id, std::size_t arg_id, double val, mul_type m) {
            use_op[arg_id - to] = true;

            this->values.push_back(val);
            values_type.push_back(m);
            this->pos.push_back(res_id);
            this->pos.push_back(arg_id);

            auto& mult_loc = multiplier_loc[arg_id - to];
            if (mult_loc == passive_id<std::size_t>) {
                mult_loc = this->values.size() - 1;
            }
            else {
                mult_loc = passive_id<std::size_t>;
            }
        };

        auto add_or_chain_pos_val = [this, &multiplier_loc, &use_op, &values_type, to](
                                      std::size_t res_id, std::size_t arg_id, double val, mul_type m) {
            auto& mult_loc_res = multiplier_loc[res_id - to];
            if (mult_loc_res != passive_id<std::size_t>) {
                this->values[mult_loc_res] *= val;
                if (m == mul_type::ANY || values_type[mult_loc_res] == mul_type::ANY) {
                    values_type[mult_loc_res] = mul_type::ANY;
                }
                else if (m == values_type[mult_loc_res]) {
                    values_type[mult_loc_res] = mul_type::ONE;
                }
                else {
                    values_type[mult_loc_res] = mul_type::MINUS_ONE;
                }

                this->pos[(mult_loc_res * 2) + 1] = arg_id;
            }
            else {
                use_op[arg_id - to] = true;

                this->values.push_back(val);
                values_type.push_back(m);
                this->pos.push_back(res_id);
                this->pos.push_back(arg_id);
            }

            auto& mult_loc = multiplier_loc[arg_id - to];
            if (mult_loc == passive_id<std::size_t>) {
                mult_loc = this->values.size() - 1;
            }
            else {
                mult_loc = passive_id<std::size_t>;
            }
        };

        for (std::size_t op_idx = from; op_idx-- > to;) {
            OpCode const& op = ops[op_idx];
            bool const use_this_op =
              this->node_location_on_buffer[op_idx] != passive_id<std::size_t> || use_op[op_idx - to];

            switch (op) {
                case OpCode::REG_INPUT: {
                    id_idx -= 1;
                    break;
                }
                case OpCode::REG_OUTPUT: {
                    id_idx -= 2;
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        add_pos_val(res_id, arg_id, 1.0, mul_type::ONE);
                    }
                    break;
                }
                case OpCode::ADD: {
                    id_idx -= 3;
                    if (use_this_op) {
                        std::size_t const lhs_id = ids[id_idx];
                        std::size_t const rhs_id = ids[id_idx + 1];
                        std::size_t const res_id = ids[id_idx + 2];
                        add_pos_val(res_id, lhs_id, 1.0, mul_type::ONE);
                        add_pos_val(res_id, rhs_id, 1.0, mul_type::ONE);
                    }
                    break;
                }
                case OpCode::SUB: {
                    id_idx -= 3;
                    if (use_this_op) {
                        std::size_t const lhs_id = ids[id_idx];
                        std::size_t const rhs_id = ids[id_idx + 1];
                        std::size_t const res_id = ids[id_idx + 2];
                        add_pos_val(res_id, lhs_id, 1.0, mul_type::ONE);
                        add_pos_val(res_id, rhs_id, -1.0, mul_type::MINUS_ONE);
                    }
                    break;
                }
                case OpCode::MUL: {
                    val_idx -= 2;
                    id_idx -= 3;
                    if (use_this_op) {
                        double const lhs_val = vals[val_idx];
                        double const rhs_val = vals[val_idx + 1];
                        std::size_t const lhs_id = ids[id_idx];
                        std::size_t const rhs_id = ids[id_idx + 1];
                        std::size_t const res_id = ids[id_idx + 2];
                        add_pos_val(res_id, lhs_id, rhs_val, mul_type::ANY);
                        add_pos_val(res_id, rhs_id, lhs_val, mul_type::ANY);
                    }
                    break;
                }
                case OpCode::ADD_C: {
                    id_idx -= 2;
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        add_or_chain_pos_val(res_id, arg_id, 1.0, mul_type::ONE);
                    }
                    break;
                }
                case OpCode::SUB_C: {
                    id_idx -= 2;
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        add_or_chain_pos_val(res_id, arg_id, -1.0, mul_type::MINUS_ONE);
                    }
                    break;
                }
                case OpCode::MUL_C: {
                    val_idx -= 1;
                    id_idx -= 2;
                    if (use_this_op) {
                        double const der_local_1 = vals[val_idx];
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        add_or_chain_pos_val(res_id, arg_id, der_local_1, mul_type::ANY);
                    }
                    break;
                }
                case OpCode::NORM: {
                    val_idx -= 1;
                    id_idx -= 2;
                    if (use_this_op) {
                        double const der_local_1 = 2.0 * vals[val_idx];
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        add_or_chain_pos_val(res_id, arg_id, der_local_1, mul_type::ANY);
                    }
                    break;
                }
                case OpCode::INV: {
                    val_idx -= 1;
                    id_idx -= 2;
                    if (use_this_op) {
                        double const der_local_1 = -vals[val_idx] * vals[val_idx];
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        add_or_chain_pos_val(res_id, arg_id, der_local_1, mul_type::ANY);
                    }
                    break;
                }
                case OpCode::ABS: {
                    val_idx -= 1;
                    id_idx -= 2;
                    if (use_this_op) {
                        double const der_local_1 = std::copysign(1.0, vals[val_idx]);
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        add_or_chain_pos_val(res_id, arg_id, der_local_1, mul_type::ANY);
                    }
                    break;
                }
                case OpCode::EXP: {
                    val_idx -= 1;
                    id_idx -= 2;
                    if (use_this_op) {
                        double const der_local_1 = vals[val_idx];
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        add_or_chain_pos_val(res_id, arg_id, der_local_1, mul_type::ANY);
                    }
                    break;
                }
                case OpCode::LOG: {
                    val_idx -= 1;
                    id_idx -= 2;
                    if (use_this_op) {
                        double const der_local_1 = 1.0 / vals[val_idx];
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        add_or_chain_pos_val(res_id, arg_id, der_local_1, mul_type::ANY);
                    }
                    break;
                }
                case OpCode::ERF: {
                    val_idx -= 1;
                    id_idx -= 2;
                    if (use_this_op) {
                        constexpr double two_over_root_pi = 2. * std::numbers::inv_sqrtpi_v<double>;
                        double const der_local_1 = std::exp(-vals[val_idx] * vals[val_idx]) * two_over_root_pi;
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        add_or_chain_pos_val(res_id, arg_id, der_local_1, mul_type::ANY);
                    }
                    break;
                }
                case OpCode::ERFC: {
                    val_idx -= 1;
                    id_idx -= 2;
                    if (use_this_op) {
                        constexpr double minus_two_over_root_pi = -2. * std::numbers::inv_sqrtpi_v<double>;
                        double const der_local_1 = std::exp(-vals[val_idx] * vals[val_idx]) * minus_two_over_root_pi;
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        add_or_chain_pos_val(res_id, arg_id, der_local_1, mul_type::ANY);
                    }
                    break;
                }
                case OpCode::COS: {
                    val_idx -= 2;
                    id_idx -= 2;
                    if (use_this_op) {
                        double const der_local_1 = -std::sin(vals[val_idx]);
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        add_or_chain_pos_val(res_id, arg_id, der_local_1, mul_type::ANY);
                    }
                    break;
                }
                case OpCode::SQRT: {
                    val_idx -= 2;
                    id_idx -= 2;
                    if (use_this_op) {
                        double const one_over_in = 1. / vals[val_idx];
                        double const der_local_1 = 0.5 * vals[val_idx + 1] * one_over_in;
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        add_or_chain_pos_val(res_id, arg_id, der_local_1, mul_type::ANY);
                    }
                    break;
                }
                case OpCode::POW_C: {
                    val_idx -= 3;
                    id_idx -= 2;
                    if (use_this_op) {
                        double const one_over_in = 1. / vals[val_idx];
                        double const der_local_1 = vals[val_idx + 2] * vals[val_idx + 1] * one_over_in;
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        add_or_chain_pos_val(res_id, arg_id, der_local_1, mul_type::ANY);
                    }
                    break;
                }
            }
        }

        std::size_t values_idx = 0;
        std::size_t pos_idx = 0;
        std::size_t values_idx2 = 0;
        while (values_idx < this->values.size()) {
            bool is_not_last = values_idx < (this->values.size() - 1);
            bool is_bivariate = is_not_last && this->pos[(values_idx * 2)] == this->pos[(values_idx * 2) + 2];
            if (is_bivariate) {
                bool is_actually_univariate = this->pos[(values_idx * 2) + 1] == this->pos[(values_idx * 2) + 3];
                if (is_actually_univariate) {
                    std::size_t const res_id = this->pos[(values_idx * 2)];
                    std::size_t const arg_id = this->pos[(values_idx * 2) + 1];

                    auto const res_pos = this->node_location_on_buffer[res_id];
                    this->node_location_on_buffer[res_id] = passive_id<std::size_t>;
                    std::size_t& arg_pos = this->node_location_on_buffer[arg_id];

                    auto const arg_pos_data = get_loc(arg_id);
                    bool arg_is_new = (arg_pos == passive_id<std::size_t>);
                    bool arg_inplace = arg_is_new && std::get<0>(arg_pos_data);

                    double const der_local_1 = this->values[values_idx] + this->values[values_idx + 1];
                    if (arg_inplace) {
                        // res id should now be arg id, avoiding a copy and a potential buffer increase
                        arg_pos = res_pos;
                        this->pos[pos_idx++] = res_pos;
                        this->lossy_op.push_back(LossyOpCode::MUL_INPLACE);
                        this->values[values_idx2++] = der_local_1;
                    }
                    else {
                        throw;
                        if (arg_is_new) {
                            // this is a new value, we NEED to override
                            this->lossy_op.push_back(LossyOpCode::COPY);
                        }
                        else {
                            this->lossy_op.push_back(LossyOpCode::ADD);
                        }

                        // copy_mul(res_pos, arg_pos, std::get<1>(arg_pos_data), der_local_1);
                        buffer_free_positions.push_back(res_pos);
                    }
                }
                else {
                    std::size_t const res_id = this->pos[(values_idx * 2)];
                    std::size_t const lhs_id = this->pos[(values_idx * 2) + 1];
                    std::size_t const rhs_id = this->pos[(values_idx * 2) + 3];
                    double const lhs_der_val = this->values[values_idx];
                    double const rhs_der_val = this->values[values_idx + 1];

                    auto const res_pos = node_location_on_buffer[res_id];
                    node_location_on_buffer[res_id] = passive_id<std::size_t>;
                    std::size_t& lhs_pos = node_location_on_buffer[lhs_id];
                    std::size_t& rhs_pos = node_location_on_buffer[rhs_id];

                    auto const lhs_pos_data = get_loc(lhs_id);
                    auto const rhs_pos_data = get_loc(rhs_id);

                    bool lhs_is_new = (lhs_pos == passive_id<std::size_t>);
                    bool rhs_is_new = (rhs_pos == passive_id<std::size_t>);
                    bool lhs_inplace = lhs_is_new && std::get<0>(lhs_pos_data);
                    bool rhs_inplace = !lhs_inplace && rhs_is_new && std::get<0>(rhs_pos_data);

                    if (!lhs_inplace) {
                        throw;
                    }

                    if (!rhs_inplace) {
                        if (rhs_is_new) {
                            // this is a new value, we NEED to override
                            this->lossy_op.push_back(LossyOpCode::MUL_SET);
                        }
                        else {
                            this->lossy_op.push_back(LossyOpCode::MUL_ADD);
                        }
                        update_loc2(rhs_pos, std::get<1>(rhs_pos_data));
                        this->pos[pos_idx++] = res_pos;
                        this->pos[pos_idx++] = rhs_pos;
                        this->on_which_buffer.push_back(std::get<1>(rhs_pos_data));
                        this->values[values_idx2++] = rhs_der_val;
                        throw;
                    }

                    if (lhs_inplace) {
                        throw;
                        // res id should now be lhs id, avoiding a copy and a potential buffer increase
                        lhs_pos = res_pos;
                    }
                    else if (rhs_inplace) {
                        throw;
                        // res id should now be rhs id, avoiding a copy and a potential buffer increase
                        rhs_pos = res_pos;
                    }
                    else {
                        // don't forget to free res_id from the buffer!
                        buffer_free_positions.push_back(res_pos);
                    }
                }

                ++values_idx;
                ++values_idx;
            }
            else {
                std::size_t const res_id = this->pos[(values_idx * 2)];
                std::size_t const arg_id = this->pos[(values_idx * 2) + 1];

                auto const res_pos = this->node_location_on_buffer[res_id];
                this->node_location_on_buffer[res_id] = passive_id<std::size_t>;
                std::size_t& arg_pos = this->node_location_on_buffer[arg_id];

                auto const arg_pos_data = get_loc(arg_id);
                bool arg_is_new = (arg_pos == passive_id<std::size_t>);
                bool arg_inplace = arg_is_new && std::get<0>(arg_pos_data);

                double const der_local_1 = this->values[values_idx];
                if (arg_inplace) {
                    // res id should now be arg id, avoiding a copy and a potential buffer increase
                    arg_pos = res_pos;
                    if (values_type[values_idx] == mul_type::ONE) {
                        // do nothing
                    }
                    else if (values_type[values_idx] == mul_type::MINUS_ONE) {
                        this->pos[pos_idx++] = res_pos;
                        this->lossy_op.push_back(LossyOpCode::MINUS_INPLACE);
                        ++values_idx2;
                    }
                    else {
                        this->pos[pos_idx++] = res_pos;
                        this->lossy_op.push_back(LossyOpCode::MUL_INPLACE);
                        this->values[values_idx2++] = der_local_1;
                    }
                }
                else {
                    if (values_type[values_idx] == mul_type::ONE) {
                        if (arg_is_new) {
                            throw;
                            // this is a new value, we NEED to override
                            this->lossy_op.push_back(LossyOpCode::COPY);
                        }
                        else {
                            throw;
                            this->lossy_op.push_back(LossyOpCode::ADD);
                        }
                    }
                    else if (values_type[values_idx] == mul_type::MINUS_ONE) {
                        if (arg_is_new) {
                            throw;
                            // this is a new value, we NEED to override
                            this->lossy_op.push_back(LossyOpCode::COPY_MINUS);
                        }
                        else {
                            throw;
                            this->lossy_op.push_back(LossyOpCode::SUB);
                        }
                    }
                    else {
                        if (arg_is_new) {
                            throw;
                            // this is a new value, we NEED to override
                            this->lossy_op.push_back(LossyOpCode::MUL_SET);
                        }
                        else {
                            throw;
                            this->lossy_op.push_back(LossyOpCode::MUL_ADD);
                        }
                    }

                    // throw;
                    // if (arg_is_new) {
                    //     // this is a new value, we NEED to override
                    //     this->lossy_op.push_back(LossyOpCode::COPY);
                    // }
                    // else {
                    //     this->lossy_op.push_back(LossyOpCode::ADD);
                    // }

                    // copy_mul(res_pos, arg_pos, std::get<1>(arg_pos_data), der_local_1);
                    buffer_free_positions.push_back(res_pos);
                }

                ++values_idx;
            }
        }
    }

    if constexpr (Reset) {
        data.ops.resize(to);
        data.vals.resize(val_idx);
        data.ids.resize(id_idx);
        data.next_id = to;

        this->node_location_on_buffer.resize(to);

        from_prev = 0;
        to_prev = 0;
    }

    for (auto& b : buffers) {
        b.values.resize(b.size * this->m_num_lanes, 0.);
    }
    auto& buffer = buffers.back().values;

    std::size_t id_idx_l = 0;

    for (std::size_t lossy_op_idx = 0; lossy_op_idx < this->lossy_op.size(); ++lossy_op_idx) {
        auto const& op = this->lossy_op[lossy_op_idx];

        std::uint8_t const which = this->on_which_buffer[lossy_op_idx];
        std::size_t const out_pos = this->pos[id_idx_l++];
        std::size_t const in_pos = this->pos[id_idx_l++];
        double const multiplier = this->values[lossy_op_idx];
        switch (op) {
            case LossyOpCode::COPY:
            case LossyOpCode::COPY_MINUS:
            case LossyOpCode::MUL_SET:
            case LossyOpCode::MINUS_INPLACE:
            case LossyOpCode::MUL_INPLACE: {
                if constexpr (Vectorised) {
                    const double* src = &buffer[out_pos * this->m_num_lanes];
                    double* dest = &buffers[which].values[in_pos * this->m_num_lanes];
#pragma omp simd
                    for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                        dest[i] = src[i] * multiplier;
                    }
                }
                else {
                    buffers[which].values[in_pos] = buffer[out_pos] * multiplier;
                }
                break;
            }
            case LossyOpCode::ADD:
            case LossyOpCode::SUB:
            case LossyOpCode::MUL_ADD: {
                if constexpr (Vectorised) {
                    const double* src = &buffer[out_pos * this->m_num_lanes];
                    double* dest = &buffers[which].values[in_pos * this->m_num_lanes];
#pragma omp simd
                    for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                        dest[i] += src[i] * multiplier;
                    }
                }
                else {
                    buffers[which].values[in_pos] += buffer[out_pos] * multiplier;
                }
                break;
            }
        }
    }

    if constexpr (Reset) {
        if (!checkpoints.empty() && to == checkpoints.back()) {
            // buffers.pop_back();
            buffers.back() = buffer_t{};
        }
    }
}

} // namespace adhoc

#endif // ADHOC_BACKPROPAGATOR1LOSSYCOMPRESSED_HPP
