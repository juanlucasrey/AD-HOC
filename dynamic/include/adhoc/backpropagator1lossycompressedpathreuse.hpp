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

#ifndef ADHOC_BACKPROPAGATOR1LOSSYCOMPRESSEDPATHREUSE_HPP
#define ADHOC_BACKPROPAGATOR1LOSSYCOMPRESSEDPATHREUSE_HPP

#include "hash.hpp"
#include "passive_id.hpp"
#include "position_impl.hpp"
#include "tape_data.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numbers>
#include <vector>

namespace adhoc {

template<class Float, bool Vectorised = false>
class BackPropagatorLossyCompressedPathReuse {
  private:
    std::size_t m_num_lanes{ 1 };

    std::vector<std::size_t> node_location_on_buffer;

    struct buffer_t {
        std::vector<double> values;
        std::size_t size{ 0 };
        std::vector<std::size_t> free_positions;
    };

    std::vector<std::size_t> checkpoints{ 0 };
    std::vector<buffer_t> buffers{ buffer_t{} };

    // lossy tape
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

    struct LossyTape {

        std::vector<bool> op_multiply;
        std::vector<bool> op_invert;
        std::vector<std::size_t> pos_multiplier;
        std::vector<bool> bivariate_consolidate;
        std::vector<bool> univariate_consolidate;

        std::vector<bool> multiplier_keep_alive;

        std::vector<std::uint8_t> on_which_buffer;
        std::vector<std::size_t> pos;
        std::vector<LossyOpCode> lossy_op;
        std::vector<bool> invert_mult;
        std::vector<bool> use_op;
        std::size_t buffer_size{ 0 };
    };

    std::vector<std::map<std::size_t, LossyTape> > hash_to_lossy_tape;

  public:
    explicit BackPropagatorLossyCompressedPathReuse() = default;

    void set_checkpoint(std::size_t ops_size)
    {
        if (this->checkpoints.back() > ops_size) {
            // should not happen
            throw;
        }

        if (this->checkpoints.back() != ops_size) {
            this->checkpoints.push_back(ops_size);
            this->buffers.push_back(buffer_t{});
            this->hash_to_lossy_tape.push_back(std::map<std::size_t, LossyTape>{});
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
            auto it = std::upper_bound(this->checkpoints.cbegin(), this->checkpoints.cend(), var_id);
            auto& var_buffer = this->buffers[std::distance(this->checkpoints.cbegin(), it) - 1];
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
            auto it = std::upper_bound(this->checkpoints.cbegin(), this->checkpoints.cend(), var_id);
            auto& var_buffer = this->buffers[std::distance(this->checkpoints.cbegin(), it) - 1];
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

            auto it = std::upper_bound(this->checkpoints.cbegin(), this->checkpoints.cend(), var_id);
            auto buffed_id = static_cast<std::uint8_t>(std::distance(this->checkpoints.cbegin(), it) - 1);
            this->buffers[buffed_id].values[(var_pos * this->m_num_lanes) + lane] = deriv;
        }
        else {
            throw;
        }
    }

    void set_derivative(std::size_t /* var_id1 */, std::size_t /* var_id2 */, double /* deriv */) { throw; }

    auto get_derivative(std::size_t var_id, std::size_t lane) const -> double
    {
        if (lane < this->m_num_lanes) {
            std::size_t const var_pos = this->node_location_on_buffer[var_id];
            if (var_pos == passive_id<std::size_t>) {
                // this derivative is not on buffer.
                // this is probably nor an input nor an output.
                throw;
            }

            auto it = std::upper_bound(this->checkpoints.cbegin(), this->checkpoints.cend(), var_id);
            return this->buffers[std::distance(this->checkpoints.cbegin(), it) - 1]
              .values[(var_pos * this->m_num_lanes) + lane];
        }

        throw;
        return 0.;
    }

    auto get_derivative(std::size_t /* var_id1 */, std::size_t /* var_id2 */, std::size_t /* lane */) const -> double
    {
        return 0.;
    }

    void clear() {}

    void zero_adjoints()
    {
        for (auto& b : this->buffers) {
            std::fill(b.values.begin(), b.values.end(), 0.0);
        }
    }

    auto size_of(bool capacity = false) const -> std::size_t
    {
        std::size_t size = 0;
        size += sizeof(std::size_t); // m_num_lanes
        size += sizeof(std::size_t) * (capacity ? node_location_on_buffer.capacity() : node_location_on_buffer.size());
        size += sizeof(std::size_t) * (capacity ? this->checkpoints.capacity() : this->checkpoints.size());
        for (const auto& buffer : this->buffers) {
            size += sizeof(double) * (capacity ? buffer.values.capacity() : buffer.values.size());
            size += sizeof(std::size_t) * (capacity ? buffer.free_positions.capacity() : buffer.free_positions.size());
        }
        return size;
    }

    template<bool Reset = false, bool ResetInPlace = false, bool Log = false>
    void backpropagate_to(PositionImpl const& pos, TapeData& data);
};

template<class Float, bool Vectorised>
template<bool Reset, bool ResetInPlace, bool Log>
void
BackPropagatorLossyCompressedPathReuse<Float, Vectorised>::backpropagate_to(PositionImpl const& pos, TapeData& data)
{
    auto convert_to_lossy_tape = [](PositionImpl const& pos,
                                    TapeData const& data,
                                    std::vector<std::size_t>& node_location_on_buffer,
                                    std::vector<std::size_t> const& checkpoints,
                                    std::vector<buffer_t>& buffers,
                                    std::size_t num_lanes) -> LossyTape {
        std::size_t const to = pos.op_position;
        LossyTape result;

        std::size_t from = data.next_id;
        const auto& ops = data.ops;
        const auto& ids = data.ids;

        auto& op_multiply = result.op_multiply;
        auto& op_invert = result.op_invert;
        auto& pos_multiplier = result.pos_multiplier;
        auto& bivariate_consolidate = result.bivariate_consolidate;
        auto& univariate_consolidate = result.univariate_consolidate;

        std::vector<std::size_t> number_dependents(from - to);

        // LOOP 1: backward, to count number of dependents for each node and detect which nodes are active
        std::size_t id_idx = ids.size();
        for (std::size_t op_idx = from; op_idx-- > to;) {
            OpCode const& op = ops[op_idx];
            bool const use_this_op =
              node_location_on_buffer[op_idx] != passive_id<std::size_t> || (number_dependents[op_idx - to] > 0);

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

                    if (node_location_on_buffer[op_idx] != passive_id<std::size_t>) {
                        // if the node is already active, it measn it had at least on precedent
                        // on other buffers
                        number_dependents[op_idx - to] += 1;
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

                    if (node_location_on_buffer[op_idx] != passive_id<std::size_t>) {
                        // if the node is already active, it measn it had at least on precedent
                        // on other buffers
                        number_dependents[op_idx - to] += 1;
                    }
                    break;
                }
            }
        }

        std::vector<std::size_t> number_dependents2(number_dependents.begin(), number_dependents.end());

        std::vector<std::size_t> multiplier_loc_from;
        std::vector<bool> multiplier_keep_alive;
        std::vector<std::size_t> multiplier_origin((from - to) * 2, passive_id<std::size_t>);

        // LOOP 2: forward, to calculate multipliers after compressing induced paths
        enum class mul_type : std::uint8_t {
            ANY,
            ONE,
            MINUS_ONE,
        };
        std::vector<mul_type> values_type;

        auto combine_mul_type = [](mul_type a, mul_type b) -> mul_type {
            if (a == mul_type::ANY || b == mul_type::ANY) {
                return mul_type::ANY;
            }
            if (a == b) {
                return mul_type::ONE;
            }
            return mul_type::MINUS_ONE;
        };

        auto multiplier_set_incoming = [&]<mul_type M>(std::size_t const pos) {
            if constexpr (M == mul_type::ONE) {
                // buffer_multipliers.values[pos] = 1.0;
                values_type[pos] = mul_type::ONE;
            }
            else if constexpr (M == mul_type::MINUS_ONE) {
                // buffer_multipliers.values[pos] = -1.0;
                values_type[pos] = mul_type::MINUS_ONE;
            }
            else {
                static_assert(M == mul_type::ANY, "Invalid multiplier type");
                // buffer_multipliers.values[pos] = multiplier;
                op_multiply.push_back(false);
                pos_multiplier.push_back(pos);
                values_type[pos] = mul_type::ANY;
            }
        };

        auto multiplier_multiply_incoming = [&, combine_mul_type]<mul_type M>(std::size_t const pos) {
            if constexpr (M == mul_type::ONE) {
                // no OP
                values_type[pos] = combine_mul_type(values_type[pos], mul_type::ONE);
            }
            else if constexpr (M == mul_type::MINUS_ONE) {
                // buffer_multipliers.values[pos] = -buffer_multipliers.values[pos];
                op_multiply.push_back(true);
                op_invert.push_back(true);
                pos_multiplier.push_back(pos);
                values_type[pos] = combine_mul_type(values_type[pos], mul_type::MINUS_ONE);
            }
            else {
                static_assert(M == mul_type::ANY, "Invalid multiplier type");
                // buffer_multipliers.values[pos] *= multiplier;
                op_multiply.push_back(true);
                op_invert.push_back(false);
                pos_multiplier.push_back(pos);
                values_type[pos] = mul_type::ANY;
            }
        };

        auto multiplier_add = [&](std::size_t const pos1, std::size_t const /* pos2 */) {
            // buffer_multipliers.values[pos1] += buffer_multipliers.values[pos2];
            values_type[pos1] = mul_type::ANY;
        };

        auto multiplier_multiply = [&](std::size_t const pos1, std::size_t const /* pos2 */) {
            // buffer_multipliers.values[pos1] *= buffer_multipliers.values[pos2];
            values_type[pos1] = mul_type::ANY;
        };

        struct buffer_no_values {
            // to be removed
            std::vector<double> values;

            std::size_t size{ 0 };
            std::vector<std::size_t> free_positions;
        };

        buffer_no_values buffer_multipliers;

        auto get_mult_loc = [&]() -> std::size_t {
            if (buffer_multipliers.free_positions.empty()) {
                std::size_t result = buffer_multipliers.size;
                ++buffer_multipliers.size;
                // buffer_multipliers.values.resize(buffer_multipliers.values.size() + 1);
                values_type.resize(buffer_multipliers.size);
                multiplier_loc_from.resize(buffer_multipliers.size);
                multiplier_keep_alive.resize(buffer_multipliers.size);
                return result;
            }

            std::size_t result = buffer_multipliers.free_positions.back();
            buffer_multipliers.free_positions.pop_back();
            return result;
        };

        auto update_univariate =
          [&]<mul_type M, bool KeepAlive = false>(std::size_t const arg_id,
                                                  std::size_t const res_id /* , double const multiplier = 0 */) {
              bool const arg_is_induced_path = (arg_id >= to) && (number_dependents[arg_id - to] == 1) &&
                                               (multiplier_origin[(arg_id - to) * 2] != passive_id<std::size_t>) &&
                                               (multiplier_origin[((arg_id - to) * 2) + 1] == passive_id<std::size_t>);

              auto& mult_origin_res = multiplier_origin[(res_id - to) * 2];

              if (arg_is_induced_path) {
                  auto& mul_origin_arg = multiplier_origin[(arg_id - to) * 2];
                  multiplier_multiply_incoming.template operator()<M>(mul_origin_arg);
                  multiplier_keep_alive[mul_origin_arg] = KeepAlive;

                  std::swap(mul_origin_arg, mult_origin_res);
                  --number_dependents[arg_id - to];
              }
              else {
                  mult_origin_res = get_mult_loc();
                  multiplier_set_incoming.template operator()<M>(mult_origin_res);
                  multiplier_loc_from[mult_origin_res] = arg_id;
                  multiplier_keep_alive[mult_origin_res] = KeepAlive;
              }
          };

        auto update_bivariate = [&]<mul_type M1, mul_type M2>(
                                  std::size_t const lhs_id, std::size_t const rhs_id, std::size_t const res_id) {
            bool const lhs_is_induced_path = (lhs_id >= to) && (number_dependents[lhs_id - to] == 1) &&
                                             (multiplier_origin[(lhs_id - to) * 2] != passive_id<std::size_t>) &&
                                             (multiplier_origin[((lhs_id - to) * 2) + 1] == passive_id<std::size_t>);
            bool const rhs_is_induced_path = (rhs_id >= to) && (number_dependents[rhs_id - to] == 1) &&
                                             (multiplier_origin[(rhs_id - to) * 2] != passive_id<std::size_t>) &&
                                             (multiplier_origin[((rhs_id - to) * 2) + 1] == passive_id<std::size_t>);

            std::size_t multiplier_loc_lhs = passive_id<std::size_t>;
            if (lhs_is_induced_path) {
                auto& mul_origin_arg = multiplier_origin[(lhs_id - to) * 2];
                multiplier_multiply_incoming.template operator()<M1>(mul_origin_arg);
                std::swap(multiplier_loc_lhs, mul_origin_arg);
                --number_dependents[lhs_id - to];
            }
            else {
                multiplier_loc_lhs = get_mult_loc();
                multiplier_set_incoming.template operator()<M1>(multiplier_loc_lhs);
                multiplier_loc_from[multiplier_loc_lhs] = lhs_id;
                multiplier_keep_alive[multiplier_loc_lhs] = false;
            }

            std::size_t multiplier_loc_rhs = passive_id<std::size_t>;
            if (rhs_is_induced_path) {
                auto& mul_origin_arg = multiplier_origin[(rhs_id - to) * 2];
                multiplier_multiply_incoming.template operator()<M2>(mul_origin_arg);
                std::swap(multiplier_loc_rhs, mul_origin_arg);
                --number_dependents[rhs_id - to];
            }
            else {
                multiplier_loc_rhs = get_mult_loc();
                multiplier_set_incoming.template operator()<M2>(multiplier_loc_rhs);
                multiplier_loc_from[multiplier_loc_rhs] = rhs_id;
                multiplier_keep_alive[multiplier_loc_rhs] = false;
            }

            bool const has_induced_path = lhs_is_induced_path || rhs_is_induced_path;
            // there is potential for a bivariate operator on the same argument, we need to check if
            // this is the case and update the multiplier if so
            bool const bivariate_consolidate_this =
              has_induced_path && (multiplier_loc_from[multiplier_loc_lhs] == multiplier_loc_from[multiplier_loc_rhs]);
            bivariate_consolidate.push_back(bivariate_consolidate_this);

            if (bivariate_consolidate_this) {
                std::size_t const origin_id = multiplier_loc_from[multiplier_loc_lhs];
                pos_multiplier.push_back(multiplier_loc_lhs);
                pos_multiplier.push_back(multiplier_loc_rhs);
                buffer_multipliers.free_positions.push_back(multiplier_loc_rhs);
                multiplier_add(multiplier_loc_lhs, multiplier_loc_rhs);

                multiplier_loc_from[multiplier_loc_rhs] = passive_id<std::size_t>;

                multiplier_origin[(res_id - to) * 2] = multiplier_loc_lhs;

                if (origin_id >= to) {
                    auto& number_dependents_to_update = number_dependents[origin_id - to];
                    --number_dependents_to_update;

                    bool const has_single_origin =
                      (multiplier_origin[(origin_id - to) * 2] != passive_id<std::size_t>) &&
                      (multiplier_origin[((origin_id - to) * 2) + 1] == passive_id<std::size_t>);

                    bool const univariate_consolidate_this = (number_dependents_to_update == 1) && has_single_origin;
                    univariate_consolidate.push_back(univariate_consolidate_this);
                    if (univariate_consolidate_this) {
                        // this node now has only one dependent, we can reintroduce a
                        // multiplication chain
                        auto& coming_from = multiplier_origin[(origin_id - to) * 2];
                        pos_multiplier.push_back(coming_from);
                        buffer_multipliers.free_positions.push_back(multiplier_loc_lhs);
                        multiplier_multiply(coming_from, multiplier_loc_lhs);

                        multiplier_loc_from[multiplier_loc_lhs] = passive_id<std::size_t>;

                        multiplier_origin[(res_id - to) * 2] = coming_from;
                        coming_from = passive_id<std::size_t>;
                        --number_dependents[origin_id - to];
                    }
                }
                else {
                    univariate_consolidate.push_back(false);
                }
            }
            else {
                multiplier_origin[(res_id - to) * 2] = multiplier_loc_lhs;
                multiplier_origin[((res_id - to) * 2) + 1] = multiplier_loc_rhs;
            }
        };

        // LOOP 2: forward, to calculate multipliers after compressing induced paths
        id_idx = pos.id_position;
        // std::size_t val_idx = pos.val_position;

        for (std::size_t op_idx = to; op_idx < from; ++op_idx) {
            OpCode const& op = ops[op_idx];
            bool const use_this_op = (number_dependents[op_idx - to] > 0);

            switch (op) {
                case OpCode::REG_INPUT: {
                    id_idx += 1;
                    break;
                }
                case OpCode::REG_OUTPUT: {
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        update_univariate.template operator()<mul_type::ONE, !Reset>(arg_id, res_id);
                    }
                    id_idx += 2;
                    break;
                }
                case OpCode::ADD: {
                    if (use_this_op) {
                        std::size_t const lhs_id = ids[id_idx];
                        std::size_t const rhs_id = ids[id_idx + 1];
                        std::size_t const res_id = ids[id_idx + 2];
                        update_bivariate.template operator()<mul_type::ONE, mul_type::ONE>(lhs_id, rhs_id, res_id);
                    }
                    id_idx += 3;
                    break;
                }
                case OpCode::SUB: {
                    if (use_this_op) {
                        std::size_t const lhs_id = ids[id_idx];
                        std::size_t const rhs_id = ids[id_idx + 1];
                        std::size_t const res_id = ids[id_idx + 2];
                        update_bivariate.template operator()<mul_type::ONE, mul_type::MINUS_ONE>(
                          lhs_id, rhs_id, res_id);
                    }
                    id_idx += 3;
                    break;
                }
                case OpCode::MUL: {
                    if (use_this_op) {
                        std::size_t const lhs_id = ids[id_idx];
                        std::size_t const rhs_id = ids[id_idx + 1];
                        std::size_t const res_id = ids[id_idx + 2];
                        // double const lhs_val = vals[val_idx];
                        // double const rhs_val = vals[val_idx + 1];
                        update_bivariate.template operator()<mul_type::ANY, mul_type::ANY>(
                          lhs_id, rhs_id, res_id /* , rhs_val, lhs_val */);
                    }
                    // val_idx += 2;
                    id_idx += 3;
                    break;
                }
                case OpCode::ADD_C: {
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        update_univariate.template operator()<mul_type::ONE>(arg_id, res_id);
                    }
                    id_idx += 2;
                    break;
                }
                case OpCode::SUB_C: {
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        update_univariate.template operator()<mul_type::MINUS_ONE>(arg_id, res_id);
                    }
                    id_idx += 2;
                    break;
                }
                case OpCode::MUL_C: {
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        // double const multiplier = vals[val_idx];
                        update_univariate.template operator()<mul_type::ANY>(arg_id, res_id /* , multiplier */);
                    }
                    // val_idx += 1;
                    id_idx += 2;
                    break;
                }
                case OpCode::NORM: {
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        // double const multiplier = 2.0 * vals[val_idx];
                        update_univariate.template operator()<mul_type::ANY>(arg_id, res_id /* , multiplier */);
                    }
                    // val_idx += 1;
                    id_idx += 2;
                    break;
                }
                case OpCode::INV: {
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        // double const multiplier = -vals[val_idx] * vals[val_idx];
                        update_univariate.template operator()<mul_type::ANY>(arg_id, res_id /* , multiplier */);
                    }
                    // val_idx += 1;
                    id_idx += 2;
                    break;
                }
                case OpCode::ABS: {
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        // double const multiplier = std::copysign(1.0, vals[val_idx]);
                        update_univariate.template operator()<mul_type::ANY>(arg_id, res_id /* , multiplier */);
                    }
                    // val_idx += 1;
                    id_idx += 2;
                    break;
                }
                case OpCode::EXP: {
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        // double const multiplier = vals[val_idx];
                        update_univariate.template operator()<mul_type::ANY>(arg_id, res_id /* , multiplier */);
                    }
                    // val_idx += 1;
                    id_idx += 2;
                    break;
                }
                case OpCode::LOG: {
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        // double const multiplier = 1.0 / vals[val_idx];
                        update_univariate.template operator()<mul_type::ANY>(arg_id, res_id /* , multiplier */);
                    }
                    // val_idx += 1;
                    id_idx += 2;
                    break;
                }
                case OpCode::ERF: {
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        // constexpr double two_over_root_pi = 2. * std::numbers::inv_sqrtpi_v<double>;
                        // double const multiplier = std::exp(-vals[val_idx] * vals[val_idx]) * two_over_root_pi;
                        update_univariate.template operator()<mul_type::ANY>(arg_id, res_id /* , multiplier */);
                    }
                    // val_idx += 1;
                    id_idx += 2;
                    break;
                }
                case OpCode::ERFC: {
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        // constexpr double minus_two_over_root_pi = -2. * std::numbers::inv_sqrtpi_v<double>;
                        // double const multiplier = std::exp(-vals[val_idx] * vals[val_idx]) * minus_two_over_root_pi;
                        update_univariate.template operator()<mul_type::ANY>(arg_id, res_id /* , multiplier */);
                    }
                    // val_idx += 1;
                    id_idx += 2;
                    break;
                }
                case OpCode::COS: {
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        // double const multiplier = -std::sin(vals[val_idx]);
                        update_univariate.template operator()<mul_type::ANY>(arg_id, res_id /* , multiplier */);
                    }
                    // val_idx += 2;
                    id_idx += 2;
                    break;
                }
                case OpCode::SQRT: {
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        // double const one_over_in = 1. / vals[val_idx];
                        // double const multiplier = 0.5 * vals[val_idx + 1] * one_over_in;
                        update_univariate.template operator()<mul_type::ANY>(arg_id, res_id /* , multiplier */);
                    }
                    // val_idx += 2;
                    id_idx += 2;
                    break;
                }
                case OpCode::POW_C: {
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        // double const lhs_arg = vals[val_idx];
                        // double const rhs_arg = vals[val_idx + 1];
                        // double const multiplier = rhs_arg != 0.0 ? rhs_arg * std::pow(lhs_arg, rhs_arg - 1.) : 0.0;
                        update_univariate.template operator()<mul_type::ANY>(arg_id, res_id /* , multiplier */);
                    }
                    // val_idx += 2;
                    id_idx += 2;
                    break;
                }
            }
        }

        // calculate actual buffer values to check
        const auto& vals = data.vals;
        std::size_t val_idx = pos.val_position;
        std::size_t mult_op_idx = 0;
        std::size_t mult_pos_idx = 0;
        std::size_t mult_bc_idx = 0;
        std::size_t mult_uc_idx = 0;

        buffer_multipliers.values.resize(buffer_multipliers.size);

        for (std::size_t op_idx = to; op_idx < from; ++op_idx) {
            OpCode const& op = ops[op_idx];
            bool const use_this_op = (number_dependents2[op_idx - to] > 0);

            switch (op) {
                case OpCode::REG_INPUT: {
                    break;
                }
                case OpCode::REG_OUTPUT: {
                    break;
                }
                case OpCode::ADD: {
                    if (use_this_op) {
                        // double constexpr lhs_der = 1.0;
                        // double constexpr rhs_der = 1.0;

                        // lhs multiplier
                        // {
                        //     double constexpr multiplier = lhs_der;
                        //     bool const op_m = op_multiply[mult_op_idx++];
                        //     std::size_t const pos = pos_multiplier[mult_pos_idx++];
                        //     if (op_m) {
                        //         values_type[pos] = combine_mul_type(values_type[pos], mul_type::ONE);
                        //     }
                        //     else {
                        //         buffer_multipliers.values[pos] = multiplier;
                        //         values_type[pos] = mul_type::ONE;
                        //     }
                        // }

                        // rhs multiplier
                        // {
                        //     double constexpr multiplier = rhs_der;
                        //     bool const op_m = op_multiply[mult_op_idx++];
                        //     std::size_t const pos = pos_multiplier[mult_pos_idx++];
                        //     if (op_m) {
                        //         values_type[pos] = combine_mul_type(values_type[pos], mul_type::ONE);
                        //     }
                        //     else {
                        //         buffer_multipliers.values[pos] = multiplier;
                        //         values_type[pos] = mul_type::ONE;
                        //     }
                        // }

                        if (bivariate_consolidate[mult_bc_idx++]) {
                            std::size_t const pos_lhs = pos_multiplier[mult_pos_idx++];
                            std::size_t const pos_rhs = pos_multiplier[mult_pos_idx++];
                            buffer_multipliers.values[pos_lhs] += buffer_multipliers.values[pos_rhs];
                            // values_type[pos_lhs] = mul_type::ANY;
                            if (univariate_consolidate[mult_uc_idx++]) {
                                std::size_t const pos = pos_multiplier[mult_pos_idx++];
                                buffer_multipliers.values[pos] *= buffer_multipliers.values[pos_lhs];
                                // values_type[pos] = mul_type::ANY;
                            }
                        }
                    }
                    break;
                }
                case OpCode::SUB: {
                    if (use_this_op) {
                        double constexpr lhs_der = 1.0;
                        double constexpr rhs_der = -1.0;

                        // lhs multiplier
                        {
                            double constexpr multiplier = lhs_der;
                            bool const op_m = op_multiply[mult_op_idx++];
                            std::size_t const pos = pos_multiplier[mult_pos_idx++];
                            if (op_m) {
                                // values_type[pos] = combine_mul_type(values_type[pos], mul_type::ONE);
                            }
                            else {
                                buffer_multipliers.values[pos] = multiplier;
                                // values_type[pos] = mul_type::ONE;
                            }
                        }

                        // rhs multiplier
                        {
                            double constexpr multiplier = rhs_der;
                            bool const op_m = op_multiply[mult_op_idx++];
                            std::size_t const pos = pos_multiplier[mult_pos_idx++];
                            if (op_m) {
                                buffer_multipliers.values[pos] = -buffer_multipliers.values[pos];
                                // values_type[pos] = combine_mul_type(values_type[pos], mul_type::MINUS_ONE);
                            }
                            else {
                                buffer_multipliers.values[pos] = multiplier;
                                // values_type[pos] = mul_type::MINUS_ONE;
                            }
                        }

                        if (bivariate_consolidate[mult_bc_idx++]) {
                            std::size_t const pos_lhs = pos_multiplier[mult_pos_idx++];
                            std::size_t const pos_rhs = pos_multiplier[mult_pos_idx++];
                            buffer_multipliers.values[pos_lhs] += buffer_multipliers.values[pos_rhs];
                            // values_type[pos_lhs] = mul_type::ANY;
                            if (univariate_consolidate[mult_uc_idx++]) {
                                std::size_t const pos = pos_multiplier[mult_pos_idx++];
                                buffer_multipliers.values[pos] *= buffer_multipliers.values[pos_lhs];
                                // values_type[pos] = mul_type::ANY;
                            }
                        }
                    }
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
                                // values_type[pos] = mul_type::ANY;
                            }
                            else {
                                buffer_multipliers.values[pos] = multiplier;
                                // values_type[pos] = mul_type::ANY;
                            }
                        }

                        // rhs multiplier
                        {
                            double const multiplier = lhs_val;
                            bool const op_m = op_multiply[mult_op_idx++];
                            std::size_t const pos = pos_multiplier[mult_pos_idx++];
                            if (op_m) {
                                buffer_multipliers.values[pos] *= multiplier;
                                // values_type[pos] = mul_type::ANY;
                            }
                            else {
                                buffer_multipliers.values[pos] = multiplier;
                                // values_type[pos] = mul_type::ANY;
                            }
                        }

                        if (bivariate_consolidate[mult_bc_idx++]) {
                            std::size_t const pos_lhs = pos_multiplier[mult_pos_idx++];
                            std::size_t const pos_rhs = pos_multiplier[mult_pos_idx++];
                            buffer_multipliers.values[pos_lhs] += buffer_multipliers.values[pos_rhs];
                            // values_type[pos_lhs] = mul_type::ANY;
                            if (univariate_consolidate[mult_uc_idx++]) {
                                std::size_t const pos = pos_multiplier[mult_pos_idx++];
                                buffer_multipliers.values[pos] *= buffer_multipliers.values[pos_lhs];
                                // values_type[pos] = mul_type::ANY;
                            }
                        }
                    }
                    val_idx += 2;
                    break;
                }
                case OpCode::ADD_C: {
                    if (use_this_op) {
                        double constexpr multiplier = 1.;
                        bool const op_m = op_multiply[mult_op_idx++];
                        std::size_t const pos = pos_multiplier[mult_pos_idx++];
                        if (op_m) {
                            // values_type[pos] = combine_mul_type(values_type[pos], mul_type::ONE);
                        }
                        else {
                            buffer_multipliers.values[pos] = multiplier;
                            // values_type[pos] = mul_type::ONE;
                        }
                    }
                    break;
                }
                case OpCode::SUB_C: {
                    if (use_this_op) {
                        double constexpr multiplier = -1.;
                        bool const op_m = op_multiply[mult_op_idx++];
                        std::size_t const pos = pos_multiplier[mult_pos_idx++];
                        if (op_m) {
                            buffer_multipliers.values[pos] = -buffer_multipliers.values[pos];
                            // values_type[pos] = combine_mul_type(values_type[pos], mul_type::MINUS_ONE);
                        }
                        else {
                            buffer_multipliers.values[pos] = multiplier;
                            // values_type[pos] = mul_type::MINUS_ONE;
                        }
                    }
                    break;
                }
                case OpCode::MUL_C: {
                    if (use_this_op) {
                        double const multiplier = vals[val_idx];
                        bool const op_m = op_multiply[mult_op_idx++];
                        std::size_t const pos = pos_multiplier[mult_pos_idx++];
                        if (op_m) {
                            buffer_multipliers.values[pos] *= multiplier;
                            // values_type[pos] = mul_type::ANY;
                        }
                        else {
                            buffer_multipliers.values[pos] = multiplier;
                            // values_type[pos] = mul_type::ANY;
                        }
                    }
                    val_idx += 1;
                    break;
                }
                case OpCode::NORM: {
                    if (use_this_op) {
                        double const multiplier = 2.0 * vals[val_idx];
                        bool const op_m = op_multiply[mult_op_idx++];
                        std::size_t const pos = pos_multiplier[mult_pos_idx++];
                        if (op_m) {
                            buffer_multipliers.values[pos] *= multiplier;
                            // values_type[pos] = mul_type::ANY;
                        }
                        else {
                            buffer_multipliers.values[pos] = multiplier;
                            // values_type[pos] = mul_type::ANY;
                        }
                    }
                    val_idx += 1;
                    break;
                }
                case OpCode::INV: {
                    if (use_this_op) {
                        double const multiplier = -vals[val_idx] * vals[val_idx];
                        bool const op_m = op_multiply[mult_op_idx++];
                        std::size_t const pos = pos_multiplier[mult_pos_idx++];
                        if (op_m) {
                            buffer_multipliers.values[pos] *= multiplier;
                            // values_type[pos] = mul_type::ANY;
                        }
                        else {
                            buffer_multipliers.values[pos] = multiplier;
                            // values_type[pos] = mul_type::ANY;
                        }
                    }
                    val_idx += 1;
                    break;
                }
                case OpCode::ABS: {
                    if (use_this_op) {
                        double const multiplier = std::copysign(1.0, vals[val_idx]);
                        bool const op_m = op_multiply[mult_op_idx++];
                        std::size_t const pos = pos_multiplier[mult_pos_idx++];
                        if (op_m) {
                            buffer_multipliers.values[pos] *= multiplier;
                            // values_type[pos] = mul_type::ANY;
                        }
                        else {
                            buffer_multipliers.values[pos] = multiplier;
                            // values_type[pos] = mul_type::ANY;
                        }
                    }
                    val_idx += 1;
                    break;
                }
                case OpCode::EXP: {
                    if (use_this_op) {
                        double const multiplier = vals[val_idx];
                        bool const op_m = op_multiply[mult_op_idx++];
                        std::size_t const pos = pos_multiplier[mult_pos_idx++];
                        if (op_m) {
                            buffer_multipliers.values[pos] *= multiplier;
                            // values_type[pos] = mul_type::ANY;
                        }
                        else {
                            buffer_multipliers.values[pos] = multiplier;
                            // values_type[pos] = mul_type::ANY;
                        }
                    }
                    val_idx += 1;
                    break;
                }
                case OpCode::LOG: {
                    if (use_this_op) {
                        double const multiplier = 1.0 / vals[val_idx];
                        bool const op_m = op_multiply[mult_op_idx++];
                        std::size_t const pos = pos_multiplier[mult_pos_idx++];
                        if (op_m) {
                            buffer_multipliers.values[pos] *= multiplier;
                            // values_type[pos] = mul_type::ANY;
                        }
                        else {
                            buffer_multipliers.values[pos] = multiplier;
                            // values_type[pos] = mul_type::ANY;
                        }
                    }
                    val_idx += 1;
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
                            // values_type[pos] = mul_type::ANY;
                        }
                        else {
                            buffer_multipliers.values[pos] = multiplier;
                            // values_type[pos] = mul_type::ANY;
                        }
                    }
                    val_idx += 1;
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
                            // values_type[pos] = mul_type::ANY;
                        }
                        else {
                            buffer_multipliers.values[pos] = multiplier;
                            // values_type[pos] = mul_type::ANY;
                        }
                    }
                    val_idx += 1;
                    break;
                }
                case OpCode::COS: {
                    if (use_this_op) {
                        double const multiplier = -std::sin(vals[val_idx]);
                        bool const op_m = op_multiply[mult_op_idx++];
                        std::size_t const pos = pos_multiplier[mult_pos_idx++];
                        if (op_m) {
                            buffer_multipliers.values[pos] *= multiplier;
                            // values_type[pos] = mul_type::ANY;
                        }
                        else {
                            buffer_multipliers.values[pos] = multiplier;
                            // values_type[pos] = mul_type::ANY;
                        }
                    }
                    val_idx += 2;
                    break;
                }
                case OpCode::SQRT: {
                    if (use_this_op) {
                        double const one_over_in = 1. / vals[val_idx];
                        double const multiplier = 0.5 * vals[val_idx + 1] * one_over_in;
                        bool const op_m = op_multiply[mult_op_idx++];
                        std::size_t const pos = pos_multiplier[mult_pos_idx++];
                        if (op_m) {
                            buffer_multipliers.values[pos] *= multiplier;
                            // values_type[pos] = mul_type::ANY;
                        }
                        else {
                            buffer_multipliers.values[pos] = multiplier;
                            // values_type[pos] = mul_type::ANY;
                        }
                    }
                    val_idx += 2;
                    break;
                }
                case OpCode::POW_C: {
                    if (use_this_op) {
                        double const lhs_arg = vals[val_idx];
                        double const rhs_arg = vals[val_idx + 1];
                        double const multiplier = rhs_arg != 0.0 ? rhs_arg * std::pow(lhs_arg, rhs_arg - 1.) : 0.0;
                        bool const op_m = op_multiply[mult_op_idx++];
                        std::size_t const pos = pos_multiplier[mult_pos_idx++];
                        if (op_m) {
                            buffer_multipliers.values[pos] *= multiplier;
                            // values_type[pos] = mul_type::ANY;
                        }
                        else {
                            buffer_multipliers.values[pos] = multiplier;
                            // values_type[pos] = mul_type::ANY;
                        }
                    }
                    val_idx += 2;
                    break;
                }
            }
        }

        id_idx = ids.size();

        auto& buffer_free_positions = buffers.back().free_positions;

        auto& this_on_which_buffer = result.on_which_buffer;
        auto& this_pos = result.pos;
        auto& this_lossy_op = result.lossy_op;

        auto copy = [&](std::size_t const out_pos, std::size_t const in_pos) {
            this_lossy_op.push_back(LossyOpCode::COPY);
            this_pos.push_back(out_pos);
            this_pos.push_back(in_pos);
        };

        auto copy_minus = [&](std::size_t const out_pos, std::size_t const in_pos) {
            this_lossy_op.push_back(LossyOpCode::COPY_MINUS);
            this_pos.push_back(out_pos);
            this_pos.push_back(in_pos);
        };

        auto add = [&](std::size_t const out_pos, std::size_t const in_pos, std::uint8_t const which) {
            this_lossy_op.push_back(LossyOpCode::ADD);
            this_pos.push_back(out_pos);
            this_pos.push_back(in_pos);
            this_on_which_buffer.push_back(which);
        };

        auto sub = [&](std::size_t const out_pos, std::size_t const in_pos, std::uint8_t const which) {
            this_lossy_op.push_back(LossyOpCode::SUB);
            this_pos.push_back(out_pos);
            this_pos.push_back(in_pos);
            this_on_which_buffer.push_back(which);
        };

        auto minus_inplace = [&](std::size_t const pos) {
            this_lossy_op.push_back(LossyOpCode::MINUS_INPLACE);
            this_pos.push_back(pos);
        };

        auto mul_inplace = [&](std::size_t const pos) {
            this_lossy_op.push_back(LossyOpCode::MUL_INPLACE);
            this_pos.push_back(pos);
        };

        auto mul_add = [&](std::size_t const out_pos, std::size_t const in_pos, std::uint8_t const which) {
            this_lossy_op.push_back(LossyOpCode::MUL_ADD);
            this_pos.push_back(out_pos);
            this_pos.push_back(in_pos);
            this_on_which_buffer.push_back(which);
        };

        auto mul_set = [&](std::size_t const out_pos, std::size_t const in_pos) {
            this_lossy_op.push_back(LossyOpCode::MUL_SET);
            this_pos.push_back(out_pos);
            this_pos.push_back(in_pos);
        };

        // output is: is_new, is_current, buffer_id, position
        auto get_loc = [&](std::size_t id) -> std::tuple<bool, std::uint8_t> {
            auto it = std::upper_bound(checkpoints.begin(), checkpoints.end(), id);
            auto buffer_id = static_cast<std::uint8_t>(std::distance(checkpoints.cbegin(), it) - 1);
            return { it == checkpoints.end(), buffer_id };
        };

        auto update_loc = [&](std::size_t& arg_pos, std::uint8_t buffer_id) {
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

        auto reset_loc = [&](std::size_t arg_pos, std::uint8_t buffer_id) {
            if (buffers[buffer_id].values.size() > (arg_pos * num_lanes)) {
                if constexpr (Vectorised) {
                    double* dest = &buffers[buffer_id].values[arg_pos * num_lanes];
#pragma omp simd
                    for (std::size_t i = 0; i < num_lanes; ++i) {
                        dest[i] = 0.;
                    }
                }
                else {
                    buffers[buffer_id].values[arg_pos] = 0.;
                }
            }
        };

        for (std::size_t op_idx = from; op_idx-- > to;) {
            OpCode const& op = ops[op_idx];
            bool const use_this_op = node_location_on_buffer[op_idx] != passive_id<std::size_t>;

            switch (op) {
                case OpCode::REG_INPUT: {
                    id_idx -= 1;
                    if constexpr (Reset) {
                        std::size_t const id = ids[id_idx];
                        std::size_t& pos = node_location_on_buffer[id];
                        buffer_free_positions.push_back(pos);
                        pos = passive_id<std::size_t>;
                    }
                    break;
                }
                case OpCode::REG_OUTPUT: {
                    id_idx -= 2;
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];

                        std::size_t& res_pos = node_location_on_buffer[res_id];
                        std::size_t& arg_pos = node_location_on_buffer[arg_id];

                        auto const arg_pos_data = get_loc(arg_id);
                        bool const arg_is_new = (arg_pos == passive_id<std::size_t>);
                        bool const arg_is_in_current_buffer = std::get<0>(arg_pos_data);
                        bool const arg_inplace = arg_is_new && arg_is_in_current_buffer;

                        if constexpr (Reset) {
                            if (arg_inplace) {
                                // res id should now be lhs id, avoiding a copy and a potential buffer increase
                                std::swap(arg_pos, res_pos);
                            }
                            else {
                                update_loc(arg_pos, std::get<1>(arg_pos_data));
                                add(res_pos, arg_pos, std::get<1>(arg_pos_data));
                                buffer_free_positions.push_back(res_pos);
                                res_pos = passive_id<std::size_t>;
                                if (!arg_is_in_current_buffer && arg_is_new) {
                                    // if the argument is not in the current buffer and it's not new, we need to reset
                                    // its value to 0
                                    reset_loc(arg_pos, std::get<1>(arg_pos_data));
                                }
                            }
                        }
                        else {
                            update_loc(arg_pos, std::get<1>(arg_pos_data));

                            // we don't reset so we don't free res_pos,
                            // it will be potentially used in the next operations as an active node location
                            if (arg_inplace) {
                                copy(res_pos, arg_pos);
                            }
                            else {
                                add(res_pos, arg_pos, std::get<1>(arg_pos_data));
                                if (!arg_is_in_current_buffer && arg_is_new) {
                                    // if the argument is not in the current buffer and it's not new, we need to reset
                                    // its value to 0
                                    reset_loc(arg_pos, std::get<1>(arg_pos_data));
                                }
                            }
                        }
                    }
                    break;
                }
                case OpCode::ADD: {
                    id_idx -= 3;
                    if (use_this_op) {
                        std::size_t const lhs_id = ids[id_idx];
                        std::size_t const rhs_id = ids[id_idx + 1];
                        std::size_t const res_id = ids[id_idx + 2];

                        std::size_t& res_pos = node_location_on_buffer[res_id];
                        std::size_t& lhs_pos = node_location_on_buffer[lhs_id];
                        std::size_t& rhs_pos = node_location_on_buffer[rhs_id];

                        auto const lhs_pos_data = get_loc(lhs_id);
                        auto const rhs_pos_data = get_loc(rhs_id);

                        bool const lhs_is_new = (lhs_pos == passive_id<std::size_t>);
                        bool const rhs_is_new = (rhs_pos == passive_id<std::size_t>);
                        bool const lhs_is_in_current_buffer = std::get<0>(lhs_pos_data);
                        bool const rhs_is_in_current_buffer = std::get<0>(rhs_pos_data);
                        bool const lhs_inplace = lhs_is_new && lhs_is_in_current_buffer;
                        bool const rhs_inplace_pre = rhs_is_new && rhs_is_in_current_buffer;
                        bool const rhs_inplace = !lhs_inplace && rhs_inplace_pre;

                        if (!lhs_inplace) {
                            update_loc(lhs_pos, std::get<1>(lhs_pos_data));
                            add(res_pos, lhs_pos, std::get<1>(lhs_pos_data));
                            if (!lhs_is_in_current_buffer && lhs_is_new) {
                                // if the argument is not in the current buffer and it's not new, we need to reset its
                                // value to 0
                                reset_loc(lhs_pos, std::get<1>(lhs_pos_data));
                            }
                        }

                        if (!rhs_inplace) {
                            update_loc(rhs_pos, std::get<1>(rhs_pos_data));
                            if (rhs_inplace_pre) {
                                copy(res_pos, rhs_pos);
                            }
                            else {
                                add(res_pos, rhs_pos, std::get<1>(rhs_pos_data));
                                if (!rhs_is_in_current_buffer && rhs_is_new) {
                                    // if the argument is not in the current buffer and it's not new, we need to reset
                                    // its value to 0
                                    reset_loc(rhs_pos, std::get<1>(rhs_pos_data));
                                }
                            }
                        }

                        if (lhs_inplace) {
                            // res id should now be lhs id, avoiding a copy and a potential buffer increase
                            std::swap(lhs_pos, res_pos);
                        }
                        else if (rhs_inplace) {
                            // res id should now be rhs id, avoiding a copy and a potential buffer increase
                            std::swap(rhs_pos, res_pos);
                        }
                        else {
                            // don't forget to free res_id from the buffer!
                            buffer_free_positions.push_back(res_pos);
                            res_pos = passive_id<std::size_t>;
                        }
                    }
                    break;
                }
                case OpCode::SUB: {
                    id_idx -= 3;
                    if (use_this_op) {
                        std::size_t const lhs_id = ids[id_idx];
                        std::size_t const rhs_id = ids[id_idx + 1];
                        std::size_t const res_id = ids[id_idx + 2];

                        std::size_t& res_pos = node_location_on_buffer[res_id];
                        std::size_t& lhs_pos = node_location_on_buffer[lhs_id];
                        std::size_t& rhs_pos = node_location_on_buffer[rhs_id];

                        auto const lhs_pos_data = get_loc(lhs_id);
                        auto const rhs_pos_data = get_loc(rhs_id);

                        bool const lhs_is_new = (lhs_pos == passive_id<std::size_t>);
                        bool const rhs_is_new = (rhs_pos == passive_id<std::size_t>);
                        bool const lhs_is_in_current_buffer = std::get<0>(lhs_pos_data);
                        bool const rhs_is_in_current_buffer = std::get<0>(rhs_pos_data);
                        bool const lhs_inplace = lhs_is_new && lhs_is_in_current_buffer;
                        bool const rhs_inplace_pre = rhs_is_new && rhs_is_in_current_buffer;
                        bool const rhs_inplace = !lhs_inplace && rhs_inplace_pre;

                        if (!lhs_inplace) {
                            update_loc(lhs_pos, std::get<1>(lhs_pos_data));
                            add(res_pos, lhs_pos, std::get<1>(lhs_pos_data));
                            if (!lhs_is_in_current_buffer && lhs_is_new) {
                                // if the argument is not in the current buffer and it's not new, we need to reset its
                                // value to 0
                                reset_loc(lhs_pos, std::get<1>(lhs_pos_data));
                            }
                        }

                        if (!rhs_inplace) {
                            update_loc(rhs_pos, std::get<1>(rhs_pos_data));
                            if (rhs_inplace_pre) {
                                copy_minus(res_pos, rhs_pos);
                            }
                            else {
                                sub(res_pos, rhs_pos, std::get<1>(rhs_pos_data));
                                if (!rhs_is_in_current_buffer && rhs_is_new) {
                                    // if the argument is not in the current buffer and it's not new, we need to reset
                                    // its value to 0
                                    reset_loc(rhs_pos, std::get<1>(rhs_pos_data));
                                }
                            }
                        }

                        if (lhs_inplace) {
                            // res id should now be lhs id, avoiding a copy and a potential buffer increase
                            std::swap(lhs_pos, res_pos);
                        }
                        else if (rhs_inplace) {
                            // this is a subtraction, so we need to negate the value in the buffer
                            minus_inplace(res_pos);
                            // res id should now be rhs id, avoiding a copy and a potential buffer increase
                            std::swap(rhs_pos, res_pos);
                        }
                        else {
                            // don't forget to free res_id from the buffer!
                            buffer_free_positions.push_back(res_pos);
                            res_pos = passive_id<std::size_t>;
                        }
                    }
                    break;
                }
                case OpCode::MUL: {
                    id_idx -= 3;
                    result.use_op.push_back(use_this_op);
                    if (use_this_op) {
                        std::size_t const lhs_id = ids[id_idx];
                        std::size_t const rhs_id = ids[id_idx + 1];
                        std::size_t const res_id = ids[id_idx + 2];

                        std::size_t& res_pos = node_location_on_buffer[res_id];
                        std::size_t& lhs_pos = node_location_on_buffer[lhs_id];
                        std::size_t& rhs_pos = node_location_on_buffer[rhs_id];

                        auto const lhs_pos_data = get_loc(lhs_id);
                        auto const rhs_pos_data = get_loc(rhs_id);

                        bool const lhs_is_new = (lhs_pos == passive_id<std::size_t>);
                        bool const rhs_is_new = (rhs_pos == passive_id<std::size_t>);
                        bool const lhs_is_in_current_buffer = std::get<0>(lhs_pos_data);
                        bool const rhs_is_in_current_buffer = std::get<0>(rhs_pos_data);
                        bool const lhs_inplace = lhs_is_new && lhs_is_in_current_buffer;
                        bool const rhs_inplace_pre = rhs_is_new && rhs_is_in_current_buffer;
                        bool const rhs_inplace = !lhs_inplace && rhs_inplace_pre;

                        // this is the only case when rhs_multiplier is used BEFORE lhs_multiplier.
                        result.invert_mult.push_back(lhs_inplace);

                        if (!lhs_inplace) {
                            update_loc(lhs_pos, std::get<1>(lhs_pos_data));
                            mul_add(res_pos, lhs_pos, std::get<1>(lhs_pos_data));
                            if (!lhs_is_in_current_buffer && lhs_is_new) {
                                // if the argument is not in the current buffer and it's not new, we need to reset its
                                // value to 0
                                reset_loc(lhs_pos, std::get<1>(lhs_pos_data));
                            }
                        }

                        if (!rhs_inplace) {
                            update_loc(rhs_pos, std::get<1>(rhs_pos_data));
                            if (rhs_inplace_pre) {
                                mul_set(res_pos, rhs_pos);
                            }
                            else {
                                mul_add(res_pos, rhs_pos, std::get<1>(rhs_pos_data));
                                if (!rhs_is_in_current_buffer && rhs_is_new) {
                                    // if the argument is not in the current buffer and it's not new, we need to reset
                                    // its value to 0
                                    reset_loc(rhs_pos, std::get<1>(rhs_pos_data));
                                }
                            }
                        }

                        if (lhs_inplace) {
                            mul_inplace(res_pos);
                            // res id should now be lhs id, avoiding a copy and a potential buffer increase
                            std::swap(lhs_pos, res_pos);
                        }
                        else if (rhs_inplace) {
                            mul_inplace(res_pos);
                            // res id should now be rhs id, avoiding a copy and a potential buffer increase
                            std::swap(rhs_pos, res_pos);
                        }
                        else {
                            // don't forget to free res_id from the buffer!
                            buffer_free_positions.push_back(res_pos);
                            res_pos = passive_id<std::size_t>;
                        }
                    }
                    break;
                }
                case OpCode::ADD_C: {
                    id_idx -= 2;
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];

                        std::size_t& res_pos = node_location_on_buffer[res_id];
                        std::size_t& arg_pos = node_location_on_buffer[arg_id];

                        auto const arg_pos_data = get_loc(arg_id);
                        bool const arg_is_new = (arg_pos == passive_id<std::size_t>);
                        bool const arg_is_in_current_buffer = std::get<0>(arg_pos_data);
                        bool const arg_inplace = arg_is_new && arg_is_in_current_buffer;

                        if (arg_inplace) {
                            // res id should now be lhs id, avoiding a copy and a potential buffer increase
                            std::swap(arg_pos, res_pos);
                        }
                        else {
                            update_loc(arg_pos, std::get<1>(arg_pos_data));
                            add(res_pos, arg_pos, std::get<1>(arg_pos_data));
                            buffer_free_positions.push_back(res_pos);
                            res_pos = passive_id<std::size_t>;
                            if (!arg_is_in_current_buffer && arg_is_new) {
                                // if the argument is not in the current buffer and it's not new, we need to reset its
                                // value to 0
                                reset_loc(arg_pos, std::get<1>(arg_pos_data));
                            }
                        }
                    }
                    break;
                }
                case OpCode::SUB_C: {
                    id_idx -= 2;
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];

                        std::size_t& res_pos = node_location_on_buffer[res_id];
                        std::size_t& arg_pos = node_location_on_buffer[arg_id];

                        auto const arg_pos_data = get_loc(arg_id);
                        bool const arg_is_new = (arg_pos == passive_id<std::size_t>);
                        bool const arg_is_in_current_buffer = std::get<0>(arg_pos_data);
                        bool const arg_inplace = arg_is_new && arg_is_in_current_buffer;

                        if (arg_inplace) {
                            // this is a subtraction, so we need to negate the value in the buffer
                            minus_inplace(res_pos);
                            // res id should now be lhs id, avoiding a copy and a potential buffer increase
                            std::swap(arg_pos, res_pos);
                        }
                        else {
                            update_loc(arg_pos, std::get<1>(arg_pos_data));
                            sub(res_pos, arg_pos, std::get<1>(arg_pos_data));
                            buffer_free_positions.push_back(res_pos);
                            res_pos = passive_id<std::size_t>;
                            if (!arg_is_in_current_buffer && arg_is_new) {
                                // if the argument is not in the current buffer and it's not new, we need to reset its
                                // value to 0
                                reset_loc(arg_pos, std::get<1>(arg_pos_data));
                            }
                        }
                    }
                    break;
                }
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
                    result.use_op.push_back(use_this_op);
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];

                        std::size_t& res_pos = node_location_on_buffer[res_id];
                        std::size_t& arg_pos = node_location_on_buffer[arg_id];

                        auto const arg_pos_data = get_loc(arg_id);
                        bool const arg_is_new = (arg_pos == passive_id<std::size_t>);
                        bool const arg_is_in_current_buffer = std::get<0>(arg_pos_data);
                        bool const arg_inplace = arg_is_new && arg_is_in_current_buffer;

                        if (arg_inplace) {
                            mul_inplace(res_pos);
                            // res id should now be arg id, avoiding a copy and a potential buffer increase
                            std::swap(arg_pos, res_pos);
                        }
                        else {
                            update_loc(arg_pos, std::get<1>(arg_pos_data));
                            mul_add(res_pos, arg_pos, std::get<1>(arg_pos_data));
                            buffer_free_positions.push_back(res_pos);
                            res_pos = passive_id<std::size_t>;
                            if (!arg_is_in_current_buffer && arg_is_new) {
                                // if the argument is not in the current buffer and it's not new, we need to reset its
                                // value to 0
                                reset_loc(arg_pos, std::get<1>(arg_pos_data));
                            }
                        }
                    }
                    break;
                }
            }
        }
        result.buffer_size = buffers.back().size;

        return result;
    };

    class ValueFetcher {
      public:
        explicit ValueFetcher(TapeData const& data_in, LossyTape const& lossy_tape_in, std::size_t to_in)
          : data(data_in)
          , lossy_tape(lossy_tape_in)
          , op_idx(data.next_id)
          , to(to_in)
          , val_idx(data.vals.size())
        {
        }

        auto operator()() -> double
        {
            if (use_cache) {
                use_cache = false;
                return cached_value;
            }

            const auto& ops = data.ops;

            for (; op_idx-- > to;) {
                OpCode const& op = ops[op_idx];
                // bool const use_this_op = this->node_location_on_buffer[op_idx] != passive_id<std::size_t>;

                switch (op) {
                    case OpCode::REG_INPUT:
                    case OpCode::REG_OUTPUT:
                    case OpCode::ADD:
                    case OpCode::SUB:
                    case OpCode::ADD_C:
                    case OpCode::SUB_C: {
                        break;
                    }
                    case OpCode::MUL: {
                        val_idx -= 2;
                        bool const use_this_op = lossy_tape.use_op[use_op_idx];
                        use_op_idx += 1;
                        if (use_this_op) {

                            double const lhs_val = data.vals[val_idx];
                            double const rhs_val = data.vals[val_idx + 1];
                            bool const this_invert_idx = lossy_tape.invert_mult[invert_mult_idx];
                            invert_mult_idx += 1;

                            use_cache = true;
                            if (this_invert_idx) {
                                cached_value = rhs_val;
                                return lhs_val;
                            }
                            else {
                                cached_value = lhs_val;
                                return rhs_val;
                            }
                        }
                        break;
                    }
                    case OpCode::MUL_C: {
                        val_idx -= 1;
                        bool const use_this_op = lossy_tape.use_op[use_op_idx];
                        use_op_idx += 1;
                        if (use_this_op) {
                            double const der_local_1 = data.vals[val_idx];
                            return der_local_1;
                        }
                        break;
                    }
                    case OpCode::NORM: {
                        val_idx -= 1;
                        bool const use_this_op = lossy_tape.use_op[use_op_idx];
                        use_op_idx += 1;
                        if (use_this_op) {
                            double const der_local_1 = 2.0 * data.vals[val_idx];
                            return der_local_1;
                        }
                        break;
                    }
                    case OpCode::INV: {
                        val_idx -= 1;
                        bool const use_this_op = lossy_tape.use_op[use_op_idx];
                        use_op_idx += 1;
                        if (use_this_op) {
                            double const der_local_1 = -data.vals[val_idx] * data.vals[val_idx];
                            return der_local_1;
                        }
                        break;
                    }
                    case OpCode::ABS: {
                        val_idx -= 1;
                        bool const use_this_op = lossy_tape.use_op[use_op_idx];
                        use_op_idx += 1;
                        if (use_this_op) {
                            double const der_local_1 = std::copysign(1.0, data.vals[val_idx]);
                            return der_local_1;
                        }
                        break;
                    }
                    case OpCode::EXP: {
                        val_idx -= 1;
                        bool const use_this_op = lossy_tape.use_op[use_op_idx];
                        use_op_idx += 1;
                        if (use_this_op) {
                            double const der_local_1 = data.vals[val_idx];
                            return der_local_1;
                        }
                        break;
                    }
                    case OpCode::LOG: {
                        val_idx -= 1;
                        bool const use_this_op = lossy_tape.use_op[use_op_idx];
                        use_op_idx += 1;
                        if (use_this_op) {
                            double const der_local_1 = 1.0 / data.vals[val_idx];
                            return der_local_1;
                        }
                        break;
                    }
                    case OpCode::ERF: {
                        val_idx -= 1;
                        bool const use_this_op = lossy_tape.use_op[use_op_idx];
                        use_op_idx += 1;
                        if (use_this_op) {
                            constexpr double two_over_root_pi = 2. * std::numbers::inv_sqrtpi_v<double>;
                            double const der_local_1 =
                              std::exp(-data.vals[val_idx] * data.vals[val_idx]) * two_over_root_pi;
                            return der_local_1;
                        }
                        break;
                    }
                    case OpCode::ERFC: {
                        val_idx -= 1;
                        bool const use_this_op = lossy_tape.use_op[use_op_idx];
                        use_op_idx += 1;
                        if (use_this_op) {
                            constexpr double minus_two_over_root_pi = -2. * std::numbers::inv_sqrtpi_v<double>;
                            double const der_local_1 =
                              std::exp(-data.vals[val_idx] * data.vals[val_idx]) * minus_two_over_root_pi;
                            return der_local_1;
                        }
                        break;
                    }
                    case OpCode::COS: {
                        val_idx -= 2;
                        bool const use_this_op = lossy_tape.use_op[use_op_idx];
                        use_op_idx += 1;
                        if (use_this_op) {
                            double const der_local_1 = -std::sin(data.vals[val_idx]);
                            return der_local_1;
                        }
                        break;
                    }
                    case OpCode::SQRT: {
                        val_idx -= 2;
                        bool const use_this_op = lossy_tape.use_op[use_op_idx];
                        use_op_idx += 1;
                        if (use_this_op) {
                            double const one_over_in = 1. / data.vals[val_idx];
                            double const der_local_1 = 0.5 * data.vals[val_idx + 1] * one_over_in;
                            return der_local_1;
                        }
                        break;
                    }
                    case OpCode::POW_C: {
                        val_idx -= 2;
                        bool const use_this_op = lossy_tape.use_op[use_op_idx];
                        use_op_idx += 1;
                        if (use_this_op) {
                            double const lhs_arg = data.vals[val_idx];
                            double const rhs_arg = data.vals[val_idx + 1];
                            double const der_local_1 = rhs_arg != 0.0 ? rhs_arg * std::pow(lhs_arg, rhs_arg - 1.) : 0.0;
                            return der_local_1;
                        }
                        break;
                    }
                }
            }

            // should not reach here
            throw;
            return 0.;
        }

      private:
        TapeData const& data;
        LossyTape const& lossy_tape;

        std::size_t op_idx{ 0 };
        std::size_t to{ 0 };

        std::size_t val_idx{ 0 };
        std::size_t invert_mult_idx{ 0 };
        std::size_t use_op_idx{ 0 };

        bool use_cache{ false };
        double cached_value{ 0. };
    };

    std::size_t to = pos.op_position;
    std::size_t from = data.next_id;

    const auto& ops = data.ops;
    const auto& vals = data.vals;
    const auto& ids = data.ids;

    std::size_t val_idx = vals.size();
    std::size_t id_idx = ids.size();

    this->node_location_on_buffer.resize(ops.size(), passive_id<std::size_t>);

    if (from == this->checkpoints.back()) {
        this->checkpoints.pop_back();
        this->buffers.pop_back();
        this->hash_to_lossy_tape.pop_back();
    }

    if (this->checkpoints.size() > 1 && to != this->checkpoints.back()) {
        // this is now allowed. if we have a buffer after the first one,
        // it necessarily has to be in MC context, so backpropagation
        // has to be all the way to the next checkpoint
        // why? because otherwise code becomes too complex
        throw;
    }

    // we assume the first section is not part of a MC simulation,
    // so the first section lossy opcode won't be stored.
    // if (this->buffers.size() > 1) {
    if (true) {
        this->hash_to_lossy_tape.push_back(std::map<std::size_t, LossyTape>{});

        auto const h = hash(pos, data);

        bool const tape_exists = hash_to_lossy_tape.back().find(h) != hash_to_lossy_tape.back().end();
        auto& lossy_tape = hash_to_lossy_tape.back()[h];
        if (!tape_exists) {
            // we need to compute the lossy tape for this section.
            lossy_tape =
              convert_to_lossy_tape(pos, data, node_location_on_buffer, checkpoints, buffers, this->m_num_lanes);
        }
        else {
            buffers.back().size = lossy_tape.buffer_size;

            // check equality of lossy tapes
            // auto new_lossy_tape = convert_to_lossy_tape<Reset>(to, data, node_location_on_buffer, checkpoints,
            // buffers);

            // auto equal_lossy_tape = [](LossyTape const& lt1, LossyTape const& lt2) -> bool {
            //     if (lt1.lossy_op != lt2.lossy_op) {
            //         return false;
            //     }

            //     if (lt1.on_which_buffer != lt2.on_which_buffer) {
            //         return false;
            //     }

            //     if (lt1.pos != lt2.pos) {
            //         return false;
            //     }

            //     return true;
            // };

            // if (!equal_lossy_tape(lossy_tape, new_lossy_tape)) {
            //     throw;
            // }
        }

        ValueFetcher value_fetcher(data, lossy_tape, to);

        auto const& lossy_op = lossy_tape.lossy_op;
        auto const& on_which_buffer = lossy_tape.on_which_buffer;
        auto const& pos = lossy_tape.pos;
        auto& buffer_vals = buffers.back().values;

        for (auto& b : buffers) {
            b.values.resize(b.size * this->m_num_lanes, 0.);
        }

        std::size_t id_idx_l = 0;
        std::size_t on_which_buffer_idx = 0;

        for (std::size_t lossy_op_idx = 0; lossy_op_idx < lossy_op.size(); ++lossy_op_idx) {
            auto const& op = lossy_op[lossy_op_idx];

            switch (op) {
                case LossyOpCode::COPY: {
                    std::size_t const out_pos = pos[id_idx_l++];
                    std::size_t const in_pos = pos[id_idx_l++];
                    if constexpr (Vectorised) {
                        const double* src = &buffer_vals[out_pos * this->m_num_lanes];
                        double* dest = &buffer_vals[in_pos * this->m_num_lanes];
#pragma omp simd
                        for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                            dest[i] = src[i];
                        }
                    }
                    else {
                        buffer_vals[in_pos] = buffer_vals[out_pos];
                    }
                    break;
                }
                case LossyOpCode::COPY_MINUS: {
                    std::size_t const out_pos = pos[id_idx_l++];
                    std::size_t const in_pos = pos[id_idx_l++];
                    if constexpr (Vectorised) {
                        const double* src = &buffer_vals[out_pos * this->m_num_lanes];
                        double* dest = &buffer_vals[in_pos * this->m_num_lanes];
#pragma omp simd
                        for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                            dest[i] = -src[i];
                        }
                    }
                    else {
                        buffer_vals[in_pos] = -buffer_vals[out_pos];
                    }
                    break;
                }
                case LossyOpCode::ADD: {
                    std::uint8_t const which = on_which_buffer[on_which_buffer_idx++];
                    std::size_t const out_pos = pos[id_idx_l++];
                    std::size_t const in_pos = pos[id_idx_l++];
                    if constexpr (Vectorised) {
                        const double* src = &buffer_vals[out_pos * this->m_num_lanes];
                        double* dest = &buffers[which].values[in_pos * this->m_num_lanes];
#pragma omp simd
                        for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                            dest[i] += src[i];
                        }
                    }
                    else {
                        buffers[which].values[in_pos] += buffer_vals[out_pos];
                    }
                    break;
                }
                case LossyOpCode::SUB: {
                    std::uint8_t const which = on_which_buffer[on_which_buffer_idx++];
                    std::size_t const out_pos = pos[id_idx_l++];
                    std::size_t const in_pos = pos[id_idx_l++];
                    if constexpr (Vectorised) {
                        const double* src = &buffer_vals[out_pos * this->m_num_lanes];
                        double* dest = &buffers[which].values[in_pos * this->m_num_lanes];
#pragma omp simd
                        for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                            dest[i] -= src[i];
                        }
                    }
                    else {
                        buffers[which].values[in_pos] -= buffer_vals[out_pos];
                    }
                    break;
                }
                case LossyOpCode::MINUS_INPLACE: {
                    std::size_t const out_pos = pos[id_idx_l++];
                    if constexpr (Vectorised) {
                        double* dest = &buffer_vals[out_pos * this->m_num_lanes];
#pragma omp simd
                        for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                            dest[i] = -dest[i];
                        }
                    }
                    else {
                        buffer_vals[out_pos] = -buffer_vals[out_pos];
                    }
                    break;
                }
                case LossyOpCode::MUL_INPLACE: {
                    std::size_t const out_pos = pos[id_idx_l++];
                    double const multiplier = value_fetcher();
                    if constexpr (Vectorised) {
                        double* dest = &buffer_vals[out_pos * this->m_num_lanes];
#pragma omp simd
                        for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                            dest[i] *= multiplier;
                        }
                    }
                    else {
                        buffer_vals[out_pos] *= multiplier;
                    }
                    break;
                }
                case LossyOpCode::MUL_ADD: {
                    std::uint8_t const which = on_which_buffer[on_which_buffer_idx++];
                    std::size_t const out_pos = pos[id_idx_l++];
                    std::size_t const in_pos = pos[id_idx_l++];
                    double const multiplier = value_fetcher();
                    if constexpr (Vectorised) {
                        const double* src = &buffer_vals[out_pos * this->m_num_lanes];
                        double* dest = &buffers[which].values[in_pos * this->m_num_lanes];
#pragma omp simd
                        for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                            dest[i] += src[i] * multiplier;
                        }
                    }
                    else {
                        buffers[which].values[in_pos] += buffer_vals[out_pos] * multiplier;
                    }
                    break;
                }
                case LossyOpCode::MUL_SET: {
                    std::size_t const out_pos = pos[id_idx_l++];
                    std::size_t const in_pos = pos[id_idx_l++];
                    double const multiplier = value_fetcher();
                    if constexpr (Vectorised) {
                        const double* src = &buffer_vals[out_pos * this->m_num_lanes];
                        double* dest = &buffer_vals[in_pos * this->m_num_lanes];
#pragma omp simd
                        for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                            dest[i] = src[i] * multiplier;
                        }
                    }
                    else {
                        buffer_vals[in_pos] = buffer_vals[out_pos] * multiplier;
                    }
                    break;
                }
            }
        }
    }
    else {

        if (buffers.size() != 1) {
            throw;
        }

        auto& buffer_vals = this->buffers[0].values;
        auto& buffer_free_positions = buffers[0].free_positions;

        auto copy = [&](std::size_t const out_pos, std::size_t const in_pos) {
            if constexpr (Vectorised) {
                const double* src = &buffer_vals[out_pos * this->m_num_lanes];
                double* dest = &buffer_vals[in_pos * this->m_num_lanes];
#pragma omp simd
                for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                    dest[i] = src[i];
                }
            }
            else {
                buffer_vals[in_pos] = buffer_vals[out_pos];
            }
        };

        auto copy_minus = [&](std::size_t const out_pos, std::size_t const in_pos) {
            if constexpr (Vectorised) {
                const double* src = &buffer_vals[out_pos * this->m_num_lanes];
                double* dest = &buffer_vals[in_pos * this->m_num_lanes];
#pragma omp simd
                for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                    dest[i] = -src[i];
                }
            }
            else {
                buffer_vals[in_pos] = -buffer_vals[out_pos];
            }
        };

        auto add = [&](std::size_t const out_pos, std::size_t const in_pos) {
            if constexpr (Vectorised) {
                const double* src = &buffer_vals[out_pos * this->m_num_lanes];
                double* dest = &buffer_vals[in_pos * this->m_num_lanes];
#pragma omp simd
                for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                    dest[i] += src[i];
                }
            }
            else {
                buffer_vals[in_pos] += buffer_vals[out_pos];
            }
        };

        auto sub = [&](std::size_t const out_pos, std::size_t const in_pos) {
            if constexpr (Vectorised) {
                const double* src = &buffer_vals[out_pos * this->m_num_lanes];
                double* dest = &buffer_vals[in_pos * this->m_num_lanes];
#pragma omp simd
                for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                    dest[i] -= src[i];
                }
            }
            else {
                buffer_vals[in_pos] -= buffer_vals[out_pos];
            }
        };

        auto minus_inplace = [&](std::size_t const pos) {
            if constexpr (Vectorised) {
                double* dest = &buffer_vals[pos * this->m_num_lanes];
#pragma omp simd
                for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                    dest[i] = -dest[i];
                }
            }
            else {
                buffer_vals[pos] = -buffer_vals[pos];
            }
        };

        auto mul_inplace = [&](std::size_t const pos, double const multiplier) {
            if constexpr (Vectorised) {
                double* dest = &buffer_vals[pos * this->m_num_lanes];
#pragma omp simd
                for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                    dest[i] *= multiplier;
                }
            }
            else {
                buffer_vals[pos] *= multiplier;
            }
        };

        auto mul_add = [&](std::size_t const out_pos, std::size_t const in_pos, double const multiplier) {
            if constexpr (Vectorised) {
                const double* src = &buffer_vals[out_pos * this->m_num_lanes];
                double* dest = &buffer_vals[in_pos * this->m_num_lanes];
#pragma omp simd
                for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                    dest[i] += src[i] * multiplier;
                }
            }
            else {
                buffer_vals[in_pos] += buffer_vals[out_pos] * multiplier;
            }
        };

        auto mul_set = [&](std::size_t const out_pos, std::size_t const in_pos, double const multiplier) {
            if constexpr (Vectorised) {
                const double* src = &buffer_vals[out_pos * this->m_num_lanes];
                double* dest = &buffer_vals[in_pos * this->m_num_lanes];
#pragma omp simd
                for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                    dest[i] = src[i] * multiplier;
                }
            }
            else {
                buffer_vals[in_pos] = buffer_vals[out_pos] * multiplier;
            }
        };

        auto update_loc = [this](std::size_t& arg_pos) {
            if (arg_pos == passive_id<std::size_t>) {
                auto& arg_buffer = this->buffers[0];
                if (arg_buffer.free_positions.empty()) {
                    arg_pos = arg_buffer.size;
                    ++arg_buffer.size;
                    arg_buffer.values.resize(arg_buffer.values.size() + this->m_num_lanes);
                }
                else {
                    arg_pos = arg_buffer.free_positions.back();
                    arg_buffer.free_positions.pop_back();
                }
            }
        };

        auto copy_add = [copy, add, update_loc](std::size_t res_pos, std::size_t& arg_pos) {
            bool arg_is_new = (arg_pos == passive_id<std::size_t>);
            update_loc(arg_pos);

            if (arg_is_new) {
                copy(res_pos, arg_pos);
            }
            else {
                add(res_pos, arg_pos);
            }
        };

        auto copy_sub = [copy_minus, sub, update_loc](std::size_t res_pos, std::size_t& arg_pos) {
            bool arg_is_new = (arg_pos == passive_id<std::size_t>);
            update_loc(arg_pos);

            if (arg_is_new) {
                // this is a new value, we NEED to override
                copy_minus(res_pos, arg_pos);
            }
            else {
                sub(res_pos, arg_pos);
            }
        };

        auto copy_mul = [mul_set, mul_add, update_loc](std::size_t res_pos, std::size_t& arg_pos, double multiplier) {
            bool arg_is_new = (arg_pos == passive_id<std::size_t>);
            update_loc(arg_pos);

            if (arg_is_new) {
                // this is a new value, we NEED to override
                mul_set(res_pos, arg_pos, multiplier);
            }
            else {
                mul_add(res_pos, arg_pos, multiplier);
            }
        };

        auto update_univariate =
          [mul_inplace,
           this,
           copy_mul,
           &node_location_on_buffer = this->node_location_on_buffer,
           &buffer_free_positions](std::size_t const arg_id, std::size_t const res_id, double const der_local_1) {
              std::size_t& res_pos = node_location_on_buffer[res_id];
              std::size_t& arg_pos = node_location_on_buffer[arg_id];

              bool const arg_inplace = (arg_pos == passive_id<std::size_t>);

              if (arg_inplace) {
                  mul_inplace(res_pos, der_local_1);
                  // res id should now be arg id, avoiding a copy and a potential buffer increase
                  std::swap(arg_pos, res_pos);
              }
              else {
                  copy_mul(res_pos, arg_pos, der_local_1);
                  buffer_free_positions.push_back(res_pos);
                  res_pos = passive_id<std::size_t>;
              }
          };

        for (std::size_t op_idx = from; op_idx-- > to;) {
            OpCode const& op = ops[op_idx];
            bool const use_this_op = this->node_location_on_buffer[op_idx] != passive_id<std::size_t>;

            switch (op) {
                case OpCode::REG_INPUT: {
                    id_idx -= 1;
                    if constexpr (Reset) {
                        std::size_t const id = ids[id_idx];
                        std::size_t& pos = node_location_on_buffer[id];
                        buffer_free_positions.push_back(pos);
                        pos = passive_id<std::size_t>;
                    }
                    break;
                }
                case OpCode::REG_OUTPUT: {
                    id_idx -= 2;
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];

                        std::size_t& res_pos = node_location_on_buffer[res_id];
                        std::size_t& arg_pos = node_location_on_buffer[arg_id];

                        if constexpr (Reset) {
                            bool const arg_inplace = (arg_pos == passive_id<std::size_t>);
                            if (arg_inplace) {
                                // res id should now be lhs id, avoiding a copy and a potential buffer increase
                                std::swap(arg_pos, res_pos);
                            }
                            else {
                                copy_add(res_pos, arg_pos);
                                buffer_free_positions.push_back(res_pos);
                                res_pos = passive_id<std::size_t>;
                            }
                        }
                        else {
                            copy_add(res_pos, arg_pos);
                            // we don't reset so we don't free res_pos,
                            // it will be potentially used in the next operations as an active node location
                        }
                    }
                    break;
                }
                case OpCode::ADD: {
                    id_idx -= 3;
                    if (use_this_op) {
                        std::size_t const lhs_id = ids[id_idx];
                        std::size_t const rhs_id = ids[id_idx + 1];
                        std::size_t const res_id = ids[id_idx + 2];

                        std::size_t& res_pos = node_location_on_buffer[res_id];
                        std::size_t& lhs_pos = node_location_on_buffer[lhs_id];
                        std::size_t& rhs_pos = node_location_on_buffer[rhs_id];

                        bool const lhs_inplace = (lhs_pos == passive_id<std::size_t>);
                        bool const rhs_inplace = !lhs_inplace && (rhs_pos == passive_id<std::size_t>);

                        if (!lhs_inplace) {
                            copy_add(res_pos, lhs_pos);
                        }

                        if (!rhs_inplace) {
                            copy_add(res_pos, rhs_pos);
                        }

                        if (lhs_inplace) {
                            // res id should now be lhs id, avoiding a copy and a potential buffer increase
                            std::swap(lhs_pos, res_pos);
                        }
                        else if (rhs_inplace) {
                            // res id should now be rhs id, avoiding a copy and a potential buffer increase
                            std::swap(rhs_pos, res_pos);
                        }
                        else {
                            // don't forget to free res_id from the buffer!
                            buffer_free_positions.push_back(res_pos);
                            res_pos = passive_id<std::size_t>;
                        }
                    }
                    break;
                }
                case OpCode::SUB: {
                    id_idx -= 3;
                    if (use_this_op) {
                        std::size_t const lhs_id = ids[id_idx];
                        std::size_t const rhs_id = ids[id_idx + 1];
                        std::size_t const res_id = ids[id_idx + 2];

                        std::size_t& res_pos = node_location_on_buffer[res_id];
                        std::size_t& lhs_pos = node_location_on_buffer[lhs_id];
                        std::size_t& rhs_pos = node_location_on_buffer[rhs_id];

                        bool const lhs_inplace = (lhs_pos == passive_id<std::size_t>);
                        bool const rhs_inplace = !lhs_inplace && (rhs_pos == passive_id<std::size_t>);

                        if (!lhs_inplace) {
                            copy_add(res_pos, lhs_pos);
                        }

                        if (!rhs_inplace) {
                            copy_sub(res_pos, rhs_pos);
                        }

                        if (lhs_inplace) {
                            // res id should now be lhs id, avoiding a copy and a potential buffer increase
                            std::swap(lhs_pos, res_pos);
                        }
                        else if (rhs_inplace) {
                            // this is a subtraction, so we need to negate the value in the buffer
                            minus_inplace(res_pos);
                            // res id should now be rhs id, avoiding a copy and a potential buffer increase
                            std::swap(rhs_pos, res_pos);
                        }
                        else {
                            // don't forget to free res_id from the buffer!
                            buffer_free_positions.push_back(res_pos);
                            res_pos = passive_id<std::size_t>;
                        }
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

                        std::size_t& res_pos = node_location_on_buffer[res_id];
                        std::size_t& lhs_pos = node_location_on_buffer[lhs_id];
                        std::size_t& rhs_pos = node_location_on_buffer[rhs_id];

                        bool const lhs_inplace = (lhs_pos == passive_id<std::size_t>);
                        bool const rhs_inplace = !lhs_inplace && (rhs_pos == passive_id<std::size_t>);

                        if (!lhs_inplace) {
                            copy_mul(res_pos, lhs_pos, rhs_val);
                        }

                        if (!rhs_inplace) {
                            copy_mul(res_pos, rhs_pos, lhs_val);
                        }

                        if (lhs_inplace) {
                            mul_inplace(res_pos, rhs_val);
                            // res id should now be lhs id, avoiding a copy and a potential buffer increase
                            std::swap(lhs_pos, res_pos);
                        }
                        else if (rhs_inplace) {
                            mul_inplace(res_pos, lhs_val);
                            // res id should now be rhs id, avoiding a copy and a potential buffer increase
                            std::swap(rhs_pos, res_pos);
                        }
                        else {
                            // don't forget to free res_id from the buffer!
                            buffer_free_positions.push_back(res_pos);
                            res_pos = passive_id<std::size_t>;
                        }
                    }
                    break;
                }
                case OpCode::ADD_C: {
                    id_idx -= 2;
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];

                        std::size_t& res_pos = node_location_on_buffer[res_id];
                        std::size_t& arg_pos = node_location_on_buffer[arg_id];

                        bool const arg_inplace = (arg_pos == passive_id<std::size_t>);

                        if (arg_inplace) {
                            // res id should now be lhs id, avoiding a copy and a potential buffer increase
                            std::swap(arg_pos, res_pos);
                        }
                        else {
                            copy_add(res_pos, arg_pos);
                            buffer_free_positions.push_back(res_pos);
                            res_pos = passive_id<std::size_t>;
                        }
                    }
                    break;
                }
                case OpCode::SUB_C: {
                    id_idx -= 2;
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];

                        std::size_t& res_pos = node_location_on_buffer[res_id];
                        std::size_t& arg_pos = node_location_on_buffer[arg_id];

                        bool const arg_inplace = (arg_pos == passive_id<std::size_t>);

                        if (arg_inplace) {
                            // this is a subtraction, so we need to negate the value in the buffer
                            minus_inplace(res_pos);
                            // res id should now be lhs id, avoiding a copy and a potential buffer increase
                            std::swap(arg_pos, res_pos);
                        }
                        else {
                            copy_sub(res_pos, arg_pos);
                            buffer_free_positions.push_back(res_pos);
                            res_pos = passive_id<std::size_t>;
                        }
                    }
                    break;
                }
                case OpCode::MUL_C: {
                    val_idx -= 1;
                    id_idx -= 2;
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        double const der_local_1 = vals[val_idx];
                        update_univariate(arg_id, res_id, der_local_1);
                    }
                    break;
                }
                case OpCode::NORM: {
                    val_idx -= 1;
                    id_idx -= 2;
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        double const der_local_1 = 2.0 * vals[val_idx];
                        update_univariate(arg_id, res_id, der_local_1);
                    }
                    break;
                }
                case OpCode::INV: {
                    val_idx -= 1;
                    id_idx -= 2;
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        double const der_local_1 = -vals[val_idx] * vals[val_idx];
                        update_univariate(arg_id, res_id, der_local_1);
                    }
                    break;
                }
                case OpCode::ABS: {
                    val_idx -= 1;
                    id_idx -= 2;
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        double const der_local_1 = std::copysign(1.0, vals[val_idx]);
                        update_univariate(arg_id, res_id, der_local_1);
                    }
                    break;
                }
                case OpCode::EXP: {
                    val_idx -= 1;
                    id_idx -= 2;
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        double const der_local_1 = vals[val_idx];
                        update_univariate(arg_id, res_id, der_local_1);
                    }
                    break;
                }
                case OpCode::LOG: {
                    val_idx -= 1;
                    id_idx -= 2;
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        double const der_local_1 = 1.0 / vals[val_idx];
                        update_univariate(arg_id, res_id, der_local_1);
                    }
                    break;
                }
                case OpCode::ERF: {
                    val_idx -= 1;
                    id_idx -= 2;
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        constexpr double two_over_root_pi = 2. * std::numbers::inv_sqrtpi_v<double>;
                        double const der_local_1 = std::exp(-vals[val_idx] * vals[val_idx]) * two_over_root_pi;
                        update_univariate(arg_id, res_id, der_local_1);
                    }
                    break;
                }
                case OpCode::ERFC: {
                    val_idx -= 1;
                    id_idx -= 2;
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        constexpr double minus_two_over_root_pi = -2. * std::numbers::inv_sqrtpi_v<double>;
                        double const der_local_1 = std::exp(-vals[val_idx] * vals[val_idx]) * minus_two_over_root_pi;
                        update_univariate(arg_id, res_id, der_local_1);
                    }
                    break;
                }
                case OpCode::COS: {
                    val_idx -= 2;
                    id_idx -= 2;
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        double const der_local_1 = -std::sin(vals[val_idx]);
                        update_univariate(arg_id, res_id, der_local_1);
                    }
                    break;
                }
                case OpCode::SQRT: {
                    val_idx -= 2;
                    id_idx -= 2;
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        double const one_over_in = 1. / vals[val_idx];
                        double const der_local_1 = 0.5 * vals[val_idx + 1] * one_over_in;
                        update_univariate(arg_id, res_id, der_local_1);
                    }
                    break;
                }
                case OpCode::POW_C: {
                    val_idx -= 2;
                    id_idx -= 2;
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        double const lhs_arg = vals[val_idx];
                        double const rhs_arg = vals[val_idx + 1];
                        double const der_local_1 = rhs_arg != 0.0 ? rhs_arg * std::pow(lhs_arg, rhs_arg - 1.) : 0.0;
                        update_univariate(arg_id, res_id, der_local_1);
                    }
                    break;
                }
            }
        }
    }

    if constexpr (Reset) {
        if (to == this->checkpoints.back()) {
            this->buffers.back() = {};
        }

        reset(pos, data);
        this->node_location_on_buffer.resize(pos.op_position);
    }
}

} // namespace adhoc

#endif // ADHOC_BACKPROPAGATOR1LOSSYCOMPRESSEDPATHREUSE_HPP
