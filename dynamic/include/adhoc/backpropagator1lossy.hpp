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

#ifndef ADHOC_BACKPROPAGATOR1LOSSY_HPP
#define ADHOC_BACKPROPAGATOR1LOSSY_HPP

#include "passive_id.hpp"
#include "tape_data.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

namespace adhoc {

template<class Float, bool Vectorised = false>
class BackPropagatorLossy {
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

  public:
    explicit BackPropagatorLossy() = default;

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
        size += 1 * sizeof(std::size_t); // m_num_lanes
        size += sizeof(std::size_t) * (capacity ? node_location_on_buffer.capacity() : node_location_on_buffer.size());
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
BackPropagatorLossy<Float, Vectorised>::backpropagate_to(std::size_t to, TapeData& data)
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

    this->node_location_on_buffer.resize(ops.size(), passive_id<std::size_t>);

    auto& buffer_vals = buffers.back().values;
    auto& buffer_free_positions = buffers.back().free_positions;

    auto copy = [&](std::uint8_t const which, std::size_t const out_pos, std::size_t const in_pos) {
        if constexpr (Vectorised) {
            const double* src = &buffer_vals[out_pos * this->m_num_lanes];
            double* dest = &buffers[which].values[in_pos * this->m_num_lanes];
#pragma omp simd
            for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                dest[i] = src[i];
            }
        }
        else {
            buffers[which].values[in_pos] = buffer_vals[out_pos];
        }
    };

    auto copy_minus = [&](std::uint8_t const which, std::size_t const out_pos, std::size_t const in_pos) {
        if constexpr (Vectorised) {
            const double* src = &buffer_vals[out_pos * this->m_num_lanes];
            double* dest = &buffers[which].values[in_pos * this->m_num_lanes];
#pragma omp simd
            for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                dest[i] = -src[i];
            }
        }
        else {
            buffers[which].values[in_pos] = -buffer_vals[out_pos];
        }
    };

    auto add = [&](std::uint8_t const which, std::size_t const out_pos, std::size_t const in_pos) {
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
    };

    auto sub = [&](std::uint8_t const which, std::size_t const out_pos, std::size_t const in_pos) {
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

    auto mul_add =
      [&](std::uint8_t const which, std::size_t const out_pos, std::size_t const in_pos, double const multiplier) {
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
      };

    auto mul_set =
      [&](std::uint8_t const which, std::size_t const out_pos, std::size_t const in_pos, double const multiplier) {
          if constexpr (Vectorised) {
              const double* src = &buffer_vals[out_pos * this->m_num_lanes];
              double* dest = &buffers[which].values[in_pos * this->m_num_lanes];
#pragma omp simd
              for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                  dest[i] = src[i] * multiplier;
              }
          }
          else {
              buffers[which].values[in_pos] = buffer_vals[out_pos] * multiplier;
          }
      };

    // output is: is_new, is_current, buffer_id, position
    auto const& checkpoints_c = this->checkpoints;
    auto get_loc = [checkpoints_c](std::size_t id) -> std::tuple<bool, std::uint8_t> {
        auto it = std::upper_bound(checkpoints_c.begin(), checkpoints_c.end(), id);
        auto buffer_id = static_cast<std::uint8_t>(std::distance(checkpoints_c.cbegin(), it) - 1);
        return { it == checkpoints_c.end(), buffer_id };
    };

    auto update_loc = [this, &buffer_vals, &buffers = this->buffers](std::size_t& arg_pos, std::uint8_t buffer_id) {
        if (arg_pos == passive_id<std::size_t>) {
            auto& arg_buffer = buffers[buffer_id];
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

    auto copy_add = [copy, add, update_loc](std::size_t res_pos, std::size_t& arg_pos, std::uint8_t buffer_id) {
        bool arg_is_new = (arg_pos == passive_id<std::size_t>);
        update_loc(arg_pos, buffer_id);

        if (arg_is_new) {
            copy(buffer_id, res_pos, arg_pos);
        }
        else {
            add(buffer_id, res_pos, arg_pos);
        }
    };

    auto copy_sub = [copy_minus, sub, update_loc](std::size_t res_pos, std::size_t& arg_pos, std::uint8_t buffer_id) {
        bool arg_is_new = (arg_pos == passive_id<std::size_t>);
        update_loc(arg_pos, buffer_id);

        if (arg_is_new) {
            // this is a new value, we NEED to override
            copy_minus(buffer_id, res_pos, arg_pos);
        }
        else {
            sub(buffer_id, res_pos, arg_pos);
        }
    };

    auto copy_mul = [mul_set, mul_add, update_loc](
                      std::size_t res_pos, std::size_t& arg_pos, std::uint8_t buffer_id, double multiplier) {
        bool arg_is_new = (arg_pos == passive_id<std::size_t>);
        update_loc(arg_pos, buffer_id);

        if (arg_is_new) {
            // this is a new value, we NEED to override
            mul_set(buffer_id, res_pos, arg_pos, multiplier);
        }
        else {
            mul_add(buffer_id, res_pos, arg_pos, multiplier);
        }
    };

    auto update_univariate =
      [mul_inplace,
       this,
       copy_mul,
       get_loc,
       &node_location_on_buffer = this->node_location_on_buffer,
       &buffer_free_positions](std::size_t const arg_id, std::size_t const res_id, double const der_local_1) {
          auto const res_pos = node_location_on_buffer[res_id];
          node_location_on_buffer[res_id] = passive_id<std::size_t>;
          std::size_t& arg_pos = node_location_on_buffer[arg_id];

          auto const arg_pos_data = get_loc(arg_id);
          bool arg_is_new = (arg_pos == passive_id<std::size_t>);
          bool arg_inplace = arg_is_new && std::get<0>(arg_pos_data);

          if (arg_inplace) {
              // res id should now be arg id, avoiding a copy and a potential buffer increase
              arg_pos = res_pos;
              mul_inplace(res_pos, der_local_1);
          }
          else {
              copy_mul(res_pos, arg_pos, std::get<1>(arg_pos_data), der_local_1);
              buffer_free_positions.push_back(res_pos);
          }
      };

    for (std::size_t op_idx = from; op_idx-- > to;) {
        OpCode const& op = ops[op_idx];
        bool const use_this_op = this->node_location_on_buffer[op_idx] != passive_id<std::size_t>;

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

                    auto const res_pos = node_location_on_buffer[res_id];
                    std::size_t& arg_pos = node_location_on_buffer[arg_id];

                    auto const arg_pos_data = get_loc(arg_id);
                    copy_add(res_pos, arg_pos, std::get<1>(arg_pos_data));
                }
                break;
            }
            case OpCode::ADD: {
                id_idx -= 3;
                if (use_this_op) {
                    std::size_t const lhs_id = ids[id_idx];
                    std::size_t const rhs_id = ids[id_idx + 1];
                    std::size_t const res_id = ids[id_idx + 2];

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
                        copy_add(res_pos, lhs_pos, std::get<1>(lhs_pos_data));
                    }

                    if (!rhs_inplace) {
                        copy_add(res_pos, rhs_pos, std::get<1>(rhs_pos_data));
                    }

                    if (lhs_inplace) {
                        // res id should now be lhs id, avoiding a copy and a potential buffer increase
                        lhs_pos = res_pos;
                    }
                    else if (rhs_inplace) {
                        // res id should now be rhs id, avoiding a copy and a potential buffer increase
                        rhs_pos = res_pos;
                    }
                    else {
                        // don't forget to free res_id from the buffer!
                        buffer_free_positions.push_back(res_pos);
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
                        copy_add(res_pos, lhs_pos, std::get<1>(lhs_pos_data));
                    }

                    if (!rhs_inplace) {
                        copy_sub(res_pos, rhs_pos, std::get<1>(rhs_pos_data));
                    }

                    if (lhs_inplace) {
                        // res id should now be lhs id, avoiding a copy and a potential buffer increase
                        lhs_pos = res_pos;
                    }
                    else if (rhs_inplace) {
                        // res id should now be rhs id, avoiding a copy and a potential buffer increase
                        rhs_pos = res_pos;
                        // however this is a subtraction, so we need to negate the value in the buffer
                        minus_inplace(rhs_pos);
                    }
                    else {
                        // don't forget to free res_id from the buffer!
                        buffer_free_positions.push_back(res_pos);
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
                        copy_mul(res_pos, lhs_pos, std::get<1>(lhs_pos_data), rhs_val);
                    }

                    if (!rhs_inplace) {
                        copy_mul(res_pos, rhs_pos, std::get<1>(rhs_pos_data), lhs_val);
                    }

                    if (lhs_inplace) {
                        // res id should now be lhs id, avoiding a copy and a potential buffer increase
                        lhs_pos = res_pos;
                        mul_inplace(lhs_pos, rhs_val);
                    }
                    else if (rhs_inplace) {
                        // res id should now be rhs id, avoiding a copy and a potential buffer increase
                        rhs_pos = res_pos;
                        mul_inplace(rhs_pos, lhs_val);
                    }
                    else {
                        // don't forget to free res_id from the buffer!
                        buffer_free_positions.push_back(res_pos);
                    }
                }
                break;
            }
            case OpCode::ADD_C: {
                id_idx -= 2;
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];

                    auto const res_pos = node_location_on_buffer[res_id];
                    node_location_on_buffer[res_id] = passive_id<std::size_t>;
                    std::size_t& arg_pos = node_location_on_buffer[arg_id];

                    auto const arg_pos_data = get_loc(arg_id);
                    bool arg_is_new = (arg_pos == passive_id<std::size_t>);
                    bool arg_inplace = arg_is_new && std::get<0>(arg_pos_data);

                    if (arg_inplace) {
                        // res id should now be lhs id, avoiding a copy and a potential buffer increase
                        arg_pos = res_pos;
                    }
                    else {
                        copy_add(res_pos, arg_pos, std::get<1>(arg_pos_data));
                        buffer_free_positions.push_back(res_pos);
                    }
                }
                break;
            }
            case OpCode::SUB_C: {
                id_idx -= 2;
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];

                    auto const res_pos = node_location_on_buffer[res_id];
                    node_location_on_buffer[res_id] = passive_id<std::size_t>;
                    std::size_t& arg_pos = node_location_on_buffer[arg_id];

                    auto const arg_pos_data = get_loc(arg_id);
                    bool arg_is_new = (arg_pos == passive_id<std::size_t>);
                    bool arg_inplace = arg_is_new && std::get<0>(arg_pos_data);

                    if (arg_inplace) {
                        // res id should now be lhs id, avoiding a copy and a potential buffer increase
                        arg_pos = res_pos;
                        // however this is a subtraction, so we need to negate the value in the buffer
                        minus_inplace(arg_pos);
                    }
                    else {
                        copy_sub(res_pos, arg_pos, std::get<1>(arg_pos_data));
                        buffer_free_positions.push_back(res_pos);
                    }
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
                    update_univariate(arg_id, res_id, der_local_1);
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
                    update_univariate(arg_id, res_id, der_local_1);
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
                    update_univariate(arg_id, res_id, der_local_1);
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
                    update_univariate(arg_id, res_id, der_local_1);
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
                    update_univariate(arg_id, res_id, der_local_1);
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
                    update_univariate(arg_id, res_id, der_local_1);
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
                    update_univariate(arg_id, res_id, der_local_1);
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
                    update_univariate(arg_id, res_id, der_local_1);
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
                    update_univariate(arg_id, res_id, der_local_1);
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
                    update_univariate(arg_id, res_id, der_local_1);
                }
                break;
            }
            case OpCode::POW_C: {
                val_idx -= 2;
                id_idx -= 2;
                if (use_this_op) {
                    double const lhs_arg = vals[val_idx];
                    double const rhs_arg = vals[val_idx + 1];
                    double const der_local_1 = rhs_arg != 0.0 ? rhs_arg * std::pow(lhs_arg, rhs_arg - 1.) : 0.0;
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    update_univariate(arg_id, res_id, der_local_1);
                }
                break;
            }
        }
    }

    if constexpr (Reset) {
        data.ops.resize(to);
        data.vals.resize(val_idx);
        data.ids.resize(id_idx);
        data.next_id = to;

        this->node_location_on_buffer.resize(to);
    }

    if constexpr (Reset) {
        if (!checkpoints.empty() && to == checkpoints.back()) {
            // buffers.pop_back();
            buffers.back() = buffer_t{};
        }
    }
}

} // namespace adhoc

#endif // ADHOC_BACKPROPAGATOR1LOSSY_HPP
