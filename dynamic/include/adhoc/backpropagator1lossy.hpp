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

#include "buffer_t.hpp"
#include "passive_id.hpp"
#include "position_impl.hpp"
#include "tape_data.hpp"

#include <algorithm>
#include <cmath>
#include <concepts>
#include <numbers>
#include <vector>

namespace adhoc {

template<std::floating_point Float, bool Vectorised = false>
class BackPropagatorLossy {
  private:
    std::vector<std::size_t> node_location_on_buffer;
    buffer_t<Float, Vectorised> buffer;

  public:
    explicit BackPropagatorLossy() = default;

    void set_checkpoint(std::size_t /* ops_size */) {}

    void set_lanes(std::size_t num_lanes) { this->buffer = buffer_t<Float, Vectorised>{ num_lanes }; }

    auto get_lanes() const -> std::size_t { return this->buffer.lanes(); }

    void reserve_input(std::size_t count_registered) { this->buffer.reserve(this->buffer.size() + count_registered); }

    void reserve_output(std::size_t count_registered) { this->buffer.reserve(this->buffer.size() + count_registered); }

    void register_variable(std::size_t var_id)
    {
        this->node_location_on_buffer.resize(std::max(this->node_location_on_buffer.size(), var_id + 1),
                                             passive_id<std::size_t>);
        std::size_t& var_pos = this->node_location_on_buffer[var_id];
        if (var_pos == passive_id<std::size_t>) {
            var_pos = this->buffer.template get_new_loc<true>();
        }
    }

    void register_output_variable(std::size_t var_id, std::size_t /* ops_size */)
    {
        this->node_location_on_buffer.resize(std::max(this->node_location_on_buffer.size(), var_id + 1),
                                             passive_id<std::size_t>);
        std::size_t& var_pos = this->node_location_on_buffer[var_id];
        if (var_pos == passive_id<std::size_t>) {
            var_pos = this->buffer.template get_new_loc<true>();
        }
    }

    void set_derivative(std::size_t var_id, double deriv, std::size_t /* ops_size */, std::size_t lane = 0)
    {

        std::size_t const var_pos = this->node_location_on_buffer[var_id];
        if (var_pos == passive_id<std::size_t>) {
            // this derivative is not on buffer.
            // this is probably nor an input nor an output.
            throw;
        }

        if constexpr (Vectorised) {
            auto val = this->buffer[var_pos];

            if (lane >= val.size()) {
                // lane value too large
                throw;
            }

            val[lane] = static_cast<Float>(deriv);
        }
        else {
            if (lane != 0) {
                // This backpropagator is not designed for multiple lanes
                throw;
            }
            this->buffer[var_pos] = static_cast<Float>(deriv);
        }
    }

    void set_derivative(std::size_t /* var_id1 */, std::size_t /* var_id2 */, double /* deriv */) { throw; }

    auto get_derivative(std::size_t var_id, std::size_t lane) const -> double
    {

        std::size_t const var_pos = this->node_location_on_buffer[var_id];
        if (var_pos == passive_id<std::size_t>) {
            // this derivative is not on buffer.
            // this is probably nor an input nor an output.
            throw;
        }

        if constexpr (Vectorised) {
            auto const val = this->buffer[var_pos];

            if (lane >= val.size()) {
                // lane value too large
                throw;
            }

            return static_cast<double>(val[lane]);
        }
        else {
            if (lane != 0) {
                // This backpropagator is not designed for multiple lanes
                throw;
            }
            return static_cast<double>(this->buffer[var_pos]);
        }
    }

    auto get_derivative(std::size_t /* var_id1 */, std::size_t /* var_id2 */, std::size_t /* lane */) const -> double
    {
        return 0.;
    }

    void clear() {}

    void zero_adjoints() { this->buffer.zero(); }

    auto size_of(bool capacity = false) const -> std::size_t
    {
        std::size_t size = 0;
        size += sizeof(std::size_t) * (capacity ? node_location_on_buffer.capacity() : node_location_on_buffer.size());
        size += this->buffer.size_of(capacity);
        return size;
    }

    void reset(PositionImpl const& /* pos */)
    {
        // empty for now
    }

    template<bool Reset, class TapeDataType>
    void backpropagate_to(PositionImpl const& pos, TapeDataType const& data);
};

template<std::floating_point Float, bool Vectorised>
template<bool Reset, class TapeDataType>
void
BackPropagatorLossy<Float, Vectorised>::backpropagate_to(PositionImpl const& pos, TapeDataType const& data)
{
    std::size_t to = pos.op_position;
    std::size_t from = data.next_id;

    const auto& ops = data.ops;
    const auto& vals = data.vals;
    const auto& ids = data.ids;

    std::size_t val_idx = vals.size();
    std::size_t id_idx = ids.size();

    this->node_location_on_buffer.resize(ops.size(), passive_id<std::size_t>);

    auto copy = [&](std::size_t const out_pos, std::size_t const in_pos) {
        if constexpr (Vectorised) {
            auto const src = this->buffer[out_pos];
            auto dest = this->buffer[in_pos];
#pragma omp simd
            for (std::size_t i = 0; i < dest.size(); ++i) {
                dest[i] = src[i];
            }
        }
        else {
            this->buffer[in_pos] = this->buffer[out_pos];
        }
    };

    auto copy_minus = [&](std::size_t const out_pos, std::size_t const in_pos) {
        if constexpr (Vectorised) {
            auto const src = this->buffer[out_pos];
            auto dest = this->buffer[in_pos];
#pragma omp simd
            for (std::size_t i = 0; i < dest.size(); ++i) {
                dest[i] = -src[i];
            }
        }
        else {
            this->buffer[in_pos] = -this->buffer[out_pos];
        }
    };

    auto add = [&](std::size_t const out_pos, std::size_t const in_pos) {
        if constexpr (Vectorised) {
            auto const src = this->buffer[out_pos];
            auto dest = this->buffer[in_pos];
#pragma omp simd
            for (std::size_t i = 0; i < dest.size(); ++i) {
                dest[i] += src[i];
            }
        }
        else {
            this->buffer[in_pos] += this->buffer[out_pos];
        }
    };

    auto sub = [&](std::size_t const out_pos, std::size_t const in_pos) {
        if constexpr (Vectorised) {
            auto const src = this->buffer[out_pos];
            auto dest = this->buffer[in_pos];
#pragma omp simd
            for (std::size_t i = 0; i < dest.size(); ++i) {
                dest[i] -= src[i];
            }
        }
        else {
            this->buffer[in_pos] -= this->buffer[out_pos];
        }
    };

    auto minus_inplace = [&](std::size_t const pos) {
        if constexpr (Vectorised) {
            auto dest = this->buffer[pos];
#pragma omp simd
            for (std::size_t i = 0; i < dest.size(); ++i) {
                dest[i] = -dest[i];
            }
        }
        else {
            this->buffer[pos] = -this->buffer[pos];
        }
    };

    auto mul_inplace = [&](std::size_t const pos, Float const multiplier) {
        if constexpr (Vectorised) {
            auto dest = this->buffer[pos];
#pragma omp simd
            for (std::size_t i = 0; i < dest.size(); ++i) {
                dest[i] *= multiplier;
            }
        }
        else {
            this->buffer[pos] *= multiplier;
        }
    };

    auto mul_add = [&](std::size_t const out_pos, std::size_t const in_pos, Float const multiplier) {
        if constexpr (Vectorised) {
            auto const src = this->buffer[out_pos];
            auto dest = this->buffer[in_pos];
#pragma omp simd
            for (std::size_t i = 0; i < dest.size(); ++i) {
                dest[i] += src[i] * multiplier;
            }
        }
        else {
            this->buffer[in_pos] += this->buffer[out_pos] * multiplier;
        }
    };

    auto mul_set = [&](std::size_t const out_pos, std::size_t const in_pos, Float const multiplier) {
        if constexpr (Vectorised) {
            auto const src = this->buffer[out_pos];
            auto dest = this->buffer[in_pos];
#pragma omp simd
            for (std::size_t i = 0; i < dest.size(); ++i) {
                dest[i] = src[i] * multiplier;
            }
        }
        else {
            this->buffer[in_pos] = this->buffer[out_pos] * multiplier;
        }
    };

    auto update_univariate = [&](std::size_t const arg_id, std::size_t const res_id, Float const der_local_1) {
        std::size_t& res_pos = node_location_on_buffer[res_id];
        std::size_t& arg_pos = node_location_on_buffer[arg_id];

        bool const arg_inplace = (arg_pos == passive_id<std::size_t>);

        if (arg_inplace) {
            mul_inplace(res_pos, der_local_1);
            // res id should now be arg id, avoiding a copy and a potential buffer increase
            std::swap(arg_pos, res_pos);
        }
        else {
            mul_add(res_pos, arg_pos, der_local_1);
            this->buffer.free_loc(res_pos);
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
                    this->buffer.free_loc(pos);
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

                    bool const arg_inplace = (arg_pos == passive_id<std::size_t>);
                    if constexpr (Reset) {
                        if (arg_inplace) {
                            // res id should now be lhs id, avoiding a copy and a potential buffer increase
                            std::swap(arg_pos, res_pos);
                        }
                        else {
                            add(res_pos, arg_pos);
                            this->buffer.free_loc(res_pos);
                            res_pos = passive_id<std::size_t>;
                        }
                    }
                    else {
                        if (arg_inplace) {
                            arg_pos = this->buffer.get_new_loc();
                            copy(res_pos, arg_pos);
                        }
                        else {
                            add(res_pos, arg_pos);
                        }
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
                    bool const rhs_is_new = (rhs_pos == passive_id<std::size_t>);
                    bool const rhs_inplace = !lhs_inplace && rhs_is_new;

                    if (!lhs_inplace) {
                        add(res_pos, lhs_pos);
                    }

                    if (!rhs_inplace) {
                        if (rhs_is_new) {
                            rhs_pos = this->buffer.get_new_loc();
                            copy(res_pos, rhs_pos);
                        }
                        else {
                            add(res_pos, rhs_pos);
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
                        this->buffer.free_loc(res_pos);
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
                    bool const rhs_is_new = (rhs_pos == passive_id<std::size_t>);
                    bool const rhs_inplace = !lhs_inplace && rhs_is_new;

                    if (!lhs_inplace) {
                        add(res_pos, lhs_pos);
                    }

                    if (!rhs_inplace) {
                        if (rhs_is_new) {
                            rhs_pos = this->buffer.get_new_loc();
                            copy_minus(res_pos, rhs_pos);
                        }
                        else {
                            sub(res_pos, rhs_pos);
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
                        this->buffer.free_loc(res_pos);
                        res_pos = passive_id<std::size_t>;
                    }
                }
                break;
            }
            case OpCode::MUL: {
                val_idx -= 2;
                id_idx -= 3;
                if (use_this_op) {
                    auto const lhs_val = static_cast<Float>(vals[val_idx]);
                    auto const rhs_val = static_cast<Float>(vals[val_idx + 1]);
                    std::size_t const lhs_id = ids[id_idx];
                    std::size_t const rhs_id = ids[id_idx + 1];
                    std::size_t const res_id = ids[id_idx + 2];

                    std::size_t& res_pos = node_location_on_buffer[res_id];
                    std::size_t& lhs_pos = node_location_on_buffer[lhs_id];
                    std::size_t& rhs_pos = node_location_on_buffer[rhs_id];

                    bool const lhs_inplace = (lhs_pos == passive_id<std::size_t>);
                    bool const rhs_is_new = (rhs_pos == passive_id<std::size_t>);
                    bool const rhs_inplace = !lhs_inplace && rhs_is_new;

                    if (!lhs_inplace) {
                        mul_add(res_pos, lhs_pos, rhs_val);
                    }

                    if (!rhs_inplace) {
                        if (rhs_is_new) {
                            // this is a new value, we NEED to override
                            rhs_pos = this->buffer.get_new_loc();
                            mul_set(res_pos, rhs_pos, lhs_val);
                        }
                        else {
                            mul_add(res_pos, rhs_pos, lhs_val);
                        }
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
                        this->buffer.free_loc(res_pos);
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
                        add(res_pos, arg_pos);
                        this->buffer.free_loc(res_pos);
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
                        sub(res_pos, arg_pos);
                        this->buffer.free_loc(res_pos);
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
                    auto const der_local_1 = static_cast<Float>(vals[val_idx]);
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
                    auto const val = static_cast<Float>(vals[val_idx]);
                    Float const der_local_1 = Float{ 2.0 } * val;
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
                    auto const val = static_cast<Float>(vals[val_idx]);
                    Float const der_local_1 = -val * val;
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
                    auto const val = static_cast<Float>(vals[val_idx]);
                    Float const der_local_1 = std::copysign(Float{ 1.0 }, val);
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
                    auto const der_local_1 = static_cast<Float>(vals[val_idx]);
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
                    auto const val = static_cast<Float>(vals[val_idx]);
                    Float const der_local_1 = Float{ 1.0 } / val;
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
                    constexpr Float two_over_root_pi = Float{ 2. } * std::numbers::inv_sqrtpi_v<Float>;
                    auto const val = static_cast<Float>(vals[val_idx]);
                    Float const der_local_1 = std::exp(-val * val) * two_over_root_pi;
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
                    constexpr Float minus_two_over_root_pi = Float{ -2. } * std::numbers::inv_sqrtpi_v<Float>;
                    auto const val = static_cast<Float>(vals[val_idx]);
                    Float const der_local_1 = std::exp(-val * val) * minus_two_over_root_pi;
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
                    auto const val = static_cast<Float>(vals[val_idx]);
                    Float const der_local_1 = -std::sin(val);
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
                    auto const val1 = static_cast<Float>(vals[val_idx]);
                    auto const val2 = static_cast<Float>(vals[val_idx + 1]);
                    Float const one_over_in = Float{ 1. } / val1;
                    Float const der_local_1 = Float{ 0.5 } * val2 * one_over_in;
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
                    auto const lhs_arg = static_cast<Float>(vals[val_idx]);
                    auto const rhs_arg = static_cast<Float>(vals[val_idx + 1]);
                    Float const der_local_1 =
                      rhs_arg != Float{ 0.0 } ? rhs_arg * std::pow(lhs_arg, rhs_arg - Float{ 1. }) : Float{ 0.0 };
                    update_univariate(arg_id, res_id, der_local_1);
                }
                break;
            }
        }
    }

    if constexpr (Reset) {
        this->node_location_on_buffer.resize(pos.op_position);
    }
}

} // namespace adhoc

#endif // ADHOC_BACKPROPAGATOR1LOSSY_HPP
