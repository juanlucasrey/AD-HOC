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

#ifndef ADHOC_BACKPROPAGATOR1_HPP
#define ADHOC_BACKPROPAGATOR1_HPP

#include "position_impl.hpp"
#include "tape_data.hpp"
#include "vector_size_of.hpp"

#include <cmath>
#include <concepts>
#include <iostream>
#include <numbers>
#include <span>
#include <vector>

namespace adhoc {

template<std::floating_point Float, bool Vectorised = false>
class BackPropagator {
  private:
    std::size_t m_num_lanes{ 1 };
    std::vector<Float> derivatives;
    std::vector<char> use_op;

  public:
    explicit BackPropagator() = default;

    void set_checkpoint(std::size_t /* ops_size */) {}
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
    void reserve_input(std::size_t /* count_registered */) {}
    void reserve_output(std::size_t /* count_registered */) {}
    void register_variable(std::size_t /* var_id */) {}
    void register_output_variable(std::size_t /* var_id */, std::size_t /* ops_size */) {}

    void set_derivative(std::size_t var_id, double deriv, std::size_t ops_size, std::size_t lane = 0)
    {
        if (lane < this->m_num_lanes) {
            this->derivatives.resize(ops_size * this->m_num_lanes);
            this->derivatives[var_id * this->m_num_lanes + lane] = static_cast<Float>(deriv);
        }
        else {
            throw;
        }
    }

    void set_derivative(std::size_t /* var_id1 */, std::size_t /* var_id2 */, double /* deriv */) { throw; }

    auto get_derivative(std::size_t var_id, std::size_t lane) const -> double
    {
        if (lane < this->m_num_lanes) {
            return static_cast<double>(this->derivatives[var_id * this->m_num_lanes + lane]);
        }

        throw;
        return 0.;
    }

    auto get_derivative(std::size_t /* var_id1 */, std::size_t /* var_id2 */, std::size_t /* lane */) const -> double
    {
        return 0.;
    }

    void clear()
    {
        this->derivatives.clear();
        this->use_op.clear();
    }

    void zero_adjoints() { std::fill(this->derivatives.begin(), this->derivatives.end(), 0.0); }

    auto size_of(bool capacity = false) const -> std::size_t
    {
        std::size_t size = 0;
        size += sizeof(std::size_t); // m_num_lanes
        size += vector_size_of(this->derivatives, capacity);
        size += vector_size_of(this->use_op, capacity);
        return size;
    }

    void reset(PositionImpl const& pos)
    {
        std::size_t to = pos.op_position;
        this->derivatives.resize(to * this->m_num_lanes);
        this->use_op.resize(to);
    }

    template<bool Reset, class TapeDataType>
    void backpropagate_to(PositionImpl const& pos, TapeDataType const& data);
};

template<std::floating_point Float, bool Vectorised>
template<bool Reset, class TapeDataType>
void
BackPropagator<Float, Vectorised>::backpropagate_to(PositionImpl const& pos, TapeDataType const& data)
{
    std::size_t to = pos.op_position;
    std::size_t from = data.next_id;

    const auto& ops = data.ops;
    const auto& vals = data.vals;
    const auto& ids = data.ids;

    this->derivatives.resize(ops.size() * this->m_num_lanes);
    this->use_op.resize(ops.size());

    auto add_derivative = [&](std::size_t arg_id, std::size_t res_id) {
        this->use_op[arg_id] = true;
        if constexpr (Vectorised) {
            auto const src = std::span<Float const>(&this->derivatives[res_id * this->m_num_lanes], this->m_num_lanes);
            auto dest = std::span<Float>(&this->derivatives[arg_id * this->m_num_lanes], this->m_num_lanes);
#pragma omp simd
            for (std::size_t i = 0; i < dest.size(); ++i) {
                dest[i] += src[i];
            }
        }
        else {
            this->derivatives[arg_id] += this->derivatives[res_id];
        }
    };

    auto sub_derivative = [&](std::size_t arg_id, std::size_t res_id) {
        this->use_op[arg_id] = true;
        if constexpr (Vectorised) {
            auto const src = std::span<Float const>(&this->derivatives[res_id * this->m_num_lanes], this->m_num_lanes);
            auto dest = std::span<Float>(&this->derivatives[arg_id * this->m_num_lanes], this->m_num_lanes);
#pragma omp simd
            for (std::size_t i = 0; i < dest.size(); ++i) {
                dest[i] -= src[i];
            }
        }
        else {
            this->derivatives[arg_id] -= this->derivatives[res_id];
        }
    };

    auto mul_add_derivative = [&](std::size_t arg_id, std::size_t res_id, Float multiplier) {
        this->use_op[arg_id] = true;
        if constexpr (Vectorised) {
            auto const src = std::span<Float const>(&this->derivatives[res_id * this->m_num_lanes], this->m_num_lanes);
            auto dest = std::span<Float>(&this->derivatives[arg_id * this->m_num_lanes], this->m_num_lanes);
#pragma omp simd
            for (std::size_t i = 0; i < dest.size(); ++i) {
                dest[i] += src[i] * multiplier;
            }
        }
        else {
            this->derivatives[arg_id] += this->derivatives[res_id] * multiplier;
        }
    };

    std::size_t val_idx = vals.size();
    std::size_t id_idx = ids.size();
    for (std::size_t op_idx = from; op_idx-- > to;) {
        OpCode const& op = ops[op_idx];
        bool const use_this_op = this->use_op[op_idx];

        switch (op) {
            case OpCode::REG_INPUT: {
                id_idx -= 1;
                break;
            }
            case OpCode::REG_OUTPUT: {
                id_idx -= 2;
                std::size_t const arg_id = ids[id_idx];
                std::size_t const res_id = ids[id_idx + 1];
                add_derivative(arg_id, res_id);
                break;
            }
            case OpCode::ADD: {
                id_idx -= 3;
                if (use_this_op) {
                    std::size_t const lhs_id = ids[id_idx];
                    std::size_t const rhs_id = ids[id_idx + 1];
                    std::size_t const res_id = ids[id_idx + 2];
                    add_derivative(lhs_id, res_id);
                    add_derivative(rhs_id, res_id);
                }
                break;
            }
            case OpCode::SUB: {
                id_idx -= 3;
                if (use_this_op) {
                    std::size_t const lhs_id = ids[id_idx];
                    std::size_t const rhs_id = ids[id_idx + 1];
                    std::size_t const res_id = ids[id_idx + 2];
                    add_derivative(lhs_id, res_id);
                    sub_derivative(rhs_id, res_id);
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
                    mul_add_derivative(lhs_id, res_id, rhs_val);
                    mul_add_derivative(rhs_id, res_id, lhs_val);
                }
                break;
            }
            case OpCode::ADD_C: {
                id_idx -= 2;
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    add_derivative(arg_id, res_id);
                }
                break;
            }
            case OpCode::SUB_C: {
                id_idx -= 2;
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    sub_derivative(arg_id, res_id);
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
                    mul_add_derivative(arg_id, res_id, der_local_1);
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
                    mul_add_derivative(arg_id, res_id, der_local_1);
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
                    mul_add_derivative(arg_id, res_id, der_local_1);
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
                    mul_add_derivative(arg_id, res_id, der_local_1);
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
                    mul_add_derivative(arg_id, res_id, der_local_1);
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
                    mul_add_derivative(arg_id, res_id, der_local_1);
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
                    mul_add_derivative(arg_id, res_id, der_local_1);
                }
                break;
            }
            case OpCode::ERFC: {
                val_idx -= 1;
                id_idx -= 2;
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    constexpr Float minus_two_over_root_pi = -Float{ 2. } * std::numbers::inv_sqrtpi_v<Float>;
                    auto const val = static_cast<Float>(vals[val_idx]);
                    Float const der_local_1 = std::exp(-val * val) * minus_two_over_root_pi;
                    mul_add_derivative(arg_id, res_id, der_local_1);
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
                    mul_add_derivative(arg_id, res_id, der_local_1);
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
                    Float const one_over_in = Float{ 1.0 } / val1;
                    Float const der_local_1 = Float{ 0.5 } * val2 * one_over_in;
                    mul_add_derivative(arg_id, res_id, der_local_1);
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
                    mul_add_derivative(arg_id, res_id, der_local_1);
                }
                break;
            }
        }
    }

    // we'll see later what system drives this boolean
    constexpr bool Log = false;
    if constexpr (Log) {
        // calculate the percentage of operations saved using the vector use_op
        std::size_t avoided_ops = 0;
        for (std::size_t i = 0; i < this->use_op.size(); ++i) {
            avoided_ops += !this->use_op[i];
        }

        double percentage = static_cast<double>(avoided_ops) / static_cast<double>(this->use_op.size());
        std::cout.precision(std::numeric_limits<double>::max_digits10);
        std::cout << "backpropagation saved operations: " << percentage << std::endl;
    }
}

} // namespace adhoc

#endif // ADHOC_BACKPROPAGATOR1_HPP
