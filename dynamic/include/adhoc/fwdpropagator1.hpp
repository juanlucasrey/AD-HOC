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

#ifndef ADHOC_FWDPROPAGATOR1_HPP
#define ADHOC_FWDPROPAGATOR1_HPP

#include "position_impl.hpp"
#include "tape_data.hpp"
#include "vector_size_of.hpp"

#include <algorithm>
#include <cmath>
#include <concepts>
#include <numbers>
#include <span>
#include <vector>

namespace adhoc {

template<std::floating_point Float, bool Vectorised = false>
class FwdPropagator {
  private:
    std::size_t m_num_lanes{ 1 };
    std::size_t m_num_inputs{ 0 };
    std::size_t m_num_outputs{ 0 };
    std::vector<Float> derivatives;

    PositionImpl pos_local;

    std::vector<std::size_t> input_ids;
    // std::vector<Float> input_derivatives;

    std::vector<std::size_t> output_ids;
    std::vector<Float> output_derivatives;

  public:
    explicit FwdPropagator(std::size_t n_inputs, std::size_t n_outputs, std::size_t num_lanes)
      : m_num_lanes(num_lanes)
      , m_num_inputs(n_inputs)
      , m_num_outputs(n_outputs)
    {
        // input_derivatives.resize(n_inputs * num_lanes, 0.0);
        // // fill with one seed per input
        // for (std::size_t i = 0; i < n_inputs; ++i) {
        //     input_derivatives[(i * num_lanes) + i] = 1.0;
        // }

        output_derivatives.resize(n_outputs * num_lanes, 0.0);
    }

    void set_checkpoint(std::size_t /* ops_size */) {}
    void set_lanes(std::size_t /* num_lanes */) {}
    auto get_lanes() const -> std::size_t { return this->m_num_lanes; }
    void reserve_input(std::size_t /* count_registered */) {}
    void reserve_output(std::size_t /* count_registered */) {}
    void register_variable(std::size_t var_id) { input_ids.push_back(var_id); }
    void register_output_variable(std::size_t var_id, std::size_t /* ops_size */) { output_ids.push_back(var_id); }

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
            auto const it = std::find(output_ids.begin(), output_ids.end(), var_id);
            if (it != output_ids.end()) {
                auto const distance = std::distance(output_ids.begin(), it);
                return static_cast<double>(output_derivatives[(distance * this->m_num_lanes) + lane]);
            }
            else {
                throw;
            }
        }

        throw;
        return 0.;
    }

    auto get_derivative2(std::size_t idx_input, std::size_t idx_output) const -> double
    {
        if (idx_output >= m_num_outputs) {
            throw;
        }

        if (idx_input >= this->m_num_lanes) {
            throw;
        }

        return static_cast<double>(output_derivatives[(idx_output * this->m_num_lanes) + idx_input]);
    }

    auto get_derivative(std::size_t /* var_id1 */, std::size_t /* var_id2 */, std::size_t /* lane */) const -> double
    {
        return 0.;
    }

    void clear() { this->derivatives.clear(); }

    void zero_adjoints() { std::fill(this->derivatives.begin(), this->derivatives.end(), 0.0); }

    auto size_of(bool capacity = false) const -> std::size_t
    {
        std::size_t size = 0;
        size += sizeof(std::size_t); // m_num_lanes
        size += vector_size_of(this->derivatives, capacity);
        return size;
    }

    template<bool Reset, bool ResetInPlace, bool Log, class TapeDataType>
    void backpropagate_to(PositionImpl const& pos, TapeDataType& data);
};

template<std::floating_point Float, bool Vectorised>
template<bool Reset, bool ResetInPlace, bool Log, class TapeDataType>
void
FwdPropagator<Float, Vectorised>::backpropagate_to(PositionImpl const& pos, TapeDataType& data)
{
    std::size_t to = pos.op_position;
    std::size_t from = data.next_id;

    const auto& ops = data.ops;
    const auto& vals = data.vals;
    const auto& ids = data.ids;

    this->derivatives.resize(ops.size() * this->m_num_lanes);

    auto add_derivative = [&](std::size_t arg_id, std::size_t res_id) {
        if constexpr (Vectorised) {
            auto const src = std::span<Float const>(&this->derivatives[arg_id * this->m_num_lanes], this->m_num_lanes);
            auto dest = std::span<Float>(&this->derivatives[res_id * this->m_num_lanes], this->m_num_lanes);
#pragma omp simd
            for (std::size_t i = 0; i < dest.size(); ++i) {
                dest[i] += src[i];
            }
        }
        else {
            this->derivatives[res_id] += this->derivatives[arg_id];
        }
    };

    auto add_derivative_output = [&](std::size_t src_id, std::size_t dest_id) {
        if constexpr (Vectorised) {
            auto const src = std::span<Float const>(&this->derivatives[src_id * this->m_num_lanes], this->m_num_lanes);
            auto dest = std::span<Float>(&output_derivatives[dest_id * this->m_num_lanes], this->m_num_lanes);
#pragma omp simd
            for (std::size_t i = 0; i < src.size(); ++i) {
                dest[i] += src[i];
            }
        }
        else {
            output_derivatives[dest_id] += this->derivatives[src_id];
        }
    };

    //     auto add_derivative_input = [&](std::size_t src_id, std::size_t dest_id) {
    //         if constexpr (Vectorised) {
    //             auto const src = std::span<Float const>(&input_derivatives[src_id * this->m_num_lanes],
    //             this->m_num_lanes); auto dest = std::span<Float>(&this->derivatives[dest_id * this->m_num_lanes],
    //             this->m_num_lanes);
    // #pragma omp simd
    //             for (std::size_t i = 0; i < dest.size(); ++i) {
    //                 dest[i] += src[i];
    //             }
    //         }
    //         else {
    //             this->derivatives[dest_id] += input_derivatives[src_id];
    //         }
    //     };

    auto sub_derivative = [&](std::size_t arg_id, std::size_t res_id) {
        if constexpr (Vectorised) {
            auto const src = std::span<Float const>(&this->derivatives[arg_id * this->m_num_lanes], this->m_num_lanes);
            auto dest = std::span<Float>(&this->derivatives[res_id * this->m_num_lanes], this->m_num_lanes);
#pragma omp simd
            for (std::size_t i = 0; i < dest.size(); ++i) {
                dest[i] -= src[i];
            }
        }
        else {
            this->derivatives[res_id] -= this->derivatives[arg_id];
        }
    };

    auto mul_add_derivative = [&](std::size_t arg_id, std::size_t res_id, Float multiplier) {
        if constexpr (Vectorised) {
            auto const src = std::span<Float const>(&this->derivatives[arg_id * this->m_num_lanes], this->m_num_lanes);
            auto dest = std::span<Float>(&this->derivatives[res_id * this->m_num_lanes], this->m_num_lanes);
#pragma omp simd
            for (std::size_t i = 0; i < dest.size(); ++i) {
                dest[i] += src[i] * multiplier;
            }
        }
        else {
            this->derivatives[res_id] += this->derivatives[arg_id] * multiplier;
        }
    };

    for (; pos_local.op_position < from; ++pos_local.op_position) {
        OpCode const& op = ops[pos_local.op_position];

        switch (op) {
            case OpCode::REG_INPUT: {
                std::size_t const arg_id = ids[pos_local.id_position];

                auto const it = std::find(input_ids.begin(), input_ids.end(), arg_id);
                if (it == input_ids.end()) {
                    throw;
                }
                // auto const distance = static_cast<std::size_t>(std::distance(input_ids.begin(), it));

                // add_derivative_input(distance, arg_id);
                pos_local.id_position += 1;
                break;
            }
            case OpCode::REG_OUTPUT: {
                std::size_t const arg_id = ids[pos_local.id_position];
                std::size_t const res_id = ids[pos_local.id_position + 1];
                add_derivative(arg_id, res_id);

                auto const it = std::find(output_ids.begin(), output_ids.end(), res_id);
                if (it == output_ids.end()) {
                    throw;
                }
                auto const distance = static_cast<std::size_t>(std::distance(output_ids.begin(), it));
                add_derivative_output(res_id, distance);
                pos_local.id_position += 2;
                break;
            }
            case OpCode::ADD: {
                std::size_t const lhs_id = ids[pos_local.id_position];
                std::size_t const rhs_id = ids[pos_local.id_position + 1];
                std::size_t const res_id = ids[pos_local.id_position + 2];
                add_derivative(lhs_id, res_id);
                add_derivative(rhs_id, res_id);
                pos_local.id_position += 3;
                break;
            }
            case OpCode::SUB: {
                std::size_t const lhs_id = ids[pos_local.id_position];
                std::size_t const rhs_id = ids[pos_local.id_position + 1];
                std::size_t const res_id = ids[pos_local.id_position + 2];
                add_derivative(lhs_id, res_id);
                sub_derivative(rhs_id, res_id);
                pos_local.id_position += 3;
                break;
            }
            case OpCode::MUL: {
                auto const lhs_val = static_cast<Float>(vals[pos_local.val_position]);
                auto const rhs_val = static_cast<Float>(vals[pos_local.val_position + 1]);
                std::size_t const lhs_id = ids[pos_local.id_position];
                std::size_t const rhs_id = ids[pos_local.id_position + 1];
                std::size_t const res_id = ids[pos_local.id_position + 2];
                mul_add_derivative(lhs_id, res_id, rhs_val);
                mul_add_derivative(rhs_id, res_id, lhs_val);
                pos_local.val_position += 2;
                pos_local.id_position += 3;
                break;
            }
            case OpCode::ADD_C: {
                std::size_t const arg_id = ids[pos_local.id_position];
                std::size_t const res_id = ids[pos_local.id_position + 1];
                add_derivative(arg_id, res_id);
                pos_local.id_position += 2;
                break;
            }
            case OpCode::SUB_C: {
                std::size_t const arg_id = ids[pos_local.id_position];
                std::size_t const res_id = ids[pos_local.id_position + 1];
                sub_derivative(arg_id, res_id);
                pos_local.id_position += 2;
                break;
            }
            case OpCode::MUL_C: {
                std::size_t const arg_id = ids[pos_local.id_position];
                std::size_t const res_id = ids[pos_local.id_position + 1];
                auto const der_local_1 = static_cast<Float>(vals[pos_local.val_position]);
                mul_add_derivative(arg_id, res_id, der_local_1);
                pos_local.val_position += 1;
                pos_local.id_position += 2;
                break;
            }
            case OpCode::NORM: {
                std::size_t const arg_id = ids[pos_local.id_position];
                std::size_t const res_id = ids[pos_local.id_position + 1];
                auto const val = static_cast<Float>(vals[pos_local.val_position]);
                Float const der_local_1 = Float{ 2.0 } * val;
                mul_add_derivative(arg_id, res_id, der_local_1);
                pos_local.val_position += 1;
                pos_local.id_position += 2;
                break;
            }
            case OpCode::INV: {
                std::size_t const arg_id = ids[pos_local.id_position];
                std::size_t const res_id = ids[pos_local.id_position + 1];
                auto const val = static_cast<Float>(vals[pos_local.val_position]);
                Float const der_local_1 = -val * val;
                mul_add_derivative(arg_id, res_id, der_local_1);
                pos_local.val_position += 1;
                pos_local.id_position += 2;
                break;
            }
            case OpCode::ABS: {
                std::size_t const arg_id = ids[pos_local.id_position];
                std::size_t const res_id = ids[pos_local.id_position + 1];
                auto const val = static_cast<Float>(vals[pos_local.val_position]);
                Float const der_local_1 = std::copysign(Float{ 1.0 }, val);
                mul_add_derivative(arg_id, res_id, der_local_1);
                pos_local.val_position += 1;
                pos_local.id_position += 2;
                break;
            }
            case OpCode::EXP: {
                std::size_t const arg_id = ids[pos_local.id_position];
                std::size_t const res_id = ids[pos_local.id_position + 1];
                auto const der_local_1 = static_cast<Float>(vals[pos_local.val_position]);
                mul_add_derivative(arg_id, res_id, der_local_1);
                pos_local.val_position += 1;
                pos_local.id_position += 2;
                break;
            }
            case OpCode::LOG: {
                std::size_t const arg_id = ids[pos_local.id_position];
                std::size_t const res_id = ids[pos_local.id_position + 1];
                auto const val = static_cast<Float>(vals[pos_local.val_position]);
                Float const der_local_1 = Float{ 1.0 } / val;
                mul_add_derivative(arg_id, res_id, der_local_1);
                pos_local.val_position += 1;
                pos_local.id_position += 2;
                break;
            }
            case OpCode::ERF: {
                std::size_t const arg_id = ids[pos_local.id_position];
                std::size_t const res_id = ids[pos_local.id_position + 1];
                constexpr Float two_over_root_pi = Float{ 2. } * std::numbers::inv_sqrtpi_v<Float>;
                auto const val = static_cast<Float>(vals[pos_local.val_position]);
                Float const der_local_1 = std::exp(-val * val) * two_over_root_pi;
                mul_add_derivative(arg_id, res_id, der_local_1);
                pos_local.val_position += 1;
                pos_local.id_position += 2;
                break;
            }
            case OpCode::ERFC: {
                std::size_t const arg_id = ids[pos_local.id_position];
                std::size_t const res_id = ids[pos_local.id_position + 1];
                constexpr Float minus_two_over_root_pi = -Float{ 2. } * std::numbers::inv_sqrtpi_v<Float>;
                auto const val = static_cast<Float>(vals[pos_local.val_position]);
                Float const der_local_1 = std::exp(-val * val) * minus_two_over_root_pi;
                mul_add_derivative(arg_id, res_id, der_local_1);
                pos_local.val_position += 1;
                pos_local.id_position += 2;
                break;
            }
            case OpCode::COS: {
                std::size_t const arg_id = ids[pos_local.id_position];
                std::size_t const res_id = ids[pos_local.id_position + 1];
                auto const val = static_cast<Float>(vals[pos_local.val_position]);
                Float const der_local_1 = -std::sin(val);
                mul_add_derivative(arg_id, res_id, der_local_1);
                pos_local.val_position += 2;
                pos_local.id_position += 2;
                break;
            }
            case OpCode::SQRT: {

                std::size_t const arg_id = ids[pos_local.id_position];
                std::size_t const res_id = ids[pos_local.id_position + 1];
                auto const val1 = static_cast<Float>(vals[pos_local.val_position]);
                auto const val2 = static_cast<Float>(vals[pos_local.val_position + 1]);
                Float const one_over_in = Float{ 1.0 } / val1;
                Float const der_local_1 = Float{ 0.5 } * val2 * one_over_in;
                mul_add_derivative(arg_id, res_id, der_local_1);
                pos_local.val_position += 2;
                pos_local.id_position += 2;
                break;
            }
            case OpCode::POW_C: {
                std::size_t const arg_id = ids[pos_local.id_position];
                std::size_t const res_id = ids[pos_local.id_position + 1];
                auto const lhs_arg = static_cast<Float>(vals[pos_local.val_position]);
                auto const rhs_arg = static_cast<Float>(vals[pos_local.val_position + 1]);
                Float const der_local_1 =
                  rhs_arg != Float{ 0.0 } ? rhs_arg * std::pow(lhs_arg, rhs_arg - Float{ 1. }) : Float{ 0.0 };
                mul_add_derivative(arg_id, res_id, der_local_1);
                pos_local.val_position += 2;
                pos_local.id_position += 2;
                break;
            }
        }
    }

    if constexpr (Reset) {
        this->derivatives.resize(to * this->m_num_lanes);
        pos_local = pos;
        reset(pos, data);
        output_ids.clear();
    }
}

} // namespace adhoc

#endif // ADHOC_FWDPROPAGATOR1_HPP
