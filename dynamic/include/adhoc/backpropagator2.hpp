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

#ifndef ADHOC_BACKPROPAGATOR2_HPP
#define ADHOC_BACKPROPAGATOR2_HPP

#include "maps.hpp"
#include "passive_id.hpp"
#include "position_impl.hpp"
#include "tape_data.hpp"

#include <array>
#include <cmath>
#include <concepts>
#include <numbers>
#include <vector>

namespace adhoc {

template<std::floating_point Float, MapType maptype, bool Vectorised = false>
class BackPropagator2 {
    static constexpr std::size_t SIMD_WIDTH = 8;

  private:
    // this works with vector instead of array, but its 3x slower!!

    // using der_container_t = std::conditional_t<Vectorised, std::vector<double>, double>;
    using der_container_t = std::conditional_t<Vectorised, std::array<double, 8>, double>;

    std::size_t m_num_lanes{ 1 };
    std::vector<map_t<maptype, std::size_t, der_container_t> > derivatives;
    std::vector<char> use_op;

  public:
    explicit BackPropagator2() = default;

    void set_checkpoint(std::size_t /* ops_size */) {}
    void set_lanes(std::size_t num_lanes)
    {
        if constexpr (!Vectorised) {
            if (num_lanes != 1) {
                // This backpropagator is not designed for multiple lanes
                throw;
            }
        }
        else {
            if (num_lanes != SIMD_WIDTH) {
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
            this->derivatives.resize(ops_size);

            auto& derivatives_local = this->derivatives[var_id][passive_id<std::size_t>];
            if constexpr (Vectorised) {
                // is using vector this needs to be uncommented
                // derivatives_local.resize(this->m_num_lanes, 0.0);
                derivatives_local[lane] = deriv;
            }
            else {
                derivatives_local = deriv;
            }
        }
        else {
            throw;
        }
    }

    void set_derivative(std::size_t var_id1, std::size_t var_id2, double deriv, std::size_t lane = 0)
    {
        if (lane < this->m_num_lanes) {
            // this->derivatives.resize(ops_size);

            if (var_id1 < var_id2 && var_id2 != passive_id<std::size_t>) {
                std::swap(var_id1, var_id2);
            }

            auto& derivatives_local = this->derivatives[var_id1][var_id2];
            if constexpr (Vectorised) {
                // is using vector this needs to be uncommented
                // derivatives_local.resize(this->m_num_lanes, 0.0);
                derivatives_local[lane] = deriv;
            }
            else {
                derivatives_local = deriv;
            }
        }
        else {
            throw;
        }
    }

    auto get_derivative(std::size_t var_id, std::size_t lane) const -> double
    {
        if (lane < this->m_num_lanes) {
            auto const& idmap = this->derivatives[var_id];
            auto const it = idmap.find(passive_id<std::size_t>);
            if (it != idmap.end()) {
                if constexpr (Vectorised) {
                    return it->second[lane];
                }
                else {
                    return it->second;
                }
            }
            return 0.;
        }

        throw;
        return 0.;
    }

    auto get_derivative(std::size_t var_id1, std::size_t var_id2, std::size_t lane) const -> double
    {
        if (lane < this->m_num_lanes) {
            if (var_id1 < var_id2) {
                std::swap(var_id1, var_id2);
            }

            auto const& idmap = this->derivatives[var_id1];
            auto const it = idmap.find(var_id2);
            if (it != idmap.end()) {
                if constexpr (Vectorised) {
                    return it->second[lane];
                }
                else {
                    return it->second;
                }
            }

            return 0.;
        }

        throw;
        return 0.;
    }

    void clear()
    {
        this->derivatives.clear();
        this->use_op.clear();
    }

    void zero_adjoints() { this->derivatives.clear(); }

    auto size_of(bool /* capacity */ = false) const -> std::size_t
    {
        std::size_t size = 0;
        // TODO
        return size;
    }

    template<bool Reset, bool ResetInPlace, class TapeDataType>
    void backpropagate_to(PositionImpl const& pos, TapeDataType& data);
};

template<std::floating_point Float, MapType maptype, bool Vectorised>
template<bool Reset, bool ResetInPlace, class TapeDataType>
void
BackPropagator2<Float, maptype, Vectorised>::backpropagate_to(PositionImpl const& pos, TapeDataType& data)
{
    std::size_t to = pos.op_position;
    std::size_t from = data.next_id;

    const auto& ops = data.ops;
    const auto& vals = data.vals;
    const auto& ids = data.ids;

    this->derivatives.resize(ops.size());
    this->use_op.resize(ops.size());

    auto add_derivative = [&](std::size_t id1, std::size_t id2, auto const& value) {
        if (id1 < id2 && id2 != passive_id<std::size_t>) {
            std::swap(id1, id2);
        }

        if constexpr (Vectorised) {
            auto& dest_vec = this->derivatives[id1][id2];

            // is using vector this needs to be uncommented
            // if (dest_vec.empty()) {
            //     dest_vec.resize(this->m_num_lanes, 0.0);
            // }

            double* dest = dest_vec.data();
            const double* src = value.data();
#pragma omp simd
            for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                dest[i] += src[i];
            }
        }
        else {
            this->derivatives[id1][id2] += value;
        }
    };

    auto sub_derivative = [&](std::size_t id1, std::size_t id2, auto const& value) {
        if (id1 < id2 && id2 != passive_id<std::size_t>) {
            std::swap(id1, id2);
        }

        if constexpr (Vectorised) {
            auto& dest_vec = this->derivatives[id1][id2];

            // is using vector this needs to be uncommented
            // if (dest_vec.empty()) {
            //     dest_vec.resize(this->m_num_lanes, 0.0);
            // }
            double* dest = dest_vec.data();
            const double* src = value.data();
#pragma omp simd
            for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                dest[i] -= src[i];
            }
        }
        else {
            this->derivatives[id1][id2] -= value;
        }
    };

    auto mul_add_derivative = [&](std::size_t id1, std::size_t id2, auto const& value, double scale) {
        if (id1 < id2 && id2 != passive_id<std::size_t>) {
            std::swap(id1, id2);
        }
        if constexpr (Vectorised) {
            auto& dest_vec = this->derivatives[id1][id2];

            // is using vector this needs to be uncommented
            // if (dest_vec.empty()) {
            //     dest_vec.resize(this->m_num_lanes, 0.0);
            // }
            double* dest = dest_vec.data();
            const double* src = value.data();
#pragma omp simd
            for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                dest[i] += src[i] * scale;
            }
        }
        else {
            this->derivatives[id1][id2] += value * scale;
        }
    };

    auto update_univariate = [&]<bool HasSecondDerivative = true>(
                               std::size_t arg_id, std::size_t res_id, double der_local_1, double der_local_2 = 0.) {
        this->use_op[arg_id] = true;
        auto const& der_list = this->derivatives[res_id];
        for (auto const& der_pair : der_list) {
            std::size_t const der_id = der_pair.first;
            auto const& der_value = der_pair.second;
            if (der_id == passive_id<std::size_t>) {
                mul_add_derivative(arg_id, der_id, der_value, der_local_1);
                if constexpr (HasSecondDerivative) {
                    mul_add_derivative(arg_id, arg_id, der_value, der_local_2);
                }
            }
            else if (der_id == res_id) {
                mul_add_derivative(arg_id, arg_id, der_value, der_local_1 * der_local_1);
            }
            else {
                mul_add_derivative(arg_id, der_id, der_value, (arg_id == der_id ? 2. : 1.) * der_local_1);
            }
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

                this->use_op[arg_id] = true;

                auto const& der_list = this->derivatives[res_id];
                for (auto const& der_pair : der_list) {
                    add_derivative(arg_id, der_pair.first, der_pair.second);
                }
                break;
            }
            case OpCode::ADD: {
                id_idx -= 3;
                if (use_this_op) {
                    std::size_t const lhs_id = ids[id_idx];
                    std::size_t const rhs_id = ids[id_idx + 1];
                    std::size_t const res_id = ids[id_idx + 2];

                    this->use_op[lhs_id] = true;
                    this->use_op[rhs_id] = true;

                    auto const& der_list = this->derivatives[res_id];
                    for (auto const& der_pair : der_list) {
                        std::size_t const der_id = der_pair.first;
                        auto const& der_value = der_pair.second;

                        if (der_id == passive_id<std::size_t>) {
                            add_derivative(lhs_id, der_id, der_value);
                            add_derivative(rhs_id, der_id, der_value);
                        }
                        else if (der_id == res_id) {
                            add_derivative(lhs_id, rhs_id, der_value);
                            add_derivative(lhs_id, lhs_id, der_value);
                            add_derivative(rhs_id, rhs_id, der_value);
                        }
                        else {
                            mul_add_derivative(lhs_id, der_id, der_value, (lhs_id == der_id ? 2. : 1.));
                            mul_add_derivative(rhs_id, der_id, der_value, (rhs_id == der_id ? 2. : 1.));
                        }
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

                    this->use_op[lhs_id] = true;
                    this->use_op[rhs_id] = true;

                    auto const& der_list = this->derivatives[res_id];
                    for (auto const& der_pair : der_list) {
                        std::size_t const der_id = der_pair.first;
                        auto const& der_value = der_pair.second;

                        if (der_id == passive_id<std::size_t>) {
                            add_derivative(lhs_id, der_id, der_value);
                            sub_derivative(rhs_id, der_id, der_value);
                        }
                        else if (der_id == res_id) {
                            sub_derivative(lhs_id, rhs_id, der_value);
                            add_derivative(lhs_id, lhs_id, der_value);
                            add_derivative(rhs_id, rhs_id, der_value);
                        }
                        else {
                            mul_add_derivative(lhs_id, der_id, der_value, (lhs_id == der_id ? 2. : 1.));
                            mul_add_derivative(rhs_id, der_id, der_value, -(rhs_id == der_id ? 2. : 1.));
                        }
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

                    this->use_op[lhs_id] = true;
                    this->use_op[rhs_id] = true;

                    auto const& der_list = this->derivatives[res_id];
                    for (auto const& der_pair : der_list) {
                        std::size_t const der_id = der_pair.first;
                        auto const& der_value = der_pair.second;

                        if (der_id == passive_id<std::size_t>) {
                            mul_add_derivative(lhs_id, der_id, der_value, rhs_val);
                            mul_add_derivative(rhs_id, der_id, der_value, lhs_val);
                            add_derivative(rhs_id, lhs_id, der_value);
                        }
                        else if (der_id == res_id) {
                            mul_add_derivative(rhs_id, lhs_id, der_value, rhs_val * lhs_val);
                            mul_add_derivative(lhs_id, lhs_id, der_value, rhs_val * rhs_val);
                            mul_add_derivative(rhs_id, rhs_id, der_value, lhs_val * lhs_val);
                        }
                        else {
                            mul_add_derivative(lhs_id, der_id, der_value, (lhs_id == der_id ? 2. : 1.) * rhs_val);
                            mul_add_derivative(rhs_id, der_id, der_value, (rhs_id == der_id ? 2. : 1.) * lhs_val);
                        }
                    }
                }
                break;
            }
            case OpCode::ADD_C: {
                id_idx -= 2;
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];

                    this->use_op[arg_id] = true;

                    auto const& der_list = this->derivatives[res_id];
                    for (auto const& der_pair : der_list) {
                        std::size_t const der_id = der_pair.first;
                        auto const& der_value = der_pair.second;
                        if (der_id == passive_id<std::size_t>) {
                            add_derivative(arg_id, der_id, der_value);
                        }
                        else if (der_id == res_id) {
                            add_derivative(arg_id, arg_id, der_value);
                        }
                        else {
                            mul_add_derivative(arg_id, der_id, der_value, (arg_id == der_id ? 2. : 1.));
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

                    this->use_op[arg_id] = true;

                    auto const& der_list = this->derivatives[res_id];
                    for (auto const& der_pair : der_list) {
                        std::size_t const der_id = der_pair.first;
                        auto const& der_value = der_pair.second;
                        if (der_id == passive_id<std::size_t>) {
                            sub_derivative(arg_id, der_id, der_value);
                        }
                        else if (der_id == res_id) {
                            add_derivative(arg_id, arg_id, der_value);
                        }
                        else {
                            mul_add_derivative(arg_id, der_id, der_value, (arg_id == der_id ? 2. : 1.) * -1.0);
                        }
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
                    update_univariate.template operator()<false>(arg_id, res_id, der_local_1);
                }
                break;
            }
            case OpCode::NORM: {
                val_idx -= 1;
                id_idx -= 2;
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    double const val = vals[val_idx];
                    double const der_local_1 = 2.0 * val;
                    double const der_local_2 = 2.0;
                    update_univariate(arg_id, res_id, der_local_1, der_local_2);
                }
                break;
            }
            case OpCode::INV: {
                val_idx -= 1;
                id_idx -= 2;
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    double const val = vals[val_idx];
                    double const der_local_1 = -val * val;
                    double const der_local_2 = -2.0 * der_local_1 * val;
                    update_univariate(arg_id, res_id, der_local_1, der_local_2);
                }
                break;
            }
            case OpCode::ABS: {
                val_idx -= 1;
                id_idx -= 2;
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    double const val = vals[val_idx];
                    double const der_local_1 = std::copysign(1.0, val);
                    update_univariate.template operator()<false>(arg_id, res_id, der_local_1);
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
                    double const der_local_2 = der_local_1;
                    update_univariate(arg_id, res_id, der_local_1, der_local_2);
                }
                break;
            }
            case OpCode::LOG: {
                val_idx -= 1;
                id_idx -= 2;
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    double const val = vals[val_idx];
                    double const der_local_1 = 1.0 / val;
                    double const der_local_2 = -der_local_1 * der_local_1;
                    update_univariate(arg_id, res_id, der_local_1, der_local_2);
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
                    double const val = vals[val_idx];
                    double const der_local_1 = std::exp(-val * val) * two_over_root_pi;
                    double const der_local_2 = -2. * val * der_local_1;
                    update_univariate(arg_id, res_id, der_local_1, der_local_2);
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
                    double const val = vals[val_idx];
                    double const der_local_1 = std::exp(-val * val) * minus_two_over_root_pi;
                    double const der_local_2 = -2. * val * der_local_1;
                    update_univariate(arg_id, res_id, der_local_1, der_local_2);
                }
                break;
            }
            case OpCode::COS: {
                val_idx -= 2;
                id_idx -= 2;
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    double const val1 = vals[val_idx];
                    double const val2 = vals[val_idx + 1];
                    double const der_local_1 = -std::sin(val1);
                    double const der_local_2 = -val2;
                    update_univariate(arg_id, res_id, der_local_1, der_local_2);
                }
                break;
            }
            case OpCode::SQRT: {
                val_idx -= 2;
                id_idx -= 2;
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    double const val1 = vals[val_idx];
                    double const val2 = vals[val_idx + 1];
                    double const one_over_in = 1. / val1;
                    double const der_local_1 = 0.5 * val2 * one_over_in;
                    double const der_local_2 = -0.5 * der_local_1 * one_over_in;
                    update_univariate(arg_id, res_id, der_local_1, der_local_2);
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
                    double const der_local_2 = rhs_arg != 0.0 && rhs_arg != 1.0
                                                 ? rhs_arg * (rhs_arg - 1.0) * std::pow(lhs_arg, rhs_arg - 2.)
                                                 : 0.0;
                    update_univariate(arg_id, res_id, der_local_1, der_local_2);
                }
                break;
            }
        }
        if constexpr (Reset && ResetInPlace) {
            this->derivatives.resize(this->derivatives.size() - 1);
            this->use_op.resize(this->use_op.size() - 1);
        }
    }

    if constexpr (Reset && !ResetInPlace) {
        this->derivatives.resize(to);
        this->use_op.resize(to);
    }

    if constexpr (Reset) {
        reset(pos, data);
    }
}

} // namespace adhoc

#endif // ADHOC_BACKPROPAGATOR2_HPP
