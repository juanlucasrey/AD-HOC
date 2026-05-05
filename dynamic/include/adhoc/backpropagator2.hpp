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

#include "passive_id.hpp"
#include "tape_data.hpp"

#include <cmath>
#include <map>
#include <numbers>
#include <vector>

namespace adhoc {

template<class Float>
class BackPropagator2 {
  private:
    std::size_t m_num_lanes{ 1 };
    std::vector<std::map<std::size_t, double> > derivatives;
    std::vector<bool> use_op;

    std::map<std::size_t, std::map<std::size_t, double> > derivatives_input_map;

  public:
    explicit BackPropagator2() = default;

    void set_checkpoint(std::size_t /* ops_size */) {}
    void set_lanes(std::size_t num_lanes)
    {
        if (this->m_num_lanes != 1) {
            // This backpropagator is not designed for multiple lanes
            throw;
        }
        this->m_num_lanes = num_lanes;
    }
    auto get_lanes() const -> std::size_t { return this->m_num_lanes; }
    void reserve_input(std::size_t /* count_registered */) {}
    void reserve_output(std::size_t /* count_registered */) {}
    void register_variable(std::size_t /* var_id */) {}
    void register_output_variable(std::size_t /* var_id */, std::size_t /* ops_size */) {}

    void set_derivative(std::size_t var_id, double deriv, std::size_t /* ops_size */, std::size_t /* lane */ = 0)
    {
        this->derivatives_input_map[var_id][passive_id<std::size_t>] = deriv;
    }

    void set_derivative(std::size_t var_id1, std::size_t var_id2, double deriv)
    {
        if (var_id1 < var_id2 && var_id2 != passive_id<std::size_t>) {
            std::swap(var_id1, var_id2);
        }
        this->derivatives_input_map[var_id1][var_id2] = deriv;
    }

    auto get_derivative(std::size_t var_id, std::size_t /* lane */) const -> double
    {
        auto const& idmap = this->derivatives[var_id];
        auto const it = idmap.find(passive_id<std::size_t>);
        if (it != idmap.end()) {
            return it->second;
        }

        return 0.;
    }

    auto get_derivative(std::size_t var_id1, std::size_t var_id2, std::size_t /* lane */) const -> double
    {
        if (var_id1 < var_id2) {
            std::swap(var_id1, var_id2);
        }

        auto const& idmap = this->derivatives[var_id1];
        auto const it = idmap.find(var_id2);
        if (it != idmap.end()) {
            return it->second;
        }

        return 0.;
    }

    // Clear all derivatives
    void clear()
    {
        this->derivatives_input_map.clear();
        this->derivatives.clear();
    }

    void zero_adjoints()
    {
        this->derivatives_input_map.clear();
        this->derivatives.clear();
    }

    auto size_of(bool /* capacity */ = false) const -> std::size_t { return 0; }

    template<bool Reset = false, bool ResetInPlace = false, bool Log = false>
    void backpropagate_to(std::size_t to, TapeData& data);
};

template<class Float>
template<bool Reset, bool ResetInPlace, bool Log>
void
BackPropagator2<Float>::backpropagate_to(std::size_t to, TapeData& data)
{
    std::size_t from = data.next_id;
    const auto& ops = data.ops;
    const auto& vals = data.vals;
    const auto& ids = data.ids;

    std::size_t val_idx = 0.;
    std::size_t id_idx = 0.;
    use_op.resize(ops.size());

    val_idx = vals.size();
    id_idx = ids.size();
    for (std::size_t op_idx = from; op_idx-- > to;) {
        OpCode const& op = ops[op_idx];
        switch (op) {
            case OpCode::REG_INPUT: {
                id_idx -= 1;

                use_op[op_idx] = true;
                break;
            }
            case OpCode::REG_OUTPUT: {
                id_idx -= 2;

                use_op[ids[id_idx]] = true;
                use_op[ids[id_idx + 1]] = true;

                break;
            }
            case OpCode::ADD:
            case OpCode::SUB:
            case OpCode::MUL: {
                id_idx -= 3;
                if (use_op[op_idx]) {
                    std::size_t const lhs_id = ids[id_idx];
                    std::size_t const rhs_id = ids[id_idx + 1];

                    use_op[lhs_id] = true;
                    use_op[rhs_id] = true;
                }

                break;
            }
            case OpCode::ADD_C:
            case OpCode::SUB_C:
            case OpCode::MUL_C:
            case OpCode::NORM:
            case OpCode::ABS:
            case OpCode::EXP:
            case OpCode::LOG:
            case OpCode::ERF:
            case OpCode::ERFC:
            case OpCode::COS:
            case OpCode::SQRT:
            case OpCode::POW_C:
            case OpCode::INV: {
                id_idx -= 2;
                if (use_op[op_idx]) {
                    std::size_t const arg_id = ids[id_idx];
                    use_op[arg_id] = true;
                }

                break;
            }
        }
    }

    this->derivatives.resize(ops.size());

    auto add_derivative = [&](std::size_t id1, std::size_t id2, double value) {
        if (id1 < id2 && id2 != passive_id<std::size_t>) {
            std::swap(id1, id2);
        }
        this->derivatives[id1][id2] += value;
    };

    val_idx = vals.size();
    id_idx = ids.size();
    for (std::size_t op_idx = from; op_idx-- > to;) {
        OpCode const& op = ops[op_idx];
        if (!static_cast<bool>(use_op[op_idx])) {
            // skip operation
            switch (op) {
                case OpCode::REG_INPUT:
                    id_idx -= 1;
                    break;
                case OpCode::REG_OUTPUT:
                    id_idx -= 2;
                    break;
                case OpCode::ADD:
                case OpCode::SUB:
                    id_idx -= 3;
                    break;
                case OpCode::MUL:
                    val_idx -= 2;
                    id_idx -= 3;
                    break;
                case OpCode::ADD_C:
                case OpCode::SUB_C:
                    id_idx -= 2;
                    break;
                case OpCode::MUL_C:
                case OpCode::NORM:
                case OpCode::ABS:
                case OpCode::EXP:
                case OpCode::LOG:
                case OpCode::ERF:
                case OpCode::ERFC:
                case OpCode::INV:
                    val_idx -= 1;
                    id_idx -= 2;
                    break;
                case OpCode::SQRT:
                case OpCode::COS:
                case OpCode::POW_C:
                    val_idx -= 2;
                    id_idx -= 2;
                    break;
            }

            if constexpr (Reset && ResetInPlace) {
                this->derivatives.resize(this->derivatives.size() - 1);
            }
            continue;
        }

        switch (op) {
            case OpCode::REG_INPUT: {
                id_idx -= 1;
                break;
            }
            case OpCode::REG_OUTPUT: {
                id_idx -= 2;

                std::size_t const arg_id = ids[id_idx];
                std::size_t const res_id = ids[id_idx + 1];

                auto const& der_list = this->derivatives_input_map[res_id];
                for (auto const& der_pair : der_list) {
                    // on the output we SET the derivative
                    this->derivatives[res_id][der_pair.first] = der_pair.second;
                    // on the dependent variable we ADD the derivative
                    this->derivatives[arg_id][der_pair.first] += der_pair.second;
                }

                if constexpr (Reset) {
                    // delete the derivative input from the map to reset
                    this->derivatives_input_map.erase(res_id);
                }
                break;
            }
            case OpCode::ADD: {
                id_idx -= 3;

                std::size_t const lhs_id = ids[id_idx];
                std::size_t const rhs_id = ids[id_idx + 1];
                std::size_t const res_id = ids[id_idx + 2];

                auto const& der_list = this->derivatives[res_id];
                for (auto const& der_pair : der_list) {
                    std::size_t const der_id = der_pair.first;
                    double const der_value = der_pair.second;

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
                        add_derivative(lhs_id, der_id, (lhs_id == der_id ? 2. : 1.) * der_value);
                        add_derivative(rhs_id, der_id, (rhs_id == der_id ? 2. : 1.) * der_value);
                    }
                }
                break;
            }
            case OpCode::SUB: {
                id_idx -= 3;

                std::size_t const lhs_id = ids[id_idx];
                std::size_t const rhs_id = ids[id_idx + 1];
                std::size_t const res_id = ids[id_idx + 2];

                auto const& der_list = this->derivatives[res_id];
                for (auto const& der_pair : der_list) {
                    std::size_t const der_id = der_pair.first;
                    double const der_value = der_pair.second;

                    if (der_id == passive_id<std::size_t>) {
                        add_derivative(lhs_id, der_id, der_value);
                        add_derivative(rhs_id, der_id, -der_value);
                    }
                    else if (der_id == res_id) {
                        add_derivative(lhs_id, rhs_id, -der_value);
                        add_derivative(lhs_id, lhs_id, der_value);
                        add_derivative(rhs_id, rhs_id, der_value);
                    }
                    else {
                        add_derivative(lhs_id, der_id, (lhs_id == der_id ? 2. : 1.) * der_value);
                        add_derivative(rhs_id, der_id, -(rhs_id == der_id ? 2. : 1.) * der_value);
                    }
                }
                break;
            }
            case OpCode::MUL: {
                val_idx -= 2;
                id_idx -= 3;

                double const lhs_val = vals[val_idx];
                double const rhs_val = vals[val_idx + 1];

                std::size_t const lhs_id = ids[id_idx];
                std::size_t const rhs_id = ids[id_idx + 1];
                std::size_t const res_id = ids[id_idx + 2];

                auto const& der_list = this->derivatives[res_id];
                for (auto const& der_pair : der_list) {
                    std::size_t const der_id = der_pair.first;
                    double const der_value = der_pair.second;

                    if (der_id == passive_id<std::size_t>) {
                        add_derivative(lhs_id, der_id, rhs_val * der_value);
                        add_derivative(rhs_id, der_id, lhs_val * der_value);
                        add_derivative(rhs_id, lhs_id, der_value);
                    }
                    else if (der_id == res_id) {
                        add_derivative(rhs_id, lhs_id, rhs_val * lhs_val * der_value);
                        add_derivative(lhs_id, lhs_id, rhs_val * rhs_val * der_value);
                        add_derivative(rhs_id, rhs_id, lhs_val * lhs_val * der_value);
                    }
                    else {
                        add_derivative(lhs_id, der_id, (lhs_id == der_id ? 2. : 1.) * rhs_val * der_value);
                        add_derivative(rhs_id, der_id, (rhs_id == der_id ? 2. : 1.) * lhs_val * der_value);
                    }
                }
                break;
            }
            case OpCode::ADD_C: {
                id_idx -= 2;

                std::size_t const arg_id = ids[id_idx];
                std::size_t const res_id = ids[id_idx + 1];

                auto const& der_list = this->derivatives[res_id];
                for (auto const& der_pair : der_list) {
                    std::size_t const der_id = der_pair.first;
                    double const der_value = der_pair.second;
                    if (der_id == passive_id<std::size_t>) {
                        add_derivative(arg_id, der_id, der_value);
                    }
                    else if (der_id == res_id) {
                        add_derivative(arg_id, arg_id, der_value);
                    }
                    else {
                        add_derivative(arg_id, der_id, (arg_id == der_id ? 2. : 1.) * der_value);
                    }
                }
                break;
            }
            case OpCode::SUB_C: {
                id_idx -= 2;

                std::size_t const arg_id = ids[id_idx];
                std::size_t const res_id = ids[id_idx + 1];

                auto const& der_list = this->derivatives[res_id];
                for (auto const& der_pair : der_list) {
                    std::size_t const der_id = der_pair.first;
                    double const der_value = der_pair.second;
                    if (der_id == passive_id<std::size_t>) {
                        add_derivative(arg_id, der_id, -der_value);
                    }
                    else if (der_id == res_id) {
                        add_derivative(arg_id, arg_id, der_value);
                    }
                    else {
                        add_derivative(arg_id, der_id, (arg_id == der_id ? 2. : 1.) * -der_value);
                    }
                }
                break;
            }
            case OpCode::MUL_C: {
                val_idx -= 1;
                id_idx -= 2;

                double const arg_val = vals[val_idx];

                std::size_t const arg_id = ids[id_idx];
                std::size_t const res_id = ids[id_idx + 1];

                auto const& der_list = this->derivatives[res_id];
                for (auto const& der_pair : der_list) {
                    std::size_t const der_id = der_pair.first;
                    double const der_value = der_pair.second;
                    if (der_id == passive_id<std::size_t>) {
                        add_derivative(arg_id, der_id, der_value * arg_val);
                    }
                    else if (der_id == res_id) {
                        add_derivative(arg_id, arg_id, arg_val * arg_val * der_value);
                    }
                    else {
                        add_derivative(arg_id, der_id, (arg_id == der_id ? 2. : 1.) * arg_val * der_value);
                    }
                }
                break;
            }
            case OpCode::NORM: {
                val_idx -= 1;
                id_idx -= 2;

                double const der_local_1 = 2.0 * vals[val_idx];
                double const der_local_2 = 2.0;

                std::size_t const arg_id = ids[id_idx];
                std::size_t const res_id = ids[id_idx + 1];

                auto const& der_list = this->derivatives[res_id];
                for (auto const& der_pair : der_list) {
                    std::size_t const der_id = der_pair.first;
                    double const der_value = der_pair.second;
                    if (der_id == passive_id<std::size_t>) {
                        add_derivative(arg_id, der_id, der_value * der_local_1);
                        add_derivative(arg_id, arg_id, der_value * der_local_2);
                    }
                    else if (der_id == res_id) {
                        add_derivative(arg_id, arg_id, der_value * der_local_1 * der_local_1);
                    }
                    else {
                        add_derivative(arg_id, der_id, (arg_id == der_id ? 2. : 1.) * der_value * der_local_1);
                    }
                }
                break;
            }
            case OpCode::EXP: {
                val_idx -= 1;
                id_idx -= 2;

                double const der_exp = vals[val_idx];

                std::size_t const arg_id = ids[id_idx];
                std::size_t const res_id = ids[id_idx + 1];

                auto const& der_list = this->derivatives[res_id];
                for (auto const& der_pair : der_list) {
                    std::size_t const der_id = der_pair.first;
                    double const der_value = der_pair.second;
                    if (der_id == passive_id<std::size_t>) {
                        double const der_exp_local = der_value * der_exp;
                        add_derivative(arg_id, der_id, der_exp_local);
                        add_derivative(arg_id, arg_id, der_exp_local);
                    }
                    else if (der_id == res_id) {
                        add_derivative(arg_id, arg_id, der_exp * der_exp * der_value);
                    }
                    else {
                        add_derivative(arg_id, der_id, (arg_id == der_id ? 2. : 1.) * der_exp * der_value);
                    }
                }
                break;
            }
            case OpCode::LOG: {
                val_idx -= 1;
                id_idx -= 2;

                double const der_local_1 = 1.0 / vals[val_idx];
                double const der_local_2 = -der_local_1 * der_local_1;

                std::size_t const arg_id = ids[id_idx];
                std::size_t const res_id = ids[id_idx + 1];

                auto const& der_list = this->derivatives[res_id];
                for (auto const& der_pair : der_list) {
                    std::size_t const der_id = der_pair.first;
                    double const der_value = der_pair.second;
                    if (der_id == passive_id<std::size_t>) {
                        add_derivative(arg_id, der_id, der_value * der_local_1);
                        add_derivative(arg_id, arg_id, der_value * der_local_2);
                    }
                    else if (der_id == res_id) {
                        add_derivative(arg_id, arg_id, der_value * der_local_1 * der_local_1);
                    }
                    else {
                        add_derivative(arg_id, der_id, (arg_id == der_id ? 2. : 1.) * der_value * der_local_1);
                    }
                }
                break;
            }
            case OpCode::ERF: {
                val_idx -= 1;
                id_idx -= 2;
                constexpr double minus_two_over_root_pi = 2. * std::numbers::inv_sqrtpi_v<double>;
                double const der_local_1 = std::exp(-vals[val_idx] * vals[val_idx]) * minus_two_over_root_pi;
                double const der_local_2 = -2. * vals[val_idx] * der_local_1;
                std::size_t const arg_id = ids[id_idx];
                std::size_t const res_id = ids[id_idx + 1];

                auto const& der_list = this->derivatives[res_id];
                for (auto const& der_pair : der_list) {
                    std::size_t const der_id = der_pair.first;
                    double const der_value = der_pair.second;
                    if (der_id == passive_id<std::size_t>) {
                        add_derivative(arg_id, der_id, der_value * der_local_1);
                        add_derivative(arg_id, arg_id, der_value * der_local_2);
                    }
                    else if (der_id == res_id) {
                        add_derivative(arg_id, arg_id, der_value * der_local_1 * der_local_1);
                    }
                    else {
                        add_derivative(arg_id, der_id, (arg_id == der_id ? 2. : 1.) * der_value * der_local_1);
                    }
                }
                break;
            }
            case OpCode::ERFC: {
                val_idx -= 1;
                id_idx -= 2;
                constexpr double minus_two_over_root_pi = -2. * std::numbers::inv_sqrtpi_v<double>;
                double const der_local_1 = std::exp(-vals[val_idx] * vals[val_idx]) * minus_two_over_root_pi;
                double const der_local_2 = -2. * vals[val_idx] * der_local_1;
                std::size_t const arg_id = ids[id_idx];
                std::size_t const res_id = ids[id_idx + 1];

                auto const& der_list = this->derivatives[res_id];
                for (auto const& der_pair : der_list) {
                    std::size_t const der_id = der_pair.first;
                    double const der_value = der_pair.second;
                    if (der_id == passive_id<std::size_t>) {
                        add_derivative(arg_id, der_id, der_value * der_local_1);
                        add_derivative(arg_id, arg_id, der_value * der_local_2);
                    }
                    else if (der_id == res_id) {
                        add_derivative(arg_id, arg_id, der_value * der_local_1 * der_local_1);
                    }
                    else {
                        add_derivative(arg_id, der_id, (arg_id == der_id ? 2. : 1.) * der_value * der_local_1);
                    }
                }
                break;
            }
            case OpCode::SQRT: {
                val_idx -= 2;
                id_idx -= 2;

                double const one_over_in = 1. / vals[val_idx];

                double const der_local_1 = 0.5 * vals[val_idx + 1] * one_over_in;
                double const der_local_2 = -0.5 * der_local_1 * one_over_in;

                std::size_t const arg_id = ids[id_idx];
                std::size_t const res_id = ids[id_idx + 1];

                auto const& der_list = this->derivatives[res_id];
                for (auto const& der_pair : der_list) {
                    std::size_t const der_id = der_pair.first;
                    double const der_value = der_pair.second;
                    if (der_id == passive_id<std::size_t>) {
                        add_derivative(arg_id, der_id, der_value * der_local_1);
                        add_derivative(arg_id, arg_id, der_value * der_local_2);
                    }
                    else if (der_id == res_id) {
                        add_derivative(arg_id, arg_id, der_value * der_local_1 * der_local_1);
                    }
                    else {
                        add_derivative(arg_id, der_id, (arg_id == der_id ? 2. : 1.) * der_value * der_local_1);
                    }
                }
                break;
            }
            case OpCode::POW_C: {
                val_idx -= 2;
                id_idx -= 2;

                double const lhs_arg = vals[val_idx];
                double const rhs_arg = vals[val_idx + 1];

                double const der_local_1 = rhs_arg != 0.0 ? rhs_arg * std::pow(lhs_arg, rhs_arg - 1.) : 0.0;
                double const der_local_2 =
                  rhs_arg != 0.0 && rhs_arg != 1.0 ? rhs_arg * (rhs_arg - 1.0) * std::pow(lhs_arg, rhs_arg - 2.) : 0.0;

                std::size_t const arg_id = ids[id_idx];
                std::size_t const res_id = ids[id_idx + 1];

                auto const& der_list = this->derivatives[res_id];
                for (auto const& der_pair : der_list) {
                    std::size_t const der_id = der_pair.first;
                    double const der_value = der_pair.second;
                    if (der_id == passive_id<std::size_t>) {
                        add_derivative(arg_id, der_id, der_value * der_local_1);
                        add_derivative(arg_id, arg_id, der_value * der_local_2);
                    }
                    else if (der_id == res_id) {
                        add_derivative(arg_id, arg_id, der_value * der_local_1 * der_local_1);
                    }
                    else {
                        add_derivative(arg_id, der_id, (arg_id == der_id ? 2. : 1.) * der_value * der_local_1);
                    }
                }
                break;
            }
            case OpCode::COS: {
                val_idx -= 2;
                id_idx -= 2;

                double const der_local_1 = -std::sin(vals[val_idx]);
                double const der_local_2 = -vals[val_idx + 1];

                std::size_t const arg_id = ids[id_idx];
                std::size_t const res_id = ids[id_idx + 1];

                auto const& der_list = this->derivatives[res_id];
                for (auto const& der_pair : der_list) {
                    std::size_t const der_id = der_pair.first;
                    double const der_value = der_pair.second;
                    if (der_id == passive_id<std::size_t>) {
                        add_derivative(arg_id, der_id, der_value * der_local_1);
                        add_derivative(arg_id, arg_id, der_value * der_local_2);
                    }
                    else if (der_id == res_id) {
                        add_derivative(arg_id, arg_id, der_value * der_local_1 * der_local_1);
                    }
                    else {
                        add_derivative(arg_id, der_id, (arg_id == der_id ? 2. : 1.) * der_value * der_local_1);
                    }
                }
                break;
            }
            case OpCode::INV: {
                val_idx -= 1;
                id_idx -= 2;

                double const der_local_1 = -vals[val_idx] * vals[val_idx];
                double const der_local_2 = -2.0 * der_local_1 * vals[val_idx];

                std::size_t const arg_id = ids[id_idx];
                std::size_t const res_id = ids[id_idx + 1];

                auto const& der_list = this->derivatives[res_id];
                for (auto const& der_pair : der_list) {
                    std::size_t const der_id = der_pair.first;
                    double const der_value = der_pair.second;
                    if (der_id == passive_id<std::size_t>) {
                        add_derivative(arg_id, der_id, der_value * der_local_1);
                        add_derivative(arg_id, arg_id, der_value * der_local_2);
                    }
                    else if (der_id == res_id) {
                        add_derivative(arg_id, arg_id, der_value * der_local_1 * der_local_1);
                    }
                    else {
                        add_derivative(arg_id, der_id, (arg_id == der_id ? 2. : 1.) * der_value * der_local_1);
                    }
                }
                break;
            }
            case OpCode::ABS: {
                val_idx -= 1;
                id_idx -= 2;

                double const der_local_1 = std::copysign(1.0, vals[val_idx]);

                std::size_t const arg_id = ids[id_idx];
                std::size_t const res_id = ids[id_idx + 1];

                auto const& der_list = this->derivatives[res_id];
                for (auto const& der_pair : der_list) {
                    std::size_t const der_id = der_pair.first;
                    double const der_value = der_pair.second;
                    if (der_id == passive_id<std::size_t>) {
                        add_derivative(arg_id, der_id, der_value * der_local_1);
                    }
                    else if (der_id == res_id) {
                        add_derivative(arg_id, arg_id, der_value * der_local_1 * der_local_1);
                    }
                    else {
                        add_derivative(arg_id, der_id, (arg_id == der_id ? 2. : 1.) * der_value * der_local_1);
                    }
                }
                break;
            }
        }

        if constexpr (Reset && ResetInPlace) {
            this->derivatives.resize(this->derivatives.size() - 1);
        }
    }

    if constexpr (Reset && !ResetInPlace) {
        this->derivatives.resize(to);
    }

    if constexpr (Reset) {
        data.ops.resize(to);
        data.vals.resize(val_idx);
        data.ids.resize(id_idx);
        data.next_id = to;

        use_op.resize(to);
    }
}

} // namespace adhoc

#endif // ADHOC_BACKPROPAGATOR2_HPP
