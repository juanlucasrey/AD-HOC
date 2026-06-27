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

#include "buffer_t.hpp"
#include "passive_id.hpp"
#include "position_impl.hpp"
#include "tape_data.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numbers>
#include <optional>
#include <vector>

namespace adhoc {

template<class Float, bool Vectorised = false, bool ConsolidateLargeUnivariate = false>
class BackPropagatorLossyCompressed {
  private:
    std::vector<std::size_t> node_location_on_buffer;
    std::vector<std::size_t> checkpoints{ 0 };
    std::vector<buffer_t<double, Vectorised> > buffers{ buffer_t<double, Vectorised>{} };

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
            this->buffers.push_back(buffer_t<double, Vectorised>{ this->get_lanes() });
        }
    }
    void set_lanes(std::size_t num_lanes) { this->buffers = { buffer_t<double, Vectorised>{ num_lanes } }; }
    auto get_lanes() const -> std::size_t { return this->buffers.front().lanes(); }
    void reserve_input(std::size_t count_registered)
    {
        this->buffers.back().reserve(this->buffers.back().size() + count_registered);
    }
    void reserve_output(std::size_t count_registered)
    {
        this->buffers.back().reserve(this->buffers.back().size() + count_registered);
    }
    void register_variable(std::size_t var_id)
    {
        this->node_location_on_buffer.resize(std::max(this->node_location_on_buffer.size(), var_id + 1),
                                             passive_id<std::size_t>);
        std::size_t& var_pos = this->node_location_on_buffer[var_id];
        if (var_pos == passive_id<std::size_t>) {
            auto it = std::upper_bound(this->checkpoints.cbegin(), this->checkpoints.cend(), var_id);
            auto& var_buffer = this->buffers[std::distance(this->checkpoints.cbegin(), it) - 1];
            var_pos = var_buffer.template get_new_loc<true>();
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
            var_pos = var_buffer.template get_new_loc<true>();
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

        auto it = std::upper_bound(this->checkpoints.cbegin(), this->checkpoints.cend(), var_id);
        auto buffed_id = static_cast<std::uint8_t>(std::distance(this->checkpoints.cbegin(), it) - 1);

        if constexpr (Vectorised) {
            auto val = this->buffers[buffed_id][var_pos];

            if (lane >= val.size()) {
                // lane value too large
                throw;
            }

            val[lane] = deriv;
        }
        else {
            if (lane != 0) {
                // This backpropagator is not designed for multiple lanes
                throw;
            }
            this->buffers[buffed_id][var_pos] = deriv;
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

        auto it = std::upper_bound(this->checkpoints.cbegin(), this->checkpoints.cend(), var_id);

        if constexpr (Vectorised) {
            auto const val = this->buffers[std::distance(this->checkpoints.cbegin(), it) - 1][var_pos];

            if (lane >= val.size()) {
                // lane value too large
                throw;
            }

            return val[lane];
        }
        else {
            if (lane != 0) {
                // This backpropagator is not designed for multiple lanes
                throw;
            }
            return this->buffers[std::distance(this->checkpoints.cbegin(), it) - 1][var_pos];
        }
    }

    auto get_derivative(std::size_t /* var_id1 */, std::size_t /* var_id2 */, std::size_t /* lane */) const -> double
    {
        return 0.;
    }

    void clear() {}

    void zero_adjoints()
    {
        for (auto& b : this->buffers) {
            b.zero();
        }
    }

    auto size_of(bool capacity = false) const -> std::size_t
    {
        std::size_t size = 0;
        size += sizeof(std::size_t); // m_num_lanes
        size += sizeof(std::size_t) * (capacity ? node_location_on_buffer.capacity() : node_location_on_buffer.size());
        size += sizeof(std::size_t) * (capacity ? this->checkpoints.capacity() : this->checkpoints.size());
        for (const auto& buffer : this->buffers) {
            size += buffer.size_of(capacity);
        }
        return size;
    }

    template<bool Reset, bool ResetInPlace, bool Log, class TapeDataType>
    void backpropagate_to(PositionImpl const& pos, TapeDataType& data);
};

template<class Float, bool Vectorised, bool ConsolidateLargeUnivariate>
template<bool Reset, bool ResetInPlace, bool Log, class TapeDataType>
void
BackPropagatorLossyCompressed<Float, Vectorised, ConsolidateLargeUnivariate>::backpropagate_to(PositionImpl const& pos,
                                                                                               TapeDataType& data)
{
    std::size_t to = pos.op_position;
    std::size_t from = data.next_id;
    if (from == checkpoints.back()) {
        checkpoints.pop_back();
        this->buffers.pop_back();
    }

    const auto& ops = data.ops;
    const auto& vals = data.vals;
    const auto& ids = data.ids;

    std::size_t val_idx = vals.size();
    std::size_t id_idx = ids.size();

    this->node_location_on_buffer.resize(ops.size(), passive_id<std::size_t>);

    std::vector<std::size_t> number_dependents(from - to);

    // LOOP 1: backward, to count number of dependents for each node and detect which nodes are active
    id_idx = ids.size();
    for (std::size_t op_idx = from; op_idx-- > to;) {
        OpCode const& op = ops[op_idx];
        bool const use_this_op =
          this->node_location_on_buffer[op_idx] != passive_id<std::size_t> || (number_dependents[op_idx - to] > 0);

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

                if (this->node_location_on_buffer[op_idx] != passive_id<std::size_t>) {
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

                if (this->node_location_on_buffer[op_idx] != passive_id<std::size_t>) {
                    // if the node is already active, it measn it had at least on precedent
                    // on other buffers
                    number_dependents[op_idx - to] += 1;
                }
                break;
            }
        }
    }
    std::size_t const id_idx_start = id_idx;
    std::size_t const val_idx_start = pos.val_position;

    std::vector<std::size_t> multiplier_origin((from - to) * 2, passive_id<std::size_t>);

    // LOOP 2: forward, to calculate multipliers after compressing induced paths
    buffer_t<double> buffer_multipliers_values;

    enum class mul_type : std::uint8_t {
        ANY,
        ONE,
        MINUS_ONE,
    };

    struct multiplier_info_t {
        std::size_t position;
        mul_type value_type;
        std::size_t loc_from;
        bool keep_alive;
    };
    buffer_t<multiplier_info_t> buffer_multipliers;

    auto copy_m = [&](std::size_t const pos, double const multiplier) { buffer_multipliers_values[pos] = multiplier; };

    auto copy_minus_m = [&](std::size_t const pos, double const multiplier) {
        buffer_multipliers_values[pos] = -multiplier;
    };

    auto minus_inplace_m = [&](std::size_t const pos) {
        buffer_multipliers_values[pos] = -buffer_multipliers_values[pos];
    };

    auto mul_inplace_m = [&](std::size_t const pos, double const multiplier) {
        buffer_multipliers_values[pos] *= multiplier;
    };

    auto add_internal = [&](std::size_t const pos1, std::size_t const pos2) {
        buffer_multipliers_values[pos1] += buffer_multipliers_values[pos2];
    };

    auto add_one = [&](std::size_t const pos1) { buffer_multipliers_values[pos1] += 1.0; };

    auto sub_one = [&](std::size_t const pos1) { buffer_multipliers_values[pos1] -= 1.0; };

    auto mul_internal = [&](std::size_t const pos1, std::size_t const pos2) {
        buffer_multipliers_values[pos1] *= buffer_multipliers_values[pos2];
    };

    auto multiplier_set_incoming = [&](std::size_t const pos, double const multiplier) {
        auto const real_pos = buffer_multipliers[pos].position;
        copy_m(real_pos, multiplier);
    };

    auto multiplier_add = [&](std::size_t const pos1, std::size_t const pos2) {
        auto& info1 = buffer_multipliers[pos1];
        auto& info2 = buffer_multipliers[pos2];

        auto const mult_type1 = info1.value_type;
        auto const mult_type2 = info2.value_type;
        if (mult_type1 == mul_type::ANY && mult_type2 == mul_type::ANY) {
            auto const pos1_buffer_loc = info1.position;
            auto const pos2_buffer_loc = info2.position;
            add_internal(pos1_buffer_loc, pos2_buffer_loc);
            buffer_multipliers_values.free_loc(pos2_buffer_loc);

            info2.position = passive_id<std::size_t>;
            // info1.position = pos1_buffer_loc;
        }
        else if (mult_type1 == mul_type::ANY) {
            auto const pos1_buffer_loc = info1.position;
            if (mult_type2 == mul_type::ONE) {
                add_one(pos1_buffer_loc);
            }
            else if (mult_type2 == mul_type::MINUS_ONE) {
                sub_one(pos1_buffer_loc);
            }

            // info1.position = pos1_buffer_loc;
        }
        else if (mult_type2 == mul_type::ANY) {
            auto const pos2_buffer_loc = info2.position;
            if (mult_type1 == mul_type::ONE) {
                add_one(pos2_buffer_loc);
            }
            else if (mult_type1 == mul_type::MINUS_ONE) {
                sub_one(pos2_buffer_loc);
            }

            info1.position = pos2_buffer_loc;
        }
        else {
            auto const new_pos = buffer_multipliers_values.get_new_loc();

            if (mult_type1 == mult_type2) {
                if (mult_type1 == mul_type::ONE) {
                    add_one(new_pos);
                    add_one(new_pos);
                }
                else {
                    sub_one(new_pos);
                    sub_one(new_pos);
                }
            }

            info1.position = new_pos;
        }

        info1.value_type = mul_type::ANY;
    };

    auto multiplier_multiply = [&](std::size_t const pos1, std::size_t const pos2) {
        auto& info1 = buffer_multipliers[pos1];
        auto& info2 = buffer_multipliers[pos2];

        auto const mult_type1 = info1.value_type;
        auto const mult_type2 = info2.value_type;
        if (mult_type2 != mul_type::ANY) {
            // should not happen
            throw;
        }

        auto const pos2_buffer_loc = info2.position;
        if (mult_type1 == mul_type::ANY) {
            auto const pos1_buffer_loc = info1.position;

            mul_internal(pos1_buffer_loc, pos2_buffer_loc);
            buffer_multipliers_values.free_loc(pos2_buffer_loc);

            info2.position = passive_id<std::size_t>;
            info1.position = pos1_buffer_loc;
        }
        else if (mult_type1 == mul_type::ONE) {
            info1.position = pos2_buffer_loc;
        }
        else {
            minus_inplace_m(pos2_buffer_loc);
            info1.position = pos2_buffer_loc;
        }

        info1.value_type = mul_type::ANY;
    };

    auto get_mult_loc = [&]<mul_type M>() -> std::size_t {
        std::size_t result = buffer_multipliers.get_new_loc();
        auto& node_data = buffer_multipliers[result];
        node_data.value_type = M;

        if constexpr (M == mul_type::ANY) {
            std::size_t size_buff = buffer_multipliers_values.get_new_loc();
            node_data.position = size_buff;
        }
        else {
            node_data.position = passive_id<std::size_t>;
        }
        return result;
    };

    auto update_univariate_aux = [&]<mul_type M, bool KeepAlive = false>(std::size_t const arg_id,
                                                                         std::size_t& mult_origin_res,
                                                                         double const multiplier = 0) -> bool {
        bool const arg_is_induced_path = (arg_id >= to) && (number_dependents[arg_id - to] == 1) &&
                                         (multiplier_origin[(arg_id - to) * 2] != passive_id<std::size_t>) &&
                                         (multiplier_origin[((arg_id - to) * 2) + 1] == passive_id<std::size_t>);

        if (arg_is_induced_path) {
            auto& mul_origin_arg = multiplier_origin[(arg_id - to) * 2];
            auto& node_data = buffer_multipliers[mul_origin_arg];
            auto& origin_type = node_data.value_type;
            auto& buffer_origin = node_data.position;

            if constexpr (M == mul_type::ANY) {
                if (origin_type != mul_type::ANY) {
                    buffer_origin = buffer_multipliers_values.get_new_loc();
                    if (origin_type == mul_type::ONE) {
                        copy_m(buffer_origin, multiplier);
                    }
                    else if (origin_type == mul_type::MINUS_ONE) {
                        copy_minus_m(buffer_origin, multiplier);
                    }
                    origin_type = mul_type::ANY;
                }
                else {
                    mul_inplace_m(buffer_origin, multiplier);
                }
            }
            else {
                if (origin_type == mul_type::ANY) {
                    if constexpr (M == mul_type::MINUS_ONE) {
                        minus_inplace_m(buffer_origin);
                    }
                }
                else if (origin_type == M) {
                    origin_type = mul_type::ONE;
                }
                else {
                    origin_type = mul_type::MINUS_ONE;
                }
            }

            node_data.keep_alive = KeepAlive;

            std::swap(mul_origin_arg, mult_origin_res);
            --number_dependents[arg_id - to];
        }
        else {
            mult_origin_res = get_mult_loc.template operator()<M>();
            if constexpr (M == mul_type::ANY) {
                multiplier_set_incoming(mult_origin_res, multiplier);
            }
            auto& node_data = buffer_multipliers[mult_origin_res];
            node_data.loc_from = arg_id;
            node_data.keep_alive = KeepAlive;
        }

        return arg_is_induced_path;
    };

    auto update_univariate = [&]<mul_type M, bool KeepAlive = false>(
                               std::size_t const arg_id, std::size_t const res_id, double const multiplier = 0) {
        auto& mult_origin_res = multiplier_origin[(res_id - to) * 2];
        update_univariate_aux.template operator()<M, KeepAlive>(arg_id, mult_origin_res, multiplier);
    };

    auto update_bivariate = [&]<mul_type M1, mul_type M2>(std::size_t const lhs_id,
                                                          std::size_t const rhs_id,
                                                          std::size_t const res_id,
                                                          double const multiplier_lhs = 0,
                                                          double const multiplier_rhs = 0) {
        std::size_t multiplier_loc_lhs = passive_id<std::size_t>;
        bool const lhs_is_induced_path =
          update_univariate_aux.template operator()<M1, false>(lhs_id, multiplier_loc_lhs, multiplier_lhs);

        std::size_t multiplier_loc_rhs = passive_id<std::size_t>;
        bool const rhs_is_induced_path =
          update_univariate_aux.template operator()<M2, false>(rhs_id, multiplier_loc_rhs, multiplier_rhs);

        bool const has_induced_path = lhs_is_induced_path || rhs_is_induced_path;
        // there is potential for a bivariate operator on the same argument, we need to check if
        // this is the case and update the multiplier if so
        auto& node_data_lhs = buffer_multipliers[multiplier_loc_lhs];
        auto& node_data_rhs = buffer_multipliers[multiplier_loc_rhs];
        bool const bivariate_consolidate_this = has_induced_path && (node_data_lhs.loc_from == node_data_rhs.loc_from);

        if (bivariate_consolidate_this) {
            std::size_t const origin_id = node_data_lhs.loc_from;
            multiplier_add(multiplier_loc_lhs, multiplier_loc_rhs);
            buffer_multipliers.free_loc(multiplier_loc_rhs);
            multiplier_origin[(res_id - to) * 2] = multiplier_loc_lhs;

            if (origin_id >= to) {
                auto& number_dependents_to_update = number_dependents[origin_id - to];
                --number_dependents_to_update;

                bool const has_single_origin =
                  (multiplier_origin[(origin_id - to) * 2] != passive_id<std::size_t>) &&
                  (multiplier_origin[((origin_id - to) * 2) + 1] == passive_id<std::size_t>);

                bool const univariate_consolidate_this = (number_dependents_to_update == 1) && has_single_origin;
                if (univariate_consolidate_this) {
                    // this node now has only one dependent, we can reintroduce a
                    // multiplication chain
                    auto& coming_from = multiplier_origin[(origin_id - to) * 2];
                    multiplier_multiply(coming_from, multiplier_loc_lhs);
                    buffer_multipliers.free_loc(multiplier_loc_lhs);
                    multiplier_origin[(res_id - to) * 2] = coming_from;
                    coming_from = passive_id<std::size_t>;
                    --number_dependents[origin_id - to];
                }
            }
        }
        else {
            multiplier_origin[(res_id - to) * 2] = multiplier_loc_lhs;
            multiplier_origin[((res_id - to) * 2) + 1] = multiplier_loc_rhs;

            if constexpr (ConsolidateLargeUnivariate) {
                auto check_univariate_subtree = [&](std::size_t res_id) -> std::optional<std::size_t> {
                    std::optional<std::size_t> result = std::nullopt;

                    std::map<std::size_t, std::size_t> local_number_dependents;
                    local_number_dependents[res_id] = number_dependents[res_id - to];

                    do {
                        auto const top_id = local_number_dependents.rbegin()->first;
                        auto const top_deps = local_number_dependents.rbegin()->second;

                        if (top_id < to) {
                            // we reached a node that crosses boundaries
                            return result;
                        }

                        if (top_deps != number_dependents[top_id - to]) {
                            // a dependency is not included in this subtree, so we cannot consolidate
                            return result;
                        }

                        local_number_dependents.erase(top_id);
                        std::size_t const origin_multiplier_id_1 = multiplier_origin[(top_id - to) * 2];
                        std::size_t const origin_multiplier_id_2 = multiplier_origin[((top_id - to) * 2) + 1];
                        auto const& node_data_1 = buffer_multipliers.data()[origin_multiplier_id_1];
                        auto const& node_data_2 = buffer_multipliers.data()[origin_multiplier_id_2];

                        if (origin_multiplier_id_2 != passive_id<std::size_t>) {
                            std::size_t const id_origin1 = node_data_1.multiplier_loc_from;
                            std::size_t const id_origin2 = node_data_2.multiplier_loc_from;
                            ++local_number_dependents[id_origin1];
                            ++local_number_dependents[id_origin2];
                        }
                        else if (origin_multiplier_id_1 != passive_id<std::size_t>) {
                            std::size_t const id_origin = node_data_1.multiplier_loc_from;
                            ++local_number_dependents[id_origin];
                        }
                        else {
                            // we reached an input node, while number of dependents is larger than 1, so
                            // this cannot be a univariate operator
                            return result;
                        }

                    } while (local_number_dependents.size() > 1);

                    result = local_number_dependents.rbegin()->first;

                    return result;
                };

                auto compress_univariate_subtree = [&](std::size_t in_id, std::size_t res_id) {
                    // when compressing a univarate subtree, we do a forward pass.
                    // why? because a forward pass is simpler, especially for higher orders.

                    std::map<std::size_t, double> local_derivatives;
                    local_derivatives[res_id] = 1.0;

                    do {
                        auto const top_id = local_derivatives.rbegin()->first;
                        auto const top_der = local_derivatives.rbegin()->second;
                        local_derivatives.erase(top_id);

                        std::size_t const origin_multiplier_id_1 = multiplier_origin[(top_id - to) * 2];
                        std::size_t const origin_multiplier_id_2 = multiplier_origin[((top_id - to) * 2) + 1];
                        auto const& node_data_1 = buffer_multipliers.data()[origin_multiplier_id_1];
                        auto const& node_data_2 = buffer_multipliers.data()[origin_multiplier_id_2];

                        if (origin_multiplier_id_2 != passive_id<std::size_t>) {
                            std::size_t const id_origin1 = node_data_1.multiplier_loc_from;
                            std::size_t const id_origin2 = node_data_2.multiplier_loc_from;

                            std::size_t& loc_multipler_1 = multiplier_origin[(top_id - to) * 2];
                            auto& node_data_loc_1 = buffer_multipliers.data()[loc_multipler_1];
                            auto const mult_type1 = node_data_loc_1.values_type;
                            if (mult_type1 == mul_type::ONE) {
                                local_derivatives[id_origin1] += top_der;
                            }
                            else if (mult_type1 == mul_type::MINUS_ONE) {
                                local_derivatives[id_origin1] -= top_der;
                            }
                            else {
                                auto const pos1_buffer_loc = node_data_loc_1.position;
                                local_derivatives[id_origin1] += top_der * buffer_multipliers_values[pos1_buffer_loc];
                            }

                            std::size_t& loc_multipler_2 = multiplier_origin[((top_id - to) * 2) + 1];
                            auto& node_data_loc_2 = buffer_multipliers.data()[loc_multipler_2];
                            auto const mult_type2 = node_data_loc_2.values_type;
                            if (mult_type2 == mul_type::ONE) {
                                local_derivatives[id_origin2] += top_der;
                            }
                            else if (mult_type2 == mul_type::MINUS_ONE) {
                                local_derivatives[id_origin2] -= top_der;
                            }
                            else {
                                auto const pos2_buffer_loc = node_data_loc_2.position;
                                local_derivatives[id_origin2] += top_der * buffer_multipliers_values[pos2_buffer_loc];
                            }

                            // unecessary?
                            node_data_loc_1.multiplier_loc_from = passive_id<std::size_t>;
                            buffer_multipliers.free_loc(loc_multipler_1);
                            loc_multipler_1 = passive_id<std::size_t>;

                            // unecessary?
                            node_data_loc_2.multiplier_loc_from = passive_id<std::size_t>;
                            buffer_multipliers.free_loc(loc_multipler_2);
                            loc_multipler_2 = passive_id<std::size_t>;

                            if (id_origin1 >= to) {
                                --number_dependents[id_origin1 - to];
                            }
                            if (id_origin2 >= to) {
                                --number_dependents[id_origin2 - to];
                            }
                        }
                        else if (origin_multiplier_id_1 != passive_id<std::size_t>) {
                            std::size_t const id_origin =
                              buffer_multipliers.data()[origin_multiplier_id_1].multiplier_loc_from;

                            std::size_t& loc_multipler = multiplier_origin[(top_id - to) * 2];
                            auto& node_data_loc = buffer_multipliers.data()[loc_multipler];
                            auto const mult_type = node_data_loc.values_type;
                            if (mult_type == mul_type::ONE) {
                                local_derivatives[id_origin] += top_der;
                            }
                            else if (mult_type == mul_type::MINUS_ONE) {
                                local_derivatives[id_origin] -= top_der;
                            }
                            else {
                                auto const pos_buffer_loc = node_data_loc.position;
                                local_derivatives[id_origin] += top_der * buffer_multipliers_values[pos_buffer_loc];
                            }

                            node_data_loc.multiplier_loc_from = passive_id<std::size_t>;
                            buffer_multipliers.free_loc(loc_multipler);
                            loc_multipler = passive_id<std::size_t>;

                            if (id_origin >= to) {
                                --number_dependents[id_origin - to];
                            }
                        }
                        else {
                            // we reached an input node, while number of dependents is larger than 1, so
                            // this cannot be a univariate operator
                            throw;
                        }

                    } while (local_derivatives.size() > 1);

                    if (local_derivatives.begin()->first != in_id) {
                        // this should not happen, because we should have detected the univariate subtree correctly
                        throw;
                    }

                    std::size_t const new_pos = get_mult_loc.template operator()<mul_type::ANY>();
                    auto& node_data_loc = buffer_multipliers.data()[new_pos];

                    node_data_loc.multiplier_loc_from = in_id;
                    multiplier_set_incoming(new_pos, local_derivatives.begin()->second);

                    multiplier_origin[(res_id - to) * 2] = new_pos;
                    if (in_id >= to) {
                        ++number_dependents[in_id - to];

                        bool const has_single_origin =
                          (multiplier_origin[(in_id - to) * 2] != passive_id<std::size_t>) &&
                          (multiplier_origin[((in_id - to) * 2) + 1] == passive_id<std::size_t>);

                        bool const univariate_consolidate_this =
                          (number_dependents[in_id - to] == 1) && has_single_origin;

                        if (univariate_consolidate_this) {
                            // this node now has only one dependent, we can reintroduce a
                            // multiplication chain
                            auto& coming_from = multiplier_origin[(in_id - to) * 2];
                            multiplier_multiply(coming_from, new_pos);

                            node_data_loc.multiplier_loc_from = passive_id<std::size_t>;
                            buffer_multipliers.free_loc(new_pos);

                            multiplier_origin[(res_id - to) * 2] = coming_from;
                            coming_from = passive_id<std::size_t>;
                            --number_dependents[in_id - to];
                        }
                    }
                };

                auto in_id_opt = check_univariate_subtree(res_id);
                if (in_id_opt) {
                    compress_univariate_subtree(*in_id_opt, res_id);
                }
            }
        }
    };

    id_idx = id_idx_start;
    val_idx = val_idx_start;

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
                    update_bivariate.template operator()<mul_type::ONE, mul_type::MINUS_ONE>(lhs_id, rhs_id, res_id);
                }
                id_idx += 3;
                break;
            }
            case OpCode::MUL: {
                if (use_this_op) {
                    std::size_t const lhs_id = ids[id_idx];
                    std::size_t const rhs_id = ids[id_idx + 1];
                    std::size_t const res_id = ids[id_idx + 2];
                    double const lhs_val = vals[val_idx];
                    double const rhs_val = vals[val_idx + 1];
                    update_bivariate.template operator()<mul_type::ANY, mul_type::ANY>(
                      lhs_id, rhs_id, res_id, rhs_val, lhs_val);
                }
                val_idx += 2;
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
                    double const multiplier = vals[val_idx];
                    update_univariate.template operator()<mul_type::ANY>(arg_id, res_id, multiplier);
                }
                val_idx += 1;
                id_idx += 2;
                break;
            }
            case OpCode::NORM: {
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    double const multiplier = 2.0 * vals[val_idx];
                    update_univariate.template operator()<mul_type::ANY>(arg_id, res_id, multiplier);
                }
                val_idx += 1;
                id_idx += 2;
                break;
            }
            case OpCode::INV: {
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    double const multiplier = -vals[val_idx] * vals[val_idx];
                    update_univariate.template operator()<mul_type::ANY>(arg_id, res_id, multiplier);
                }
                val_idx += 1;
                id_idx += 2;
                break;
            }
            case OpCode::ABS: {
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    double const multiplier = std::copysign(1.0, vals[val_idx]);
                    update_univariate.template operator()<mul_type::ANY>(arg_id, res_id, multiplier);
                }
                val_idx += 1;
                id_idx += 2;
                break;
            }
            case OpCode::EXP: {
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    double const multiplier = vals[val_idx];
                    update_univariate.template operator()<mul_type::ANY>(arg_id, res_id, multiplier);
                }
                val_idx += 1;
                id_idx += 2;
                break;
            }
            case OpCode::LOG: {
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    double const multiplier = 1.0 / vals[val_idx];
                    update_univariate.template operator()<mul_type::ANY>(arg_id, res_id, multiplier);
                }
                val_idx += 1;
                id_idx += 2;
                break;
            }
            case OpCode::ERF: {
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    constexpr double two_over_root_pi = 2. * std::numbers::inv_sqrtpi_v<double>;
                    double const multiplier = std::exp(-vals[val_idx] * vals[val_idx]) * two_over_root_pi;
                    update_univariate.template operator()<mul_type::ANY>(arg_id, res_id, multiplier);
                }
                val_idx += 1;
                id_idx += 2;
                break;
            }
            case OpCode::ERFC: {
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    constexpr double minus_two_over_root_pi = -2. * std::numbers::inv_sqrtpi_v<double>;
                    double const multiplier = std::exp(-vals[val_idx] * vals[val_idx]) * minus_two_over_root_pi;
                    update_univariate.template operator()<mul_type::ANY>(arg_id, res_id, multiplier);
                }
                val_idx += 1;
                id_idx += 2;
                break;
            }
            case OpCode::COS: {
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    double const multiplier = -std::sin(vals[val_idx]);
                    update_univariate.template operator()<mul_type::ANY>(arg_id, res_id, multiplier);
                }
                val_idx += 2;
                id_idx += 2;
                break;
            }
            case OpCode::SQRT: {
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    double const one_over_in = 1. / vals[val_idx];
                    double const multiplier = 0.5 * vals[val_idx + 1] * one_over_in;
                    update_univariate.template operator()<mul_type::ANY>(arg_id, res_id, multiplier);
                }
                val_idx += 2;
                id_idx += 2;
                break;
            }
            case OpCode::POW_C: {
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    double const lhs_arg = vals[val_idx];
                    double const rhs_arg = vals[val_idx + 1];
                    double const multiplier = rhs_arg != 0.0 ? rhs_arg * std::pow(lhs_arg, rhs_arg - 1.) : 0.0;
                    update_univariate.template operator()<mul_type::ANY>(arg_id, res_id, multiplier);
                }
                val_idx += 2;
                id_idx += 2;
                break;
            }
        }
    }

    // LOOP 3: backward, to calculate derivatives on multiple lanes
    auto& buffer_vals = this->buffers.back();

    auto copy = [&](std::uint8_t const which, std::size_t const out_pos, std::size_t const in_pos) {
        if constexpr (Vectorised) {
            auto const src = buffer_vals[out_pos];
            auto dest = this->buffers[which][in_pos];
#pragma omp simd
            for (std::size_t i = 0; i < dest.size(); ++i) {
                dest[i] = src[i];
            }
        }
        else {
            this->buffers[which][in_pos] = buffer_vals[out_pos];
        }
    };

    auto copy_minus = [&](std::uint8_t const which, std::size_t const out_pos, std::size_t const in_pos) {
        if constexpr (Vectorised) {
            auto const src = buffer_vals[out_pos];
            auto dest = this->buffers[which][in_pos];
#pragma omp simd
            for (std::size_t i = 0; i < dest.size(); ++i) {
                dest[i] = -src[i];
            }
        }
        else {
            this->buffers[which][in_pos] = -buffer_vals[out_pos];
        }
    };

    auto add = [&](std::uint8_t const which, std::size_t const out_pos, std::size_t const in_pos) {
        if constexpr (Vectorised) {
            auto const src = buffer_vals[out_pos];
            auto dest = this->buffers[which][in_pos];
#pragma omp simd
            for (std::size_t i = 0; i < dest.size(); ++i) {
                dest[i] += src[i];
            }
        }
        else {
            this->buffers[which][in_pos] += buffer_vals[out_pos];
        }
    };

    auto sub = [&](std::uint8_t const which, std::size_t const out_pos, std::size_t const in_pos) {
        if constexpr (Vectorised) {
            auto const src = buffer_vals[out_pos];
            auto dest = this->buffers[which][in_pos];
#pragma omp simd
            for (std::size_t i = 0; i < dest.size(); ++i) {
                dest[i] -= src[i];
            }
        }
        else {
            this->buffers[which][in_pos] -= buffer_vals[out_pos];
        }
    };

    auto minus_inplace = [&](std::size_t const pos) {
        if constexpr (Vectorised) {
            auto dest = buffer_vals[pos];
#pragma omp simd
            for (std::size_t i = 0; i < dest.size(); ++i) {
                dest[i] = -dest[i];
            }
        }
        else {
            buffer_vals[pos] = -buffer_vals[pos];
        }
    };

    auto mul_inplace = [&](std::size_t const pos, std::size_t const multiplier_id) {
        double const multiplier = buffer_multipliers_values[multiplier_id];
        if constexpr (Vectorised) {
            auto dest = buffer_vals[pos];
#pragma omp simd
            for (std::size_t i = 0; i < dest.size(); ++i) {
                dest[i] *= multiplier;
            }
        }
        else {
            buffer_vals[pos] *= multiplier;
        }
    };

    auto mul_add = [&](std::uint8_t const which,
                       std::size_t const out_pos,
                       std::size_t const in_pos,
                       std::size_t const multiplier_id) {
        double const multiplier = buffer_multipliers_values[multiplier_id];
        if constexpr (Vectorised) {
            auto const src = buffer_vals[out_pos];
            auto dest = this->buffers[which][in_pos];
#pragma omp simd
            for (std::size_t i = 0; i < dest.size(); ++i) {
                dest[i] += src[i] * multiplier;
            }
        }
        else {
            this->buffers[which][in_pos] += buffer_vals[out_pos] * multiplier;
        }
    };

    auto mul_set = [&](std::uint8_t const which,
                       std::size_t const out_pos,
                       std::size_t const in_pos,
                       std::size_t const multiplier_id) {
        double const multiplier = buffer_multipliers_values[multiplier_id];
        if constexpr (Vectorised) {
            auto const src = buffer_vals[out_pos];
            auto dest = this->buffers[which][in_pos];
#pragma omp simd
            for (std::size_t i = 0; i < dest.size(); ++i) {
                dest[i] = src[i] * multiplier;
            }
        }
        else {
            this->buffers[which][in_pos] = buffer_vals[out_pos] * multiplier;
        }
    };

    auto copy_mul =
      [&](std::size_t res_pos, std::size_t& arg_pos, std::uint8_t buffer_id, std::size_t const multiplier_id) {
          auto const& node_data = buffer_multipliers[multiplier_id];
          auto const multiplier_type = node_data.value_type;
          bool const arg_is_new = (arg_pos == passive_id<std::size_t>);

          if (arg_is_new) {
              auto& arg_buffer = this->buffers[buffer_id];
              arg_pos = arg_buffer.get_new_loc();
              // this is a new value, we NEED to override
              if (multiplier_type == mul_type::ANY) {
                  auto const position_buffer = node_data.position;
                  mul_set(buffer_id, res_pos, arg_pos, position_buffer);
              }
              else if (multiplier_type == mul_type::ONE) {
                  copy(buffer_id, res_pos, arg_pos);
              }
              else if (multiplier_type == mul_type::MINUS_ONE) {
                  copy_minus(buffer_id, res_pos, arg_pos);
              }
          }
          else {
              if (multiplier_type == mul_type::ANY) {
                  auto const position_buffer = node_data.position;
                  mul_add(buffer_id, res_pos, arg_pos, position_buffer);
              }
              else if (multiplier_type == mul_type::ONE) {
                  add(buffer_id, res_pos, arg_pos);
              }
              else if (multiplier_type == mul_type::MINUS_ONE) {
                  sub(buffer_id, res_pos, arg_pos);
              }
          }
      };

    auto const& checkpoints_c = this->checkpoints;
    auto get_loc = [checkpoints_c](std::size_t id) -> std::tuple<bool, std::uint8_t> {
        auto it = std::upper_bound(checkpoints_c.begin(), checkpoints_c.end(), id);
        auto const buffer_id = static_cast<std::uint8_t>(std::distance(checkpoints_c.cbegin(), it) - 1);
        return { it == checkpoints_c.end(), buffer_id };
    };

    for (std::size_t op_idx = from; op_idx-- > to;) {
        std::size_t const first_multiplier_origin = multiplier_origin[(op_idx - to) * 2];
        std::size_t const second_multiplier_origin = multiplier_origin[((op_idx - to) * 2) + 1];
        if (first_multiplier_origin != passive_id<std::size_t> && second_multiplier_origin == passive_id<std::size_t>) {
            auto const& node_data_first = buffer_multipliers[first_multiplier_origin];
            std::size_t const arg_id = node_data_first.loc_from;
            std::size_t const res_id = op_idx;
            bool const keep_alive = node_data_first.keep_alive;

            std::size_t& res_pos = this->node_location_on_buffer[res_id];
            std::size_t& arg_pos = this->node_location_on_buffer[arg_id];

            auto const arg_pos_data = get_loc(arg_id);
            bool const arg_is_new = (arg_pos == passive_id<std::size_t>);
            bool const arg_inplace = arg_is_new && std::get<0>(arg_pos_data);

            if (arg_inplace && !keep_alive) {
                auto const multiplier_type = node_data_first.value_type;
                // res id should now be arg id, avoiding a copy and a potential buffer increase
                if (multiplier_type == mul_type::ANY) {
                    auto const position_buffer = node_data_first.position;
                    mul_inplace(res_pos, position_buffer);
                }
                else if (multiplier_type == mul_type::MINUS_ONE) {
                    minus_inplace(res_pos);
                }
                std::swap(res_pos, arg_pos);
            }
            else {
                copy_mul(res_pos, arg_pos, std::get<1>(arg_pos_data), first_multiplier_origin);

                if (!keep_alive) {
                    this->buffers.back().free_loc(res_pos);
                    res_pos = passive_id<std::size_t>;
                }
            }
        }
        else if (first_multiplier_origin != passive_id<std::size_t>) {
            auto const& node_data_first = buffer_multipliers[first_multiplier_origin];
            auto const& node_data_second = buffer_multipliers[second_multiplier_origin];
            std::size_t const lhs_id = node_data_first.loc_from;
            std::size_t const rhs_id = node_data_second.loc_from;
            std::size_t const res_id = op_idx;

            std::size_t& res_pos = this->node_location_on_buffer[res_id];
            std::size_t& lhs_pos = this->node_location_on_buffer[lhs_id];
            std::size_t& rhs_pos = this->node_location_on_buffer[rhs_id];

            auto const lhs_pos_data = get_loc(lhs_id);
            auto const rhs_pos_data = get_loc(rhs_id);

            bool const lhs_is_new = (lhs_pos == passive_id<std::size_t>);
            bool const rhs_is_new = (rhs_pos == passive_id<std::size_t>);
            bool const lhs_inplace = lhs_is_new && std::get<0>(lhs_pos_data);
            bool const rhs_inplace = !lhs_inplace && rhs_is_new && std::get<0>(rhs_pos_data);

            if (!lhs_inplace) {
                copy_mul(res_pos, lhs_pos, std::get<1>(lhs_pos_data), first_multiplier_origin);
            }

            if (!rhs_inplace) {
                copy_mul(res_pos, rhs_pos, std::get<1>(rhs_pos_data), second_multiplier_origin);
            }

            if (lhs_inplace) {
                auto const mul_1_type = node_data_first.value_type;
                // res id should now be lhs id, avoiding a copy and a potential buffer increase
                if (mul_1_type == mul_type::ANY) {
                    auto const position_buffer = node_data_first.position;
                    mul_inplace(res_pos, position_buffer);
                }
                else if (mul_1_type == mul_type::MINUS_ONE) {
                    minus_inplace(res_pos);
                }

                std::swap(res_pos, lhs_pos);
            }
            else if (rhs_inplace) {
                auto const mul_2_type = node_data_second.value_type;
                // res id should now be arg id, avoiding a copy and a potential buffer increase
                if (mul_2_type == mul_type::ANY) {
                    auto const position_buffer = node_data_second.position;
                    mul_inplace(res_pos, position_buffer);
                }
                else if (mul_2_type == mul_type::MINUS_ONE) {
                    minus_inplace(res_pos);
                }

                std::swap(res_pos, rhs_pos);
            }
            else {
                // don't forget to free res_id from the buffer!
                this->buffers.back().free_loc(res_pos);
                res_pos = passive_id<std::size_t>;
            }
        }
    }

    if constexpr (Reset) {
        if (!checkpoints.empty() && to == checkpoints.back()) {
            this->buffers.back() = buffer_t<double, Vectorised>{ this->get_lanes() };
        }

        reset(pos, data);
        this->node_location_on_buffer.resize(pos.op_position);
    }
}

} // namespace adhoc

#endif // ADHOC_BACKPROPAGATOR1LOSSYCOMPRESSED_HPP
