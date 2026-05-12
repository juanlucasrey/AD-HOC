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
#include <numbers>
#include <vector>

#include <map>

namespace adhoc {

template<class Float, bool Vectorised = false, bool ConsolidateLargeUnivariate = false>
class BackPropagatorLossyCompressed {
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

template<class Float, bool Vectorised, bool ConsolidateLargeUnivariate>
template<bool Reset, bool ResetInPlace, bool Log>
void
BackPropagatorLossyCompressed<Float, Vectorised, ConsolidateLargeUnivariate>::backpropagate_to(std::size_t to,
                                                                                               TapeData& data)
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

    auto& buffer_free_positions = buffers.back().free_positions;

    std::vector<std::size_t> number_dependents(from - to);

    // LOOP 1: backward, to count number of dependents for each node and detect which nodes are active
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
            case OpCode::SQRT:
            case OpCode::POW_C: {
                val_idx -= 2;
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
    std::size_t const val_idx_start = val_idx;

    buffer_t buffer_multipliers;

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

    auto multiplier_set_incoming = [&]<mul_type M>(std::size_t const pos, double const multiplier = 0) {
        if constexpr (M == mul_type::ONE) {
            buffer_multipliers.values[pos] = 1.0;
            values_type[pos] = mul_type::ONE;
        }
        else if constexpr (M == mul_type::MINUS_ONE) {
            buffer_multipliers.values[pos] = -1.0;
            values_type[pos] = mul_type::MINUS_ONE;
        }
        else {
            static_assert(M == mul_type::ANY, "Invalid multiplier type");
            buffer_multipliers.values[pos] = multiplier;
            values_type[pos] = mul_type::ANY;
        }
    };

    auto multiplier_multiply_incoming = [&, combine_mul_type]<mul_type M>(std::size_t const pos,
                                                                          double const multiplier = 0) {
        if constexpr (M == mul_type::ONE) {
            values_type[pos] = combine_mul_type(values_type[pos], mul_type::ONE);
        }
        else if constexpr (M == mul_type::MINUS_ONE) {
            buffer_multipliers.values[pos] = -buffer_multipliers.values[pos];
            values_type[pos] = combine_mul_type(values_type[pos], mul_type::MINUS_ONE);
        }
        else {
            static_assert(M == mul_type::ANY, "Invalid multiplier type");
            buffer_multipliers.values[pos] *= multiplier;
            values_type[pos] = mul_type::ANY;
        }
    };

    auto multiplier_add = [&](std::size_t const pos1, std::size_t const pos2) {
        buffer_multipliers.values[pos1] += buffer_multipliers.values[pos2];
        values_type[pos1] = mul_type::ANY;
    };

    auto multiplier_multiply = [&](std::size_t const pos1, std::size_t const pos2) {
        buffer_multipliers.values[pos1] *= buffer_multipliers.values[pos2];
        values_type[pos1] = mul_type::ANY;
    };

    id_idx = id_idx_start;
    val_idx = val_idx_start;

    for (std::size_t op_idx = to; op_idx < from; ++op_idx) {
        OpCode const& op = ops[op_idx];
        bool const use_this_op = (number_dependents[op_idx - to] > 0);

        auto get_mult_loc = [&]() -> std::size_t {
            if (buffer_multipliers.free_positions.empty()) {
                std::size_t result = buffer_multipliers.size;
                ++buffer_multipliers.size;
                buffer_multipliers.values.resize(buffer_multipliers.values.size() + 1);
                values_type.resize(values_type.size() + 1);
                multiplier_loc_from.resize(buffer_multipliers.size);
                multiplier_keep_alive.resize(buffer_multipliers.size);
                return result;
            }

            std::size_t result = buffer_multipliers.free_positions.back();
            buffer_multipliers.free_positions.pop_back();
            return result;
        };

        auto update_univariate = [&]<mul_type M, bool KeepAlive = false>(
                                   std::size_t const arg_id, std::size_t const res_id, double const multiplier = 0) {
            bool const arg_is_induced_path = (arg_id >= to) && (number_dependents[arg_id - to] == 1) &&
                                             (multiplier_origin[(arg_id - to) * 2] != passive_id<std::size_t>) &&
                                             (multiplier_origin[((arg_id - to) * 2) + 1] == passive_id<std::size_t>);

            auto& mult_origin_res = multiplier_origin[(res_id - to) * 2];

            if (arg_is_induced_path) {
                auto& mul_origin_arg = multiplier_origin[(arg_id - to) * 2];
                if constexpr (M == mul_type::ANY) {
                    multiplier_multiply_incoming.template operator()<mul_type::ANY>(mul_origin_arg, multiplier);
                }
                else {
                    multiplier_multiply_incoming.template operator()<M>(mul_origin_arg);
                }
                multiplier_keep_alive[mul_origin_arg] = KeepAlive;

                std::swap(mul_origin_arg, mult_origin_res);
                --number_dependents[arg_id - to];
            }
            else {
                mult_origin_res = get_mult_loc();
                if constexpr (M == mul_type::ANY) {
                    multiplier_set_incoming.template operator()<mul_type::ANY>(mult_origin_res, multiplier);
                }
                else {
                    multiplier_set_incoming.template operator()<M>(mult_origin_res);
                }
                multiplier_loc_from[mult_origin_res] = arg_id;
                multiplier_keep_alive[mult_origin_res] = KeepAlive;
            }
        };

        auto update_bivariate = [&]<mul_type M1, mul_type M2>(std::size_t const lhs_id,
                                                              std::size_t const rhs_id,
                                                              std::size_t const res_id,
                                                              double const multiplier_lhs = 0,
                                                              double const multiplier_rhs = 0) {
            bool const lhs_is_induced_path = (lhs_id >= to) && (number_dependents[lhs_id - to] == 1) &&
                                             (multiplier_origin[(lhs_id - to) * 2] != passive_id<std::size_t>) &&
                                             (multiplier_origin[((lhs_id - to) * 2) + 1] == passive_id<std::size_t>);
            bool const rhs_is_induced_path = (rhs_id >= to) && (number_dependents[rhs_id - to] == 1) &&
                                             (multiplier_origin[(rhs_id - to) * 2] != passive_id<std::size_t>) &&
                                             (multiplier_origin[((rhs_id - to) * 2) + 1] == passive_id<std::size_t>);

            std::size_t multiplier_loc_lhs = passive_id<std::size_t>;
            if (lhs_is_induced_path) {
                auto& mul_origin_arg = multiplier_origin[(lhs_id - to) * 2];
                if constexpr (M1 == mul_type::ANY) {
                    multiplier_multiply_incoming.template operator()<mul_type::ANY>(mul_origin_arg, multiplier_lhs);
                }
                else {
                    multiplier_multiply_incoming.template operator()<M1>(mul_origin_arg);
                }
                std::swap(multiplier_loc_lhs, mul_origin_arg);
                --number_dependents[lhs_id - to];
            }
            else {
                multiplier_loc_lhs = get_mult_loc();
                if constexpr (M1 == mul_type::ANY) {
                    multiplier_set_incoming.template operator()<mul_type::ANY>(multiplier_loc_lhs, multiplier_lhs);
                }
                else {
                    multiplier_set_incoming.template operator()<M1>(multiplier_loc_lhs);
                }
                multiplier_loc_from[multiplier_loc_lhs] = lhs_id;
                multiplier_keep_alive[multiplier_loc_lhs] = false;
            }

            std::size_t multiplier_loc_rhs = passive_id<std::size_t>;
            if (rhs_is_induced_path) {
                auto& mul_origin_arg = multiplier_origin[(rhs_id - to) * 2];
                if constexpr (M2 == mul_type::ANY) {
                    multiplier_multiply_incoming.template operator()<mul_type::ANY>(mul_origin_arg, multiplier_rhs);
                }
                else {
                    multiplier_multiply_incoming.template operator()<M2>(mul_origin_arg);
                }
                std::swap(multiplier_loc_rhs, mul_origin_arg);
                --number_dependents[rhs_id - to];
            }
            else {
                multiplier_loc_rhs = get_mult_loc();
                if constexpr (M2 == mul_type::ANY) {
                    multiplier_set_incoming.template operator()<mul_type::ANY>(multiplier_loc_rhs, multiplier_rhs);
                }
                else {
                    multiplier_set_incoming.template operator()<M2>(multiplier_loc_rhs);
                }
                multiplier_loc_from[multiplier_loc_rhs] = rhs_id;
                multiplier_keep_alive[multiplier_loc_rhs] = false;
            }

            bool has_induced_path = lhs_is_induced_path || rhs_is_induced_path;
            // there is potential for a bivariate operator on the same argument, we need to check if
            // this is the case and update the multiplier if so
            bool bivariate_consolidate_this =
              has_induced_path && (multiplier_loc_from[multiplier_loc_lhs] == multiplier_loc_from[multiplier_loc_rhs]);

            if (bivariate_consolidate_this) {
                std::size_t origin_id = multiplier_loc_from[multiplier_loc_lhs];
                multiplier_add(multiplier_loc_lhs, multiplier_loc_rhs);

                buffer_multipliers.free_positions.push_back(multiplier_loc_rhs);
                multiplier_loc_from[multiplier_loc_rhs] = passive_id<std::size_t>;

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

                        buffer_multipliers.free_positions.push_back(multiplier_loc_lhs);
                        multiplier_loc_from[multiplier_loc_lhs] = passive_id<std::size_t>;

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

                        // this while should not be needed, because allprevious nodes should have been checked for
                        // univariate subtrees
                        //  while (true) {
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

                            if (origin_multiplier_id_2 != passive_id<std::size_t>) {
                                std::size_t const id_origin1 = multiplier_loc_from[origin_multiplier_id_1];
                                std::size_t const id_origin2 = multiplier_loc_from[origin_multiplier_id_2];
                                ++local_number_dependents[id_origin1];
                                ++local_number_dependents[id_origin2];
                            }
                            else if (origin_multiplier_id_1 != passive_id<std::size_t>) {
                                std::size_t const id_origin = multiplier_loc_from[origin_multiplier_id_1];
                                ++local_number_dependents[id_origin];
                            }
                            else {
                                // we reached an input node, while number of dependents is larger than 1, so
                                // this cannot be a univariate operator
                                return result;
                            }

                        } while (local_number_dependents.size() > 1);

                        result = local_number_dependents.rbegin()->first;
                        // this while should not be needed, because allprevious nodes should have been checked for
                        // univariate subtrees
                        // }

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

                            if (origin_multiplier_id_2 != passive_id<std::size_t>) {
                                std::size_t const id_origin1 = multiplier_loc_from[origin_multiplier_id_1];
                                std::size_t const id_origin2 = multiplier_loc_from[origin_multiplier_id_2];

                                std::size_t& loc_multipler_1 = multiplier_origin[(top_id - to) * 2];
                                std::size_t& loc_multipler_2 = multiplier_origin[((top_id - to) * 2) + 1];
                                local_derivatives[id_origin1] += top_der * buffer_multipliers.values[loc_multipler_1];
                                local_derivatives[id_origin2] += top_der * buffer_multipliers.values[loc_multipler_2];

                                buffer_multipliers.free_positions.push_back(loc_multipler_1);
                                multiplier_loc_from[loc_multipler_1] = passive_id<std::size_t>;
                                loc_multipler_1 = passive_id<std::size_t>;

                                buffer_multipliers.free_positions.push_back(loc_multipler_2);
                                multiplier_loc_from[loc_multipler_2] = passive_id<std::size_t>;
                                loc_multipler_2 = passive_id<std::size_t>;

                                if (id_origin1 >= to) {
                                    --number_dependents[id_origin1 - to];
                                }
                                if (id_origin2 >= to) {
                                    --number_dependents[id_origin2 - to];
                                }
                            }
                            else if (origin_multiplier_id_1 != passive_id<std::size_t>) {
                                std::size_t const id_origin = multiplier_loc_from[origin_multiplier_id_1];

                                std::size_t& loc_multipler = multiplier_origin[(top_id - to) * 2];
                                local_derivatives[id_origin] += top_der * buffer_multipliers.values[loc_multipler];

                                buffer_multipliers.free_positions.push_back(loc_multipler);
                                multiplier_loc_from[loc_multipler] = passive_id<std::size_t>;
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

                        std::size_t new_pos = get_mult_loc();

                        multiplier_loc_from[new_pos] = in_id;
                        buffer_multipliers.values[new_pos] = local_derivatives.begin()->second;
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

                                buffer_multipliers.free_positions.push_back(new_pos);
                                multiplier_loc_from[new_pos] = passive_id<std::size_t>;

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
    auto& buffer_vals = buffers.back().values;

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

    auto copy_mul = [this, mul_set, copy, copy_minus, mul_add, add, sub](std::size_t res_pos,
                                                                         std::size_t& arg_pos,
                                                                         std::uint8_t buffer_id,
                                                                         double multiplier,
                                                                         mul_type multiplier_type) {
        bool arg_is_new = (arg_pos == passive_id<std::size_t>);

        if (arg_is_new) {
            auto& arg_buffer = this->buffers[buffer_id];
            if (arg_buffer.free_positions.empty()) {
                arg_pos = arg_buffer.size;
                ++arg_buffer.size;
                arg_buffer.values.resize(arg_buffer.values.size() + this->m_num_lanes);
            }
            else {
                arg_pos = arg_buffer.free_positions.back();
                arg_buffer.free_positions.pop_back();
            }
            // this is a new value, we NEED to override
            if (multiplier_type == mul_type::ANY) {
                mul_set(buffer_id, res_pos, arg_pos, multiplier);
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
                mul_add(buffer_id, res_pos, arg_pos, multiplier);
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
        auto buffer_id = static_cast<std::uint8_t>(std::distance(checkpoints_c.cbegin(), it) - 1);
        return { it == checkpoints_c.end(), buffer_id };
    };

    for (std::size_t op_idx = from; op_idx-- > to;) {
        std::size_t const first_multiplier_origin = multiplier_origin[(op_idx - to) * 2];
        std::size_t const second_multiplier_origin = multiplier_origin[((op_idx - to) * 2) + 1];
        if (first_multiplier_origin != passive_id<std::size_t> && second_multiplier_origin == passive_id<std::size_t>) {
            std::size_t const arg_id = multiplier_loc_from[first_multiplier_origin];
            std::size_t const res_id = op_idx;
            bool const keep_alive = multiplier_keep_alive[first_multiplier_origin];

            std::size_t& res_pos = this->node_location_on_buffer[res_id];
            std::size_t& arg_pos = this->node_location_on_buffer[arg_id];

            double const multiplier = buffer_multipliers.values[first_multiplier_origin];
            auto const multiplier_type = values_type[first_multiplier_origin];

            auto const arg_pos_data = get_loc(arg_id);
            bool arg_is_new = (arg_pos == passive_id<std::size_t>);
            bool arg_inplace = arg_is_new && std::get<0>(arg_pos_data);

            if (arg_inplace && !keep_alive) {
                if (multiplier_type != mul_type::ONE) {
                    // res id should now be arg id, avoiding a copy and a potential buffer increase
                    if (multiplier_type == mul_type::ANY) {
                        mul_inplace(res_pos, multiplier);
                    }
                    else if (multiplier_type == mul_type::MINUS_ONE) {
                        minus_inplace(res_pos);
                    }
                }
                std::swap(res_pos, arg_pos);
            }
            else {
                copy_mul(res_pos, arg_pos, std::get<1>(arg_pos_data), multiplier, multiplier_type);

                if (!keep_alive) {
                    buffer_free_positions.push_back(res_pos);
                    res_pos = passive_id<std::size_t>;
                }
            }
        }
        else if (first_multiplier_origin != passive_id<std::size_t>) {
            double const mul_1 = buffer_multipliers.values[first_multiplier_origin];
            double const mul_2 = buffer_multipliers.values[second_multiplier_origin];
            auto const mul_1_type = values_type[first_multiplier_origin];
            auto const mul_2_type = values_type[second_multiplier_origin];

            std::size_t const lhs_id = multiplier_loc_from[first_multiplier_origin];
            std::size_t const rhs_id = multiplier_loc_from[second_multiplier_origin];
            std::size_t const res_id = op_idx;

            std::size_t& res_pos = this->node_location_on_buffer[res_id];
            std::size_t& lhs_pos = this->node_location_on_buffer[lhs_id];
            std::size_t& rhs_pos = this->node_location_on_buffer[rhs_id];

            auto const lhs_pos_data = get_loc(lhs_id);
            auto const rhs_pos_data = get_loc(rhs_id);

            bool lhs_is_new = (lhs_pos == passive_id<std::size_t>);
            bool rhs_is_new = (rhs_pos == passive_id<std::size_t>);
            bool lhs_inplace = lhs_is_new && std::get<0>(lhs_pos_data);
            bool rhs_inplace = !lhs_inplace && rhs_is_new && std::get<0>(rhs_pos_data);

            if (!lhs_inplace) {
                copy_mul(res_pos, lhs_pos, std::get<1>(lhs_pos_data), mul_1, mul_1_type);
            }

            if (!rhs_inplace) {
                copy_mul(res_pos, rhs_pos, std::get<1>(rhs_pos_data), mul_2, mul_2_type);
            }

            if (lhs_inplace) {
                // res id should now be lhs id, avoiding a copy and a potential buffer increase
                if (mul_1_type != mul_type::ONE) {
                    // res id should now be arg id, avoiding a copy and a potential buffer increase
                    if (mul_1_type == mul_type::ANY) {
                        mul_inplace(res_pos, mul_1);
                    }
                    else if (mul_1_type == mul_type::MINUS_ONE) {
                        minus_inplace(res_pos);
                    }
                }

                std::swap(res_pos, lhs_pos);
            }
            else if (rhs_inplace) {
                // res id should now be rhs id, avoiding a copy and a potential buffer increase
                if (mul_2_type != mul_type::ONE) {
                    // res id should now be arg id, avoiding a copy and a potential buffer increase
                    if (mul_2_type == mul_type::ANY) {
                        mul_inplace(res_pos, mul_2);
                    }
                    else if (mul_2_type == mul_type::MINUS_ONE) {
                        minus_inplace(res_pos);
                    }
                }

                std::swap(res_pos, rhs_pos);
            }
            else {
                // don't forget to free res_id from the buffer!
                buffer_free_positions.push_back(res_pos);
                res_pos = passive_id<std::size_t>;
            }
        }
    }

    if constexpr (Reset) {
        data.ops.resize(to);
        data.vals.resize(val_idx_start);
        data.ids.resize(id_idx_start);
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

#endif // ADHOC_BACKPROPAGATOR1LOSSYCOMPRESSED_HPP
