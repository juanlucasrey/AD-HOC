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

#include "buffer_t.hpp"
#include "hash.hpp"
#include "passive_id.hpp"
#include "position_impl.hpp"
#include "tape_data.hpp"

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <iterator>
#include <map>
#include <numbers>
#include <vector>

namespace adhoc {

template<std::floating_point Float, bool Vectorised = false>
class BackPropagatorLossyCompressedPathReuse {
  private:
    std::vector<std::size_t> node_location_on_buffer;
    std::vector<std::size_t> checkpoints{ 0 };
    std::vector<buffer_t<double, Vectorised> > buffers{ buffer_t<double, Vectorised>{} };

    // lossy tape
    enum class MultiplierOpCode : std::uint8_t {
        COPY,          // data[i0] = source
        COPY_MINUS,    // data[i0] = -source
        MINUS_INPLACE, // data[i0] = -data[i0]
        MUL_INPLACE,   // data[i0] *= source
        ADD_INTERNAL,  // data[i0] += data[i1]
        ADD_ONE,       // data[i0] += 1.0
        SUB_ONE,       // data[i0] -= 1.0
        MUL_INTERNAL,  // data[i0] *= data[i1]
    };

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
        std::vector<MultiplierOpCode> m_op;
        std::vector<std::size_t> m_pos;
        std::vector<bool> use_op;
        std::size_t mult_size{ 0 };

        std::vector<LossyOpCode> lossy_op;
        std::vector<std::uint8_t> on_which_buffer;
        std::vector<std::size_t> pos;
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
            this->buffers.push_back(buffer_t<double, Vectorised>{ this->get_lanes() });
            this->hash_to_lossy_tape.push_back(std::map<std::size_t, LossyTape>{});
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

template<std::floating_point Float, bool Vectorised>
template<bool Reset, bool ResetInPlace, bool Log, class TapeDataType>
void
BackPropagatorLossyCompressedPathReuse<Float, Vectorised>::backpropagate_to(PositionImpl const& pos, TapeDataType& data)
{
    class ValueFetcher {
      public:
        explicit ValueFetcher(TapeDataType const& data_in, LossyTape const& lossy_tape_in, PositionImpl const& pos)
          : data(data_in)
          , lossy_tape(lossy_tape_in)
          , op_idx(pos.op_position)
          , from(data.next_id)
          , val_idx(pos.val_position)
        {
        }

        auto operator()() -> double
        {
            if (use_cache) {
                use_cache = false;
                return cached_value;
            }

            std::optional<double> result;
            const auto& ops = data.ops;

            while (op_idx <= from && !result.has_value()) {
                OpCode const& op = ops[op_idx++];

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
                        bool const use_this_op = lossy_tape.use_op[use_op_idx++];
                        if (use_this_op) {
                            double const lhs_val = data.vals[val_idx];
                            double const rhs_val = data.vals[val_idx + 1];
                            use_cache = true;
                            cached_value = lhs_val;
                            result = rhs_val;
                        }
                        val_idx += 2;
                        break;
                    }
                    case OpCode::MUL_C: {
                        bool const use_this_op = lossy_tape.use_op[use_op_idx++];
                        if (use_this_op) {
                            double const der_local_1 = data.vals[val_idx];
                            result = der_local_1;
                        }
                        val_idx += 1;
                        break;
                    }
                    case OpCode::NORM: {
                        bool const use_this_op = lossy_tape.use_op[use_op_idx++];
                        if (use_this_op) {
                            double const val = data.vals[val_idx];
                            double const der_local_1 = 2.0 * val;
                            result = der_local_1;
                        }
                        val_idx += 1;
                        break;
                    }
                    case OpCode::INV: {
                        bool const use_this_op = lossy_tape.use_op[use_op_idx++];
                        if (use_this_op) {
                            double const val = data.vals[val_idx];
                            double const der_local_1 = -val * val;
                            result = der_local_1;
                        }
                        val_idx += 1;
                        break;
                    }
                    case OpCode::ABS: {
                        bool const use_this_op = lossy_tape.use_op[use_op_idx++];
                        if (use_this_op) {
                            double const val = data.vals[val_idx];
                            double const der_local_1 = std::copysign(1.0, val);
                            result = der_local_1;
                        }
                        val_idx += 1;
                        break;
                    }
                    case OpCode::EXP: {
                        bool const use_this_op = lossy_tape.use_op[use_op_idx++];
                        if (use_this_op) {
                            double const der_local_1 = data.vals[val_idx];
                            result = der_local_1;
                        }
                        val_idx += 1;
                        break;
                    }
                    case OpCode::LOG: {
                        bool const use_this_op = lossy_tape.use_op[use_op_idx++];
                        if (use_this_op) {
                            double const val = data.vals[val_idx];
                            double const der_local_1 = 1.0 / val;
                            result = der_local_1;
                        }
                        val_idx += 1;
                        break;
                    }
                    case OpCode::ERF: {
                        bool const use_this_op = lossy_tape.use_op[use_op_idx++];
                        if (use_this_op) {
                            constexpr double two_over_root_pi = 2. * std::numbers::inv_sqrtpi_v<double>;
                            double const val = data.vals[val_idx];
                            double const der_local_1 = std::exp(-val * val) * two_over_root_pi;
                            result = der_local_1;
                        }
                        val_idx += 1;
                        break;
                    }
                    case OpCode::ERFC: {
                        bool const use_this_op = lossy_tape.use_op[use_op_idx++];
                        if (use_this_op) {
                            constexpr double minus_two_over_root_pi = -2. * std::numbers::inv_sqrtpi_v<double>;
                            double const val = data.vals[val_idx];
                            double const der_local_1 = std::exp(-val * val) * minus_two_over_root_pi;
                            result = der_local_1;
                        }
                        val_idx += 1;
                        break;
                    }
                    case OpCode::COS: {
                        bool const use_this_op = lossy_tape.use_op[use_op_idx++];
                        if (use_this_op) {
                            double const val = data.vals[val_idx];
                            double const der_local_1 = -std::sin(val);
                            result = der_local_1;
                        }
                        val_idx += 2;
                        break;
                    }
                    case OpCode::SQRT: {
                        bool const use_this_op = lossy_tape.use_op[use_op_idx++];
                        if (use_this_op) {
                            double const val1 = data.vals[val_idx];
                            double const val2 = data.vals[val_idx + 1];
                            double const one_over_in = 1. / val1;
                            double const der_local_1 = 0.5 * val2 * one_over_in;
                            result = der_local_1;
                        }
                        val_idx += 2;
                        break;
                    }
                    case OpCode::POW_C: {
                        bool const use_this_op = lossy_tape.use_op[use_op_idx++];
                        if (use_this_op) {
                            double const lhs_arg = data.vals[val_idx];
                            double const rhs_arg = data.vals[val_idx + 1];
                            double const der_local_1 = rhs_arg != 0.0 ? rhs_arg * std::pow(lhs_arg, rhs_arg - 1.) : 0.0;
                            result = der_local_1;
                        }
                        val_idx += 2;
                        break;
                    }
                }
            }

            if (!result.has_value()) {
                throw;
            }

            return result.value();
        }

      private:
        TapeDataType const& data;
        LossyTape const& lossy_tape;

        std::size_t op_idx{ 0 };
        std::size_t from{ 0 };

        std::size_t val_idx{ 0 };
        std::size_t use_op_idx{ 0 };

        bool use_cache{ false };
        double cached_value{ 0. };
    };

    auto convert_to_lossy_tape = [&](PositionImpl const& pos, TapeDataType const& data) -> LossyTape {
        std::size_t const to = pos.op_position;
        LossyTape result;

        std::size_t from = data.next_id;
        const auto& ops = data.ops;
        const auto& ids = data.ids;

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
                case OpCode::SUB_C: {
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
                case OpCode::SUB: {
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
        std::size_t const id_idx_start = id_idx;

        std::vector<std::size_t> multiplier_origin((from - to) * 2, passive_id<std::size_t>);

        // LOOP 2: forward, to calculate multipliers after compressing induced paths
        buffer_t buffer_multipliers_values;

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

        auto copy_m = [&](std::size_t const pos) {
            result.m_op.push_back(MultiplierOpCode::COPY);
            result.m_pos.push_back(pos);
        };

        auto copy_minus_m = [&](std::size_t const pos) {
            result.m_op.push_back(MultiplierOpCode::COPY_MINUS);
            result.m_pos.push_back(pos);
        };

        auto minus_inplace_m = [&](std::size_t const pos) {
            result.m_op.push_back(MultiplierOpCode::MINUS_INPLACE);
            result.m_pos.push_back(pos);
        };

        auto mul_inplace_m = [&](std::size_t const pos) {
            result.m_op.push_back(MultiplierOpCode::MUL_INPLACE);
            result.m_pos.push_back(pos);
        };

        auto add_internal = [&](std::size_t const pos1, std::size_t const pos2) {
            result.m_op.push_back(MultiplierOpCode::ADD_INTERNAL);
            result.m_pos.push_back(pos1);
            result.m_pos.push_back(pos2);
        };

        auto add_one = [&](std::size_t const pos1) {
            result.m_op.push_back(MultiplierOpCode::ADD_ONE);
            result.m_pos.push_back(pos1);
        };

        auto sub_one = [&](std::size_t const pos1) {
            result.m_op.push_back(MultiplierOpCode::SUB_ONE);
            result.m_pos.push_back(pos1);
        };

        auto mul_internal = [&](std::size_t const pos1, std::size_t const pos2) {
            result.m_op.push_back(MultiplierOpCode::MUL_INTERNAL);
            result.m_pos.push_back(pos1);
            result.m_pos.push_back(pos2);
        };

        auto multiplier_set_incoming = [&](std::size_t const pos) {
            auto const real_pos = buffer_multipliers[pos].position;
            copy_m(real_pos);
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
                                                                             std::size_t& mult_origin_res) -> bool {
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
                            copy_m(buffer_origin);
                        }
                        else if (origin_type == mul_type::MINUS_ONE) {
                            copy_minus_m(buffer_origin);
                        }
                        origin_type = mul_type::ANY;
                    }
                    else {
                        mul_inplace_m(buffer_origin);
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
                    multiplier_set_incoming(mult_origin_res);
                }
                auto& node_data = buffer_multipliers[mult_origin_res];
                node_data.loc_from = arg_id;
                node_data.keep_alive = KeepAlive;
            }

            return arg_is_induced_path;
        };

        auto update_univariate = [&]<mul_type M, bool KeepAlive = false>(std::size_t const arg_id,
                                                                         std::size_t const res_id) {
            auto& mult_origin_res = multiplier_origin[(res_id - to) * 2];
            update_univariate_aux.template operator()<M, KeepAlive>(arg_id, mult_origin_res);
        };

        auto update_bivariate =
          [&]<mul_type M1, mul_type M2>(std::size_t const lhs_id, std::size_t const rhs_id, std::size_t const res_id) {
              std::size_t multiplier_loc_lhs = passive_id<std::size_t>;
              bool const lhs_is_induced_path =
                update_univariate_aux.template operator()<M1, false>(lhs_id, multiplier_loc_lhs);

              std::size_t multiplier_loc_rhs = passive_id<std::size_t>;
              bool const rhs_is_induced_path =
                update_univariate_aux.template operator()<M2, false>(rhs_id, multiplier_loc_rhs);

              bool const has_induced_path = lhs_is_induced_path || rhs_is_induced_path;
              // there is potential for a bivariate operator on the same argument, we need to check if
              // this is the case and update the multiplier if so
              auto& node_data_lhs = buffer_multipliers[multiplier_loc_lhs];
              auto& node_data_rhs = buffer_multipliers[multiplier_loc_rhs];
              bool const bivariate_consolidate_this =
                has_induced_path && (node_data_lhs.loc_from == node_data_rhs.loc_from);

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
              }
          };

        id_idx = id_idx_start;

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
                    result.use_op.push_back(use_this_op);
                    if (use_this_op) {
                        std::size_t const lhs_id = ids[id_idx];
                        std::size_t const rhs_id = ids[id_idx + 1];
                        std::size_t const res_id = ids[id_idx + 2];
                        update_bivariate.template operator()<mul_type::ANY, mul_type::ANY>(lhs_id, rhs_id, res_id);
                    }
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
                    result.use_op.push_back(use_this_op);
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];
                        update_univariate.template operator()<mul_type::ANY>(arg_id, res_id);
                    }
                    id_idx += 2;
                    break;
                }
            }
        }
        result.mult_size = buffer_multipliers_values.size();

        // LOOP 3: backward, to calculate derivatives on multiple lanes
        auto& buffer_back = this->buffers.back();
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

        auto mul_inplace = [&](std::size_t const pos, std::size_t const multiplier_id) {
            this_lossy_op.push_back(LossyOpCode::MUL_INPLACE);
            this_pos.push_back(pos);
            this_pos.push_back(multiplier_id);
        };

        auto mul_add = [&](std::size_t const out_pos,
                           std::size_t const in_pos,
                           std::uint8_t const which,
                           std::size_t const multiplier_id) {
            this_lossy_op.push_back(LossyOpCode::MUL_ADD);
            this_pos.push_back(out_pos);
            this_pos.push_back(in_pos);
            this_pos.push_back(multiplier_id);
            this_on_which_buffer.push_back(which);
        };

        auto mul_set = [&](std::size_t const out_pos, std::size_t const in_pos, std::size_t const multiplier_id) {
            this_lossy_op.push_back(LossyOpCode::MUL_SET);
            this_pos.push_back(out_pos);
            this_pos.push_back(in_pos);
            this_pos.push_back(multiplier_id);
        };

        auto reset_loc = [&](std::size_t arg_pos, std::uint8_t buffer_id) {
            if (buffers[buffer_id].allocated_size() > (arg_pos)) {
                if constexpr (Vectorised) {
                    auto dest = buffers[buffer_id][arg_pos];
#pragma omp simd
                    for (std::size_t i = 0; i < dest.size(); ++i) {
                        dest[i] = 0.;
                    }
                }
                else {
                    buffers[buffer_id][arg_pos] = 0.;
                }
            }
        };

        auto copy_mul = [&](std::size_t res_pos,
                            std::size_t& arg_pos,
                            bool const in_current_buffer,
                            std::uint8_t const buffer_id,
                            std::size_t const multiplier_id) {
            auto const& node_data = buffer_multipliers[multiplier_id];
            auto const multiplier_type = node_data.value_type;
            bool const arg_is_new = (arg_pos == passive_id<std::size_t>);

            if (arg_is_new && in_current_buffer) {
                arg_pos = buffer_back.template get_new_loc<false>();
                // this is a new value, we NEED to override
                if (multiplier_type == mul_type::ANY) {
                    auto const position_buffer = node_data.position;
                    mul_set(res_pos, arg_pos, position_buffer);
                }
                else if (multiplier_type == mul_type::ONE) {
                    copy(res_pos, arg_pos);
                }
                else if (multiplier_type == mul_type::MINUS_ONE) {
                    copy_minus(res_pos, arg_pos);
                }
            }
            else {
                if (arg_is_new) {
                    auto& arg_buffer = this->buffers[buffer_id];
                    arg_pos = arg_buffer.template get_new_loc<false>();
                    reset_loc(arg_pos, buffer_id);
                }

                if (multiplier_type == mul_type::ANY) {
                    auto const position_buffer = node_data.position;
                    mul_add(res_pos, arg_pos, buffer_id, position_buffer);
                }
                else if (multiplier_type == mul_type::ONE) {
                    add(res_pos, arg_pos, buffer_id);
                }
                else if (multiplier_type == mul_type::MINUS_ONE) {
                    sub(res_pos, arg_pos, buffer_id);
                }
            }
        };

        auto get_loc = [&](std::size_t id) -> std::tuple<bool, std::uint8_t> {
            auto it = std::upper_bound(this->checkpoints.begin(), this->checkpoints.end(), id);
            auto const buffer_id = static_cast<std::uint8_t>(std::distance(this->checkpoints.begin(), it) - 1);
            return { it == this->checkpoints.end(), buffer_id };
        };

        for (std::size_t op_idx = from; op_idx-- > to;) {
            std::size_t const first_multiplier_origin = multiplier_origin[(op_idx - to) * 2];
            std::size_t const second_multiplier_origin = multiplier_origin[((op_idx - to) * 2) + 1];
            if (first_multiplier_origin != passive_id<std::size_t> &&
                second_multiplier_origin == passive_id<std::size_t>) {
                auto const& node_data_first = buffer_multipliers[first_multiplier_origin];
                std::size_t const arg_id = node_data_first.loc_from;
                std::size_t const res_id = op_idx;
                bool const keep_alive = node_data_first.keep_alive;

                std::size_t& res_pos = this->node_location_on_buffer[res_id];
                std::size_t& arg_pos = this->node_location_on_buffer[arg_id];

                auto const arg_pos_data = get_loc(arg_id);
                bool const arg_is_new = (arg_pos == passive_id<std::size_t>);
                bool const arg_is_in_current_buffer = std::get<0>(arg_pos_data);
                bool const arg_inplace = arg_is_new && arg_is_in_current_buffer;

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
                    copy_mul(
                      res_pos, arg_pos, arg_is_in_current_buffer, std::get<1>(arg_pos_data), first_multiplier_origin);

                    if (!keep_alive) {
                        buffer_back.free_loc(res_pos);
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
                bool const lhs_is_in_current_buffer = std::get<0>(lhs_pos_data);
                bool const rhs_is_in_current_buffer = std::get<0>(rhs_pos_data);
                bool const lhs_inplace = lhs_is_new && lhs_is_in_current_buffer;
                bool const rhs_inplace = !lhs_inplace && rhs_is_new && rhs_is_in_current_buffer;

                if (!lhs_inplace) {
                    copy_mul(
                      res_pos, lhs_pos, lhs_is_in_current_buffer, std::get<1>(lhs_pos_data), first_multiplier_origin);
                }

                if (!rhs_inplace) {
                    copy_mul(
                      res_pos, rhs_pos, rhs_is_in_current_buffer, std::get<1>(rhs_pos_data), second_multiplier_origin);
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
                    buffer_back.free_loc(res_pos);
                    res_pos = passive_id<std::size_t>;
                }
            }
        }

        result.buffer_size = buffer_back.size();

        for (auto& b : this->buffers) {
            b.allocate();
        }

        return result;
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
    if (this->buffers.size() > 1) {
        // if (true) {
        // this->hash_to_lossy_tape.push_back(std::map<std::size_t, LossyTape>{});

        auto const h = hash(pos, data);

        bool const tape_exists = hash_to_lossy_tape.back().find(h) != hash_to_lossy_tape.back().end();
        auto& lossy_tape = hash_to_lossy_tape.back()[h];
        if (!tape_exists) {
            // we need to compute the lossy tape for this section.
            lossy_tape = convert_to_lossy_tape(pos, data);
        }

        // loop 1: calculate multipliers
        buffers.back().resize(lossy_tape.buffer_size);
        ValueFetcher value_fetcher(data, lossy_tape, pos);
        std::vector<double> buffer_multipliers_values;
        buffer_multipliers_values.resize(lossy_tape.mult_size, 0.);

        std::size_t pos_idx = 0;
        for (std::size_t i = 0; i < lossy_tape.m_op.size(); ++i) {
            auto const& op = lossy_tape.m_op[i];
            switch (op) {
                case MultiplierOpCode::COPY: {
                    std::size_t const pos_dest = lossy_tape.m_pos[pos_idx++];
                    auto const value = value_fetcher();
                    buffer_multipliers_values[pos_dest] = value;
                    break;
                }
                case MultiplierOpCode::COPY_MINUS: {
                    std::size_t const pos_dest = lossy_tape.m_pos[pos_idx++];
                    auto const value = value_fetcher();
                    buffer_multipliers_values[pos_dest] = -value;
                    break;
                }
                case MultiplierOpCode::MINUS_INPLACE: {
                    std::size_t const pos_dest = lossy_tape.m_pos[pos_idx++];
                    buffer_multipliers_values[pos_dest] = -buffer_multipliers_values[pos_dest];
                    break;
                }
                case MultiplierOpCode::MUL_INPLACE: {
                    std::size_t const pos_dest = lossy_tape.m_pos[pos_idx++];
                    auto const value = value_fetcher();
                    buffer_multipliers_values[pos_dest] *= value;
                    break;
                }
                case MultiplierOpCode::ADD_INTERNAL: {
                    std::size_t const pos_dest = lossy_tape.m_pos[pos_idx++];
                    std::size_t const pos_src = lossy_tape.m_pos[pos_idx++];
                    buffer_multipliers_values[pos_dest] += buffer_multipliers_values[pos_src];
                    break;
                }
                case MultiplierOpCode::ADD_ONE: {
                    std::size_t const pos_dest = lossy_tape.m_pos[pos_idx++];
                    buffer_multipliers_values[pos_dest] += 1.;
                    break;
                }
                case MultiplierOpCode::SUB_ONE: {
                    std::size_t const pos_dest = lossy_tape.m_pos[pos_idx++];
                    buffer_multipliers_values[pos_dest] -= 1.;
                    break;
                }
                case MultiplierOpCode::MUL_INTERNAL: {
                    std::size_t const pos_dest = lossy_tape.m_pos[pos_idx++];
                    std::size_t const pos_src = lossy_tape.m_pos[pos_idx++];
                    buffer_multipliers_values[pos_dest] *= buffer_multipliers_values[pos_src];
                    break;
                }
            }
        }

        // loop 2: propagate
        this->buffers.back().resize(lossy_tape.buffer_size);
        auto& buffer_vals = this->buffers.back();

        auto const& lossy_op = lossy_tape.lossy_op;
        auto const& on_which_buffer = lossy_tape.on_which_buffer;
        auto const& ltpos = lossy_tape.pos;

        std::size_t id_idx_l = 0;
        std::size_t on_which_buffer_idx = 0;

        for (std::size_t i = 0; i < lossy_op.size(); ++i) {
            auto const& op = lossy_op[i];
            switch (op) {
                case LossyOpCode::COPY: {
                    std::size_t const out_pos = ltpos[id_idx_l++];
                    std::size_t const in_pos = ltpos[id_idx_l++];
                    if constexpr (Vectorised) {
                        auto const src = buffer_vals[out_pos];
                        auto dest = buffer_vals[in_pos];
#pragma omp simd
                        for (std::size_t i = 0; i < dest.size(); ++i) {
                            dest[i] = src[i];
                        }
                    }
                    else {
                        buffer_vals[in_pos] = buffer_vals[out_pos];
                    }
                    break;
                }
                case LossyOpCode::COPY_MINUS: {
                    std::size_t const out_pos = ltpos[id_idx_l++];
                    std::size_t const in_pos = ltpos[id_idx_l++];
                    if constexpr (Vectorised) {
                        auto const src = buffer_vals[out_pos];
                        auto dest = buffer_vals[in_pos];
#pragma omp simd
                        for (std::size_t i = 0; i < dest.size(); ++i) {
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
                    std::size_t const out_pos = ltpos[id_idx_l++];
                    std::size_t const in_pos = ltpos[id_idx_l++];
                    if constexpr (Vectorised) {
                        auto const src = buffer_vals[out_pos];
                        auto dest = buffers[which][in_pos];
#pragma omp simd
                        for (std::size_t i = 0; i < dest.size(); ++i) {
                            dest[i] += src[i];
                        }
                    }
                    else {
                        buffers[which][in_pos] += buffer_vals[out_pos];
                    }
                    break;
                }
                case LossyOpCode::SUB: {
                    std::uint8_t const which = on_which_buffer[on_which_buffer_idx++];
                    std::size_t const out_pos = ltpos[id_idx_l++];
                    std::size_t const in_pos = ltpos[id_idx_l++];
                    if constexpr (Vectorised) {
                        auto const src = buffer_vals[out_pos];
                        auto dest = buffers[which][in_pos];
#pragma omp simd
                        for (std::size_t i = 0; i < dest.size(); ++i) {
                            dest[i] -= src[i];
                        }
                    }
                    else {
                        buffers[which][in_pos] -= buffer_vals[out_pos];
                    }
                    break;
                }
                case LossyOpCode::MINUS_INPLACE: {
                    std::size_t const out_pos = ltpos[id_idx_l++];
                    if constexpr (Vectorised) {
                        auto dest = buffer_vals[out_pos];
#pragma omp simd
                        for (std::size_t i = 0; i < dest.size(); ++i) {
                            dest[i] = -dest[i];
                        }
                    }
                    else {
                        buffer_vals[out_pos] = -buffer_vals[out_pos];
                    }
                    break;
                }
                case LossyOpCode::MUL_INPLACE: {
                    std::size_t const out_pos = ltpos[id_idx_l++];
                    std::size_t const mul_pos = ltpos[id_idx_l++];
                    double const multiplier = buffer_multipliers_values[mul_pos];
                    if constexpr (Vectorised) {
                        auto dest = buffer_vals[out_pos];
#pragma omp simd
                        for (std::size_t i = 0; i < dest.size(); ++i) {
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
                    std::size_t const out_pos = ltpos[id_idx_l++];
                    std::size_t const in_pos = ltpos[id_idx_l++];
                    std::size_t const mul_pos = ltpos[id_idx_l++];
                    double const multiplier = buffer_multipliers_values[mul_pos];
                    if constexpr (Vectorised) {
                        auto const src = buffer_vals[out_pos];
                        auto dest = buffers[which][in_pos];
#pragma omp simd
                        for (std::size_t i = 0; i < dest.size(); ++i) {
                            dest[i] += src[i] * multiplier;
                        }
                    }
                    else {
                        buffers[which][in_pos] += buffer_vals[out_pos] * multiplier;
                    }
                    break;
                }
                case LossyOpCode::MUL_SET: {
                    std::size_t const out_pos = ltpos[id_idx_l++];
                    std::size_t const in_pos = ltpos[id_idx_l++];
                    std::size_t const mul_pos = ltpos[id_idx_l++];
                    double const multiplier = buffer_multipliers_values[mul_pos];
                    if constexpr (Vectorised) {
                        auto const src = buffer_vals[out_pos];
                        auto dest = buffer_vals[in_pos];
#pragma omp simd
                        for (std::size_t i = 0; i < dest.size(); ++i) {
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

        // LOOP 1: backward, to count number of dependents for each node and detect which nodes are active
        std::vector<std::size_t> number_dependents(from - to);
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

        auto copy_m = [&](std::size_t const pos, double const multiplier) {
            buffer_multipliers_values[pos] = multiplier;
        };

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
            bool const bivariate_consolidate_this =
              has_induced_path && (node_data_lhs.loc_from == node_data_rhs.loc_from);

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
                        double const val = vals[val_idx];
                        double const multiplier = 2.0 * val;
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
                        double const val = vals[val_idx];
                        double const multiplier = -val * val;
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
                        double const val = vals[val_idx];
                        double const multiplier = std::copysign(1.0, val);
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
                        double const val = vals[val_idx];
                        double const multiplier = 1.0 / val;
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
                        double const val = vals[val_idx];
                        double const multiplier = std::exp(-val * val) * two_over_root_pi;
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
                        double const val = vals[val_idx];
                        double const multiplier = std::exp(-val * val) * minus_two_over_root_pi;
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
                        double const val = vals[val_idx];
                        double const multiplier = -std::sin(val);
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
                        double const val1 = vals[val_idx];
                        double const val2 = vals[val_idx + 1];
                        double const one_over_in = 1. / val1;
                        double const multiplier = 0.5 * val2 * one_over_in;
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
            if (first_multiplier_origin != passive_id<std::size_t> &&
                second_multiplier_origin == passive_id<std::size_t>) {
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
    }

    if constexpr (Reset) {
        if (to == this->checkpoints.back()) {
            this->buffers.back() = buffer_t<double, Vectorised>{ this->get_lanes() };
        }

        reset(pos, data);
        this->node_location_on_buffer.resize(pos.op_position);
    }
}

} // namespace adhoc

#endif // ADHOC_BACKPROPAGATOR1LOSSYCOMPRESSEDPATHREUSE_HPP
