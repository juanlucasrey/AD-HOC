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

#ifndef ADHOC_BACKPROPAGATOR1LOSSYPATHREUSE_HPP
#define ADHOC_BACKPROPAGATOR1LOSSYPATHREUSE_HPP

#include "hash.hpp"
#include "passive_id.hpp"
#include "position_impl.hpp"
#include "tape_converter.hpp"
#include "tape_data.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

#include <map>

namespace adhoc {

template<class Float, bool Vectorised = false>
class BackPropagatorLossyPathReuse {
  private:
    std::size_t m_num_lanes{ 1 };

    std::vector<std::size_t> node_location_on_buffer;

    std::vector<std::size_t> checkpoints{ 0 };
    std::vector<buffer_t> buffers{ buffer_t{} };

    std::vector<std::map<std::size_t, LossyTape> > hash_to_lossy_tape;

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

  public:
    explicit BackPropagatorLossyPathReuse() = default;

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

            auto it = std::upper_bound(this->checkpoints.cbegin(), this->checkpoints.cend(), var_id);
            return this->buffers[std::distance(this->checkpoints.cbegin(), it) - 1]
              .values[(var_pos * this->m_num_lanes) + lane];
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
        for (auto& b : this->buffers) {
            std::fill(b.values.begin(), b.values.end(), 0.0);
        }
    }

    auto size_of(bool capacity = false) const -> std::size_t
    {
        std::size_t size = 0;
        size += 1 * sizeof(std::size_t); // m_num_lanes
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
BackPropagatorLossyPathReuse<Float, Vectorised>::backpropagate_to(PositionImpl const& pos, TapeData& data)
{
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
        auto const h = hash(pos, data);

        bool const tape_exists = hash_to_lossy_tape.back().find(h) != hash_to_lossy_tape.back().end();
        auto& lossy_tape = hash_to_lossy_tape.back()[h];
        if (!tape_exists) {
            // we need to compute the lossy tape for this section.
            lossy_tape = convert_to_lossy_tape<Reset, Vectorised>(
              to, data, node_location_on_buffer, checkpoints, buffers, this->m_num_lanes);
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

#endif // ADHOC_BACKPROPAGATOR1LOSSYPATHREUSE_HPP
