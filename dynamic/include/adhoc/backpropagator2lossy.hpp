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

#ifndef ADHOC_BACKPROPAGATOR2LOSSY_HPP
#define ADHOC_BACKPROPAGATOR2LOSSY_HPP

#include "buffer_t.hpp"
#include "maps.hpp"
#include "passive_id.hpp"
#include "position_impl.hpp"
#include "tape_data.hpp"

#include <cmath>
#include <numbers>
#include <vector>

namespace adhoc {

template<class Float, MapType maptype, bool Vectorised = false>
class BackPropagator2Lossy {
  private:
    // for now
    static_assert(has_stable_references<maptype>(), "BackPropagator2Lossy requires a map type with stable references.");

    std::vector<map_t<maptype, std::size_t> > node_location_on_buffer;
    buffer_t<double, Vectorised> buffer;

  public:
    explicit BackPropagator2Lossy() = default;

    void set_checkpoint(std::size_t /* ops_size */) {}
    void set_lanes(std::size_t num_lanes) { this->buffer = buffer_t<double, Vectorised>{ num_lanes }; }
    auto get_lanes() const -> std::size_t { return this->buffer.lanes(); }
    void reserve_input(std::size_t count_registered) { this->buffer.reserve(this->buffer.size() + count_registered); }
    void reserve_output(std::size_t count_registered) { this->buffer.reserve(this->buffer.size() + count_registered); }
    void register_variable(std::size_t var_id)
    {
        this->node_location_on_buffer.resize(std::max(this->node_location_on_buffer.size(), var_id + 1));

        auto& var_pos = this->node_location_on_buffer[var_id];
        bool const has_id = var_pos.find(passive_id<std::size_t>) != var_pos.end();

        if (!has_id) {
            var_pos[passive_id<std::size_t>] = this->buffer.template get_new_loc<true>();
        }
    }
    void register_output_variable(std::size_t var_id, std::size_t /* ops_size */)
    {
        this->node_location_on_buffer.resize(std::max(this->node_location_on_buffer.size(), var_id + 1));

        auto& ders_var = this->node_location_on_buffer[var_id];
        bool const has_id = ders_var.find(passive_id<std::size_t>) != ders_var.end();

        if (!has_id) {
            ders_var[passive_id<std::size_t>] = this->buffer.template get_new_loc<true>();
        }
    }

    void set_derivative(std::size_t var_id, double deriv, std::size_t /* ops_size */, std::size_t lane = 0)
    {
        if (this->node_location_on_buffer.size() <= var_id) {
            // this derivative is not on buffer.
            // this is probably nor an input nor an output.
            throw;
        }

        auto& ders_var = this->node_location_on_buffer[var_id];
        bool const has_id = ders_var.find(passive_id<std::size_t>) != ders_var.end();
        if (!has_id) {
            // this derivative is not on buffer.
            // this is probably nor an input nor an output.
            throw;
        }

        std::size_t var_pos = ders_var[passive_id<std::size_t>];

        if constexpr (Vectorised) {
            auto val = this->buffer[var_pos];

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
            this->buffer[var_pos] = deriv;
        }
    }

    void set_derivative(std::size_t /* var_id1 */,
                        std::size_t /* var_id2 */,
                        double /* deriv */,
                        std::size_t /* lane */ = 0)
    {
        // TODO
        throw;
    }

    auto get_derivative(std::size_t var_id, std::size_t lane) const -> double
    {
        if (this->node_location_on_buffer.size() <= var_id) {
            return 0.;
        }

        auto& ders_var = this->node_location_on_buffer[var_id];
        bool const has_id = ders_var.find(passive_id<std::size_t>) != ders_var.end();
        if (!has_id) {
            return 0.;
        }

        std::size_t var_pos = ders_var.at(passive_id<std::size_t>);

        if constexpr (Vectorised) {
            auto const val = this->buffer[var_pos];

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
            return this->buffer[var_pos];
        }
    }

    auto get_derivative(std::size_t var_id1, std::size_t var_id2, std::size_t lane) const -> double
    {
        if (var_id1 < var_id2 && var_id2 != passive_id<std::size_t>) {
            std::swap(var_id1, var_id2);
        }

        if (this->node_location_on_buffer.size() <= var_id1) {
            return 0.;
        }

        auto& ders_var = this->node_location_on_buffer[var_id1];
        bool const has_id = ders_var.find(var_id2) != ders_var.end();
        if (!has_id) {
            return 0.;
        }

        std::size_t var_pos = ders_var.at(var_id2);

        if constexpr (Vectorised) {
            auto const val = this->buffer[var_pos];

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
            return this->buffer[var_pos];
        }
    }

    void clear() { this->node_location_on_buffer.clear(); }

    void zero_adjoints()
    {
        this->node_location_on_buffer.clear();
        this->buffer.zero();
    }

    auto size_of(bool capacity = false) const -> std::size_t
    {
        std::size_t size = 0;
        size += sizeof(std::size_t) * (capacity ? node_location_on_buffer.capacity() : node_location_on_buffer.size());
        size += this->buffer.size_of(capacity);
        return size;
    }

    template<bool Reset = false, bool ResetInPlace = false, bool Log = false>
    void backpropagate_to(PositionImpl const& pos, TapeData& data);
};

template<class Float, MapType maptype, bool Vectorised>
template<bool Reset, bool ResetInPlace, bool Log>
void
BackPropagator2Lossy<Float, maptype, Vectorised>::backpropagate_to(PositionImpl const& pos, TapeData& data)
{
    std::size_t to = pos.op_position;
    std::size_t from = data.next_id;

    const auto& ops = data.ops;
    const auto& vals = data.vals;
    const auto& ids = data.ids;

    this->node_location_on_buffer.resize(ops.size());

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

    auto mul_inplace = [&](std::size_t const pos, double const multiplier) {
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

    auto mul_add = [&](std::size_t const out_pos, std::size_t const in_pos, double const multiplier) {
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

    auto mul_set = [&](std::size_t const out_pos, std::size_t const in_pos, double const multiplier) {
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

    auto get_buffer_idx = [&](std::size_t id1, std::size_t id2) -> std::size_t& {
        if (id1 < id2 && id2 != passive_id<std::size_t>) {
            std::swap(id1, id2);
        }
        auto& locations_id1 = this->node_location_on_buffer[id1];
        auto it = locations_id1.find(id2);
        bool const arg_is_new = (it == locations_id1.end());
        std::size_t& arg_pos = locations_id1[id2];
        if (arg_is_new) {
            arg_pos = passive_id<std::size_t>;
        }
        return arg_pos;
    };

    auto update_univariate1 = [&](std::size_t arg_id, std::size_t res_id, double der_local_1) {
        auto& der_list = this->node_location_on_buffer[res_id];
        for (auto& der_pair : der_list) {
            std::size_t const der_id = der_pair.first;
            std::size_t& res_pos = der_pair.second;
            if (der_id == passive_id<std::size_t>) {
                std::size_t& id_pos = get_buffer_idx(arg_id, der_id);

                bool const id_inplace = (id_pos == passive_id<std::size_t>);

                if (id_inplace) {
                    mul_inplace(res_pos, der_local_1);
                    std::swap(id_pos, res_pos);
                }
                else {
                    mul_add(res_pos, id_pos, der_local_1);
                    // don't forget to free res_id from the buffer!
                    this->buffer.free_loc(res_pos);
                    res_pos = passive_id<std::size_t>;
                }
            }
            else if (der_id == res_id) {
                std::size_t& id_pos = get_buffer_idx(arg_id, arg_id);

                bool const id_inplace = (id_pos == passive_id<std::size_t>);

                if (id_inplace) {
                    mul_inplace(res_pos, der_local_1 * der_local_1);
                    std::swap(id_pos, res_pos);
                }
                else {
                    mul_add(res_pos, id_pos, der_local_1 * der_local_1);
                    // don't forget to free res_id from the buffer!
                    this->buffer.free_loc(res_pos);
                    res_pos = passive_id<std::size_t>;
                }
            }
            else {
                std::size_t& id_pos = get_buffer_idx(arg_id, der_id);

                bool const id_inplace = (id_pos == passive_id<std::size_t>);

                if (id_inplace) {
                    mul_inplace(res_pos, (arg_id == der_id ? 2. : 1.) * der_local_1);
                    std::swap(id_pos, res_pos);
                }
                else {
                    mul_add(res_pos, id_pos, (arg_id == der_id ? 2. : 1.) * der_local_1);
                    // don't forget to free res_id from the buffer!
                    this->buffer.free_loc(res_pos);
                    res_pos = passive_id<std::size_t>;
                }
            }
        }
    };

    auto update_univariate2 = [&](std::size_t arg_id, std::size_t res_id, double der_local_1, double der_local_2) {
        auto& der_list = this->node_location_on_buffer[res_id];
        for (auto& der_pair : der_list) {
            std::size_t const der_id = der_pair.first;
            std::size_t& res_pos = der_pair.second;
            if (der_id == passive_id<std::size_t>) {
                auto& locations_id = this->node_location_on_buffer[arg_id];
                auto const it1 = locations_id.find(passive_id<std::size_t>);
                auto const it2 = locations_id.find(arg_id);

                bool const id1_is_new = (it1 == locations_id.end());
                bool const id2_is_new = (it2 == locations_id.end());

                if (!id1_is_new) {
                    std::size_t const d1_pos = it1->second;
                    mul_add(res_pos, d1_pos, der_local_1);
                }

                if (!id2_is_new) {
                    std::size_t const d2_pos = it2->second;
                    mul_add(res_pos, d2_pos, der_local_2);
                }

                if (id1_is_new) {
                    if (id2_is_new) {
                        std::size_t const new_pos = this->buffer.get_new_loc();
                        locations_id[arg_id] = new_pos;
                        mul_set(res_pos, new_pos, der_local_2);
                    }

                    mul_inplace(res_pos, der_local_1);
                    locations_id[passive_id<std::size_t>] = res_pos;
                    res_pos = passive_id<std::size_t>;
                }
                else if (id2_is_new) {
                    mul_inplace(res_pos, der_local_2);
                    locations_id[arg_id] = res_pos;
                    res_pos = passive_id<std::size_t>;
                }
                else {
                    // don't forget to free res_id from the buffer!
                    this->buffer.free_loc(res_pos);
                    res_pos = passive_id<std::size_t>;
                }
            }
            else if (der_id == res_id) {
                std::size_t& id_pos = get_buffer_idx(arg_id, arg_id);

                bool const id_inplace = (id_pos == passive_id<std::size_t>);

                if (id_inplace) {
                    mul_inplace(res_pos, der_local_1 * der_local_1);
                    std::swap(id_pos, res_pos);
                }
                else {
                    mul_add(res_pos, id_pos, der_local_1 * der_local_1);
                    // don't forget to free res_id from the buffer!
                    this->buffer.free_loc(res_pos);
                    res_pos = passive_id<std::size_t>;
                }
            }
            else {
                std::size_t& id_pos = get_buffer_idx(arg_id, der_id);

                bool const id_inplace = (id_pos == passive_id<std::size_t>);

                if (id_inplace) {
                    mul_inplace(res_pos, (arg_id == der_id ? 2. : 1.) * der_local_1);
                    std::swap(id_pos, res_pos);
                }
                else {
                    mul_add(res_pos, id_pos, (arg_id == der_id ? 2. : 1.) * der_local_1);
                    // don't forget to free res_id from the buffer!
                    this->buffer.free_loc(res_pos);
                    res_pos = passive_id<std::size_t>;
                }
            }
        }
    };

    std::size_t val_idx = vals.size();
    std::size_t id_idx = ids.size();
    for (std::size_t op_idx = from; op_idx-- > to;) {
        OpCode const& op = ops[op_idx];
        bool const use_this_op = static_cast<bool>(this->node_location_on_buffer[op_idx].size());

        switch (op) {
            case OpCode::REG_INPUT: {
                id_idx -= 1;
                if constexpr (Reset) {
                    std::size_t const id = ids[id_idx];

                    auto& der_list = this->node_location_on_buffer[id];
                    for (auto& der_pair : der_list) {
                        std::size_t& pos = der_pair.second;
                        this->buffer.free_loc(pos);
                        pos = passive_id<std::size_t>;
                    }
                }
                break;
            }
            case OpCode::REG_OUTPUT: {
                id_idx -= 2;
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];

                    auto& der_list = this->node_location_on_buffer[res_id];
                    for (auto& der_pair : der_list) {
                        std::size_t const der_id = der_pair.first;
                        std::size_t& res_pos = der_pair.second;
                        std::size_t& arg_pos = get_buffer_idx(arg_id, der_id);

                        bool const arg_inplace = (arg_pos == passive_id<std::size_t>);

                        if constexpr (Reset) {
                            if (arg_inplace) {
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
                        }
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

                    auto& der_list = this->node_location_on_buffer[res_id];
                    for (auto& der_pair : der_list) {
                        std::size_t const der_id = der_pair.first;
                        std::size_t& res_pos = der_pair.second;

                        if (der_id == passive_id<std::size_t>) {
                            std::size_t& id1_pos = get_buffer_idx(lhs_id, passive_id<std::size_t>);
                            std::size_t& id2_pos = get_buffer_idx(rhs_id, passive_id<std::size_t>);

                            bool const id1_inplace = (id1_pos == passive_id<std::size_t>);
                            bool const id2_is_new = (id2_pos == passive_id<std::size_t>);
                            bool const id2_inplace = !id1_inplace && id2_is_new;

                            if (!id1_inplace) {
                                add(res_pos, id1_pos);
                            }

                            if (!id2_inplace) {
                                if (id2_is_new) {
                                    id2_pos = this->buffer.get_new_loc();
                                    copy(res_pos, id2_pos);
                                }
                                else {
                                    add(res_pos, id2_pos);
                                }
                            }

                            if (id1_inplace) {
                                std::swap(id1_pos, res_pos);
                            }
                            else if (id2_inplace) {
                                std::swap(id2_pos, res_pos);
                            }
                            else {
                                // don't forget to free res_id from the buffer!
                                this->buffer.free_loc(res_pos);
                                res_pos = passive_id<std::size_t>;
                            }
                        }
                        else if (der_id == res_id) {
                            bool const lhs_is_smaller = (lhs_id < rhs_id);
                            std::size_t const small_id = lhs_is_smaller ? lhs_id : rhs_id;
                            std::size_t const large_id = lhs_is_smaller ? rhs_id : lhs_id;
                            auto& locations_small_id = this->node_location_on_buffer[small_id];
                            auto& locations_large_id = this->node_location_on_buffer[large_id];

                            auto const it1 = locations_small_id.find(small_id);
                            auto const it2 = locations_large_id.find(large_id);
                            auto const it3 = locations_large_id.find(small_id);

                            bool const id1_is_new = (it1 == locations_small_id.end());
                            bool const id2_is_new = (it2 == locations_large_id.end());
                            bool const id3_is_new = (it3 == locations_large_id.end());

                            if (!id1_is_new) {
                                std::size_t const d1_pos = it1->second;
                                add(res_pos, d1_pos);
                            }

                            if (!id2_is_new) {
                                std::size_t const d2_pos = it2->second;
                                add(res_pos, d2_pos);
                            }

                            if (!id3_is_new) {
                                std::size_t const d3_pos = it3->second;
                                add(res_pos, d3_pos);
                            }

                            if (id1_is_new) {
                                if (id2_is_new) {
                                    std::size_t const new_pos = this->buffer.get_new_loc();
                                    locations_large_id[large_id] = new_pos;
                                    copy(res_pos, new_pos);
                                }

                                if (id3_is_new) {
                                    std::size_t const new_pos = this->buffer.get_new_loc();
                                    locations_large_id[small_id] = new_pos;
                                    copy(res_pos, new_pos);
                                }

                                locations_small_id[small_id] = res_pos;
                                res_pos = passive_id<std::size_t>;
                            }
                            else if (id2_is_new) {
                                if (id3_is_new) {
                                    std::size_t const new_pos = this->buffer.get_new_loc();
                                    locations_large_id[small_id] = new_pos;
                                    copy(res_pos, new_pos);
                                }

                                locations_large_id[large_id] = res_pos;
                                res_pos = passive_id<std::size_t>;
                            }
                            else if (id3_is_new) {
                                locations_large_id[small_id] = res_pos;
                                res_pos = passive_id<std::size_t>;
                            }
                            else {
                                // don't forget to free res_id from the buffer!
                                this->buffer.free_loc(res_pos);
                                res_pos = passive_id<std::size_t>;
                            }
                        }
                        else {
                            std::size_t const id1 = lhs_id < der_id ? lhs_id : der_id;
                            auto& map_id1 = lhs_id < der_id ? this->node_location_on_buffer[der_id]
                                                            : this->node_location_on_buffer[lhs_id];
                            auto const it1 = map_id1.find(id1);
                            bool const id1_is_new = (it1 == map_id1.end());

                            std::size_t const id2 = rhs_id < der_id ? rhs_id : der_id;
                            auto& map_id2 = rhs_id < der_id ? this->node_location_on_buffer[der_id]
                                                            : this->node_location_on_buffer[rhs_id];
                            auto const it2 = map_id2.find(id2);
                            bool const id2_is_new = (it2 == map_id2.end());

                            if (!id1_is_new) {
                                std::size_t const d1_pos = it1->second;
                                if (lhs_id == der_id) {
                                    mul_add(res_pos, d1_pos, 2.);
                                }
                                else {
                                    add(res_pos, d1_pos);
                                }
                            }

                            if (!id2_is_new) {
                                std::size_t const d2_pos = it2->second;
                                if (rhs_id == der_id) {
                                    mul_add(res_pos, d2_pos, 2.);
                                }
                                else {
                                    add(res_pos, d2_pos);
                                }
                            }

                            if (id1_is_new) {
                                if (id2_is_new) {
                                    std::size_t const new_pos = this->buffer.get_new_loc();
                                    map_id2[id2] = new_pos;

                                    if (rhs_id == der_id) {
                                        mul_set(res_pos, new_pos, 2.);
                                    }
                                    else {
                                        copy(res_pos, new_pos);
                                    }
                                }

                                if (lhs_id == der_id) {
                                    mul_inplace(res_pos, 2.);
                                }
                                map_id1[id1] = res_pos;
                                res_pos = passive_id<std::size_t>;
                            }
                            else if (id2_is_new) {
                                if (rhs_id == der_id) {
                                    mul_inplace(res_pos, 2.);
                                }
                                map_id2[id2] = res_pos;
                                res_pos = passive_id<std::size_t>;
                            }
                            else {
                                // don't forget to free res_id from the buffer!
                                this->buffer.free_loc(res_pos);
                                res_pos = passive_id<std::size_t>;
                            }
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

                    auto& der_list = this->node_location_on_buffer[res_id];
                    for (auto& der_pair : der_list) {
                        std::size_t const der_id = der_pair.first;
                        std::size_t& res_pos = der_pair.second;

                        if (der_id == passive_id<std::size_t>) {
                            std::size_t& id1_pos = get_buffer_idx(lhs_id, passive_id<std::size_t>);
                            std::size_t& id2_pos = get_buffer_idx(rhs_id, passive_id<std::size_t>);

                            bool const id1_inplace = (id1_pos == passive_id<std::size_t>);
                            bool const id2_is_new = (id2_pos == passive_id<std::size_t>);
                            bool const id2_inplace = !id1_inplace && id2_is_new;

                            if (!id1_inplace) {
                                add(res_pos, id1_pos);
                            }

                            if (!id2_inplace) {
                                if (id2_is_new) {
                                    id2_pos = this->buffer.get_new_loc();
                                    copy_minus(res_pos, id2_pos);
                                }
                                else {
                                    sub(res_pos, id2_pos);
                                }
                            }

                            if (id1_inplace) {
                                std::swap(id1_pos, res_pos);
                            }
                            else if (id2_inplace) {
                                minus_inplace(res_pos);
                                std::swap(id2_pos, res_pos);
                            }
                            else {
                                // don't forget to free res_id from the buffer!
                                this->buffer.free_loc(res_pos);
                                res_pos = passive_id<std::size_t>;
                            }
                        }
                        else if (der_id == res_id) {
                            bool const lhs_is_smaller = (lhs_id < rhs_id);
                            std::size_t const small_id = lhs_is_smaller ? lhs_id : rhs_id;
                            std::size_t const large_id = lhs_is_smaller ? rhs_id : lhs_id;
                            auto& locations_small_id = this->node_location_on_buffer[small_id];
                            auto& locations_large_id = this->node_location_on_buffer[large_id];

                            auto const it1 = locations_small_id.find(small_id);
                            auto const it2 = locations_large_id.find(large_id);
                            auto const it3 = locations_large_id.find(small_id);

                            bool const id1_is_new = (it1 == locations_small_id.end());
                            bool const id2_is_new = (it2 == locations_large_id.end());
                            bool const id3_is_new = (it3 == locations_large_id.end());

                            if (!id1_is_new) {
                                std::size_t const d1_pos = it1->second;
                                add(res_pos, d1_pos);
                            }

                            if (!id2_is_new) {
                                std::size_t const d2_pos = it2->second;
                                add(res_pos, d2_pos);
                            }

                            if (!id3_is_new) {
                                std::size_t const d3_pos = it3->second;
                                sub(res_pos, d3_pos);
                            }

                            if (id1_is_new) {
                                if (id2_is_new) {
                                    std::size_t const new_pos = this->buffer.get_new_loc();
                                    locations_large_id[large_id] = new_pos;
                                    copy(res_pos, new_pos);
                                }

                                if (id3_is_new) {
                                    std::size_t const new_pos = this->buffer.get_new_loc();
                                    locations_large_id[small_id] = new_pos;
                                    copy_minus(res_pos, new_pos);
                                }

                                locations_small_id[small_id] = res_pos;
                                res_pos = passive_id<std::size_t>;
                            }
                            else if (id2_is_new) {
                                if (id3_is_new) {
                                    std::size_t const new_pos = this->buffer.get_new_loc();
                                    locations_large_id[small_id] = new_pos;
                                    copy_minus(res_pos, new_pos);
                                }

                                locations_large_id[large_id] = res_pos;
                                res_pos = passive_id<std::size_t>;
                            }
                            else if (id3_is_new) {
                                minus_inplace(res_pos);
                                locations_large_id[large_id] = res_pos;
                                res_pos = passive_id<std::size_t>;
                            }
                            else {
                                // don't forget to free res_id from the buffer!
                                this->buffer.free_loc(res_pos);
                                res_pos = passive_id<std::size_t>;
                            }
                        }
                        else {
                            std::size_t const id1 = lhs_id < der_id ? lhs_id : der_id;
                            auto& map_id1 = lhs_id < der_id ? this->node_location_on_buffer[der_id]
                                                            : this->node_location_on_buffer[lhs_id];
                            auto const it1 = map_id1.find(id1);
                            bool const id1_is_new = (it1 == map_id1.end());

                            std::size_t const id2 = rhs_id < der_id ? rhs_id : der_id;
                            auto& map_id2 = rhs_id < der_id ? this->node_location_on_buffer[der_id]
                                                            : this->node_location_on_buffer[rhs_id];
                            auto const it2 = map_id2.find(id2);
                            bool const id2_is_new = (it2 == map_id2.end());

                            if (!id1_is_new) {
                                std::size_t const d1_pos = it1->second;
                                if (lhs_id == der_id) {
                                    mul_add(res_pos, d1_pos, 2.);
                                }
                                else {
                                    add(res_pos, d1_pos);
                                }
                            }

                            if (!id2_is_new) {
                                std::size_t const d2_pos = it2->second;
                                if (rhs_id == der_id) {
                                    mul_add(res_pos, d2_pos, -2.);
                                }
                                else {
                                    sub(res_pos, d2_pos);
                                }
                            }

                            if (id1_is_new) {
                                if (id2_is_new) {
                                    std::size_t const new_pos = this->buffer.get_new_loc();
                                    map_id2[id2] = new_pos;

                                    if (rhs_id == der_id) {
                                        mul_set(res_pos, new_pos, -2.);
                                    }
                                    else {
                                        copy_minus(res_pos, new_pos);
                                    }
                                }

                                if (lhs_id == der_id) {
                                    mul_inplace(res_pos, 2.);
                                }
                                map_id1[id1] = res_pos;
                                res_pos = passive_id<std::size_t>;
                            }
                            else if (id2_is_new) {
                                if (rhs_id == der_id) {
                                    mul_inplace(res_pos, -2.);
                                }
                                map_id2[id2] = res_pos;
                                res_pos = passive_id<std::size_t>;
                            }
                            else {
                                // don't forget to free res_id from the buffer!
                                this->buffer.free_loc(res_pos);
                                res_pos = passive_id<std::size_t>;
                            }
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

                    auto& der_list = this->node_location_on_buffer[res_id];
                    for (auto& der_pair : der_list) {
                        std::size_t const der_id = der_pair.first;
                        std::size_t& res_pos = der_pair.second;

                        if (der_id == passive_id<std::size_t>) {
                            auto& locations_lhs_id = this->node_location_on_buffer[lhs_id];
                            auto& locations_rhs_id = this->node_location_on_buffer[rhs_id];
                            auto& locations_large_id = (lhs_id < rhs_id) ? locations_rhs_id : locations_lhs_id;
                            std::size_t const small_id = (lhs_id < rhs_id) ? lhs_id : rhs_id;

                            // mixed derivative comes first because it implies a copy
                            auto const it1 = locations_large_id.find(small_id);
                            auto const it2 = locations_lhs_id.find(passive_id<std::size_t>);
                            auto const it3 = locations_rhs_id.find(passive_id<std::size_t>);

                            bool const id1_is_new = (it1 == locations_large_id.end());
                            bool const id2_is_new = (it2 == locations_lhs_id.end());
                            bool const id3_is_new = (it3 == locations_rhs_id.end());

                            if (!id1_is_new) {
                                std::size_t const d1_pos = it1->second;
                                add(res_pos, d1_pos);
                            }

                            if (!id2_is_new) {
                                std::size_t const d2_pos = it2->second;
                                mul_add(res_pos, d2_pos, rhs_val);
                            }

                            if (!id3_is_new) {
                                std::size_t const d3_pos = it3->second;
                                mul_add(res_pos, d3_pos, lhs_val);
                            }

                            if (id1_is_new) {
                                if (id2_is_new) {
                                    std::size_t const new_pos = this->buffer.get_new_loc();
                                    locations_lhs_id[passive_id<std::size_t>] = new_pos;
                                    mul_set(res_pos, new_pos, rhs_val);
                                }

                                if (id3_is_new) {
                                    std::size_t const new_pos = this->buffer.get_new_loc();
                                    locations_rhs_id[passive_id<std::size_t>] = new_pos;
                                    mul_set(res_pos, new_pos, lhs_val);
                                }

                                locations_large_id[small_id] = res_pos;
                                res_pos = passive_id<std::size_t>;
                            }
                            else if (id2_is_new) {
                                if (id3_is_new) {
                                    std::size_t const new_pos = this->buffer.get_new_loc();
                                    locations_rhs_id[passive_id<std::size_t>] = new_pos;
                                    mul_set(res_pos, new_pos, lhs_val);
                                }

                                mul_inplace(res_pos, rhs_val);
                                locations_lhs_id[passive_id<std::size_t>] = res_pos;
                                res_pos = passive_id<std::size_t>;
                            }
                            else if (id3_is_new) {
                                mul_inplace(res_pos, lhs_val);
                                locations_rhs_id[passive_id<std::size_t>] = res_pos;
                                res_pos = passive_id<std::size_t>;
                            }
                            else {
                                // don't forget to free res_id from the buffer!
                                this->buffer.free_loc(res_pos);
                                res_pos = passive_id<std::size_t>;
                            }
                        }
                        else if (der_id == res_id) {
                            auto& locations_lhs_id = this->node_location_on_buffer[lhs_id];
                            auto& locations_rhs_id = this->node_location_on_buffer[rhs_id];
                            auto& locations_large_id = (lhs_id < rhs_id) ? locations_rhs_id : locations_lhs_id;
                            std::size_t const small_id = (lhs_id < rhs_id) ? lhs_id : rhs_id;

                            auto const it1 = locations_large_id.find(small_id);
                            auto const it2 = locations_lhs_id.find(lhs_id);
                            auto const it3 = locations_rhs_id.find(rhs_id);

                            bool const id1_is_new = (it1 == locations_large_id.end());
                            bool const id2_is_new = (it2 == locations_lhs_id.end());
                            bool const id3_is_new = (it3 == locations_rhs_id.end());

                            if (!id1_is_new) {
                                std::size_t const d1_pos = it1->second;
                                mul_add(res_pos, d1_pos, rhs_val * lhs_val);
                            }

                            if (!id2_is_new) {
                                std::size_t const d2_pos = it2->second;
                                mul_add(res_pos, d2_pos, rhs_val * rhs_val);
                            }

                            if (!id3_is_new) {
                                std::size_t const d3_pos = it3->second;
                                mul_add(res_pos, d3_pos, lhs_val * lhs_val);
                            }

                            if (id1_is_new) {
                                if (id2_is_new) {
                                    std::size_t const new_pos = this->buffer.get_new_loc();
                                    locations_lhs_id[lhs_id] = new_pos;
                                    mul_set(res_pos, new_pos, rhs_val * rhs_val);
                                }

                                if (id3_is_new) {
                                    std::size_t const new_pos = this->buffer.get_new_loc();
                                    locations_rhs_id[rhs_id] = new_pos;
                                    mul_set(res_pos, new_pos, lhs_val * lhs_val);
                                }

                                mul_inplace(res_pos, rhs_val * lhs_val);
                                locations_large_id[small_id] = res_pos;
                                res_pos = passive_id<std::size_t>;
                            }
                            else if (id2_is_new) {
                                if (id3_is_new) {
                                    std::size_t const new_pos = this->buffer.get_new_loc();
                                    locations_rhs_id[rhs_id] = new_pos;
                                    mul_set(res_pos, new_pos, lhs_val * lhs_val);
                                }

                                mul_inplace(res_pos, rhs_val * rhs_val);
                                locations_lhs_id[lhs_id] = res_pos;
                                res_pos = passive_id<std::size_t>;
                            }
                            else if (id3_is_new) {
                                mul_inplace(res_pos, lhs_val * lhs_val);
                                locations_rhs_id[rhs_id] = res_pos;
                                res_pos = passive_id<std::size_t>;
                            }
                            else {
                                // don't forget to free res_id from the buffer!
                                this->buffer.free_loc(res_pos);
                                res_pos = passive_id<std::size_t>;
                            }
                        }
                        else {
                            std::size_t const id1 = lhs_id < der_id ? lhs_id : der_id;
                            auto& map_id1 = lhs_id < der_id ? this->node_location_on_buffer[der_id]
                                                            : this->node_location_on_buffer[lhs_id];
                            auto const it1 = map_id1.find(id1);
                            bool const id1_is_new = (it1 == map_id1.end());

                            std::size_t const id2 = rhs_id < der_id ? rhs_id : der_id;
                            auto& map_id2 = rhs_id < der_id ? this->node_location_on_buffer[der_id]
                                                            : this->node_location_on_buffer[rhs_id];
                            auto const it2 = map_id2.find(id2);
                            bool const id2_is_new = (it2 == map_id2.end());

                            if (!id1_is_new) {
                                std::size_t const d1_pos = it1->second;
                                mul_add(res_pos, d1_pos, (lhs_id == der_id ? 2. : 1.) * rhs_val);
                            }

                            if (!id2_is_new) {
                                std::size_t const d2_pos = it2->second;
                                mul_add(res_pos, d2_pos, (rhs_id == der_id ? 2. : 1.) * lhs_val);
                            }

                            if (id1_is_new) {
                                if (id2_is_new) {
                                    std::size_t const new_pos = this->buffer.get_new_loc();
                                    map_id2[id2] = new_pos;
                                    mul_set(res_pos, new_pos, (rhs_id == der_id ? 2. : 1.) * lhs_val);
                                }

                                mul_inplace(res_pos, (lhs_id == der_id ? 2. : 1.) * rhs_val);
                                map_id1[id1] = res_pos;
                                res_pos = passive_id<std::size_t>;
                            }
                            else if (id2_is_new) {
                                mul_inplace(res_pos, (rhs_id == der_id ? 2. : 1.) * lhs_val);
                                map_id2[id2] = res_pos;
                                res_pos = passive_id<std::size_t>;
                            }
                            else {
                                // don't forget to free res_id from the buffer!
                                this->buffer.free_loc(res_pos);
                                res_pos = passive_id<std::size_t>;
                            }
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

                    auto& der_list = this->node_location_on_buffer[res_id];
                    for (auto& der_pair : der_list) {
                        std::size_t const der_id = der_pair.first;
                        std::size_t& res_pos = der_pair.second;
                        if (der_id == passive_id<std::size_t>) {
                            std::size_t& id_pos = get_buffer_idx(arg_id, der_id);

                            bool const id_inplace = (id_pos == passive_id<std::size_t>);

                            if (id_inplace) {
                                std::swap(id_pos, res_pos);
                            }
                            else {
                                add(res_pos, id_pos);
                                // don't forget to free res_id from the buffer!
                                this->buffer.free_loc(res_pos);
                                res_pos = passive_id<std::size_t>;
                            }
                        }
                        else if (der_id == res_id) {
                            std::size_t& id_pos = get_buffer_idx(arg_id, arg_id);

                            bool const id_inplace = (id_pos == passive_id<std::size_t>);

                            if (id_inplace) {
                                std::swap(id_pos, res_pos);
                            }
                            else {
                                add(res_pos, id_pos);
                                // don't forget to free res_id from the buffer!
                                this->buffer.free_loc(res_pos);
                                res_pos = passive_id<std::size_t>;
                            }
                        }
                        else {
                            std::size_t& id_pos = get_buffer_idx(arg_id, der_id);

                            bool const id_inplace = (id_pos == passive_id<std::size_t>);

                            if (id_inplace) {
                                if (arg_id == der_id) {
                                    mul_inplace(res_pos, 2.);
                                }
                                std::swap(id_pos, res_pos);
                            }
                            else {
                                if (arg_id == der_id) {
                                    mul_add(res_pos, id_pos, 2.);
                                }
                                else {
                                    add(res_pos, id_pos);
                                }
                                // don't forget to free res_id from the buffer!
                                this->buffer.free_loc(res_pos);
                                res_pos = passive_id<std::size_t>;
                            }
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

                    auto& der_list = this->node_location_on_buffer[res_id];
                    for (auto& der_pair : der_list) {
                        std::size_t const der_id = der_pair.first;
                        std::size_t& res_pos = der_pair.second;
                        if (der_id == passive_id<std::size_t>) {
                            std::size_t& id_pos = get_buffer_idx(arg_id, der_id);

                            bool const id_inplace = (id_pos == passive_id<std::size_t>);

                            if (id_inplace) {
                                minus_inplace(res_pos);
                                std::swap(id_pos, res_pos);
                            }
                            else {
                                sub(res_pos, id_pos);
                                // don't forget to free res_id from the buffer!
                                this->buffer.free_loc(res_pos);
                                res_pos = passive_id<std::size_t>;
                            }
                        }
                        else if (der_id == res_id) {
                            std::size_t& id_pos = get_buffer_idx(arg_id, arg_id);

                            bool const id_inplace = (id_pos == passive_id<std::size_t>);
                            if (id_inplace) {
                                std::swap(id_pos, res_pos);
                            }
                            else {
                                add(res_pos, id_pos);
                                // don't forget to free res_id from the buffer!
                                this->buffer.free_loc(res_pos);
                                res_pos = passive_id<std::size_t>;
                            }
                        }
                        else {
                            std::size_t& id_pos = get_buffer_idx(arg_id, der_id);

                            bool const id_inplace = (id_pos == passive_id<std::size_t>);

                            if (id_inplace) {
                                if (arg_id == der_id) {
                                    mul_inplace(res_pos, -2.);
                                }
                                else {
                                    minus_inplace(res_pos);
                                }
                                std::swap(id_pos, res_pos);
                            }
                            else {
                                if (arg_id == der_id) {
                                    mul_add(res_pos, id_pos, -2.);
                                }
                                else {
                                    sub(res_pos, id_pos);
                                }
                                // don't forget to free res_id from the buffer!
                                this->buffer.free_loc(res_pos);
                                res_pos = passive_id<std::size_t>;
                            }
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
                    update_univariate1(arg_id, res_id, der_local_1);
                }
                break;
            }
            case OpCode::NORM: {
                val_idx -= 1;
                id_idx -= 2;
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    double const der_local_1 = 2.0 * vals[val_idx];
                    double const der_local_2 = 2.0;
                    update_univariate2(arg_id, res_id, der_local_1, der_local_2);
                }
                break;
            }
            case OpCode::INV: {
                val_idx -= 1;
                id_idx -= 2;
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    double const der_local_1 = -vals[val_idx] * vals[val_idx];
                    double const der_local_2 = -2.0 * der_local_1 * vals[val_idx];
                    update_univariate2(arg_id, res_id, der_local_1, der_local_2);
                }
                break;
            }
            case OpCode::ABS: {
                val_idx -= 1;
                id_idx -= 2;
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    double const der_local_1 = std::copysign(1.0, vals[val_idx]);
                    update_univariate1(arg_id, res_id, der_local_1);
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
                    update_univariate2(arg_id, res_id, der_local_1, der_local_2);
                }
                break;
            }
            case OpCode::LOG: {
                val_idx -= 1;
                id_idx -= 2;
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    double const der_local_1 = 1.0 / vals[val_idx];
                    double const der_local_2 = -der_local_1 * der_local_1;
                    update_univariate2(arg_id, res_id, der_local_1, der_local_2);
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
                    double const der_local_1 = std::exp(-vals[val_idx] * vals[val_idx]) * two_over_root_pi;
                    double const der_local_2 = -2. * vals[val_idx] * der_local_1;
                    update_univariate2(arg_id, res_id, der_local_1, der_local_2);
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
                    double const der_local_1 = std::exp(-vals[val_idx] * vals[val_idx]) * minus_two_over_root_pi;
                    double const der_local_2 = -2. * vals[val_idx] * der_local_1;
                    update_univariate2(arg_id, res_id, der_local_1, der_local_2);
                }
                break;
            }
            case OpCode::COS: {
                val_idx -= 2;
                id_idx -= 2;
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    double const der_local_1 = -std::sin(vals[val_idx]);
                    double const der_local_2 = -vals[val_idx + 1];
                    update_univariate2(arg_id, res_id, der_local_1, der_local_2);
                }
                break;
            }
            case OpCode::SQRT: {
                val_idx -= 2;
                id_idx -= 2;
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    double const one_over_in = 1. / vals[val_idx];
                    double const der_local_1 = 0.5 * vals[val_idx + 1] * one_over_in;
                    double const der_local_2 = -0.5 * der_local_1 * one_over_in;
                    update_univariate2(arg_id, res_id, der_local_1, der_local_2);
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
                    update_univariate2(arg_id, res_id, der_local_1, der_local_2);
                }
                break;
            }
        }
        if constexpr (Reset && ResetInPlace) {
            this->node_location_on_buffer.resize(this->node_location_on_buffer.size() - 1);
        }
    }

    if constexpr (Reset && !ResetInPlace) {
        this->node_location_on_buffer.resize(to);
    }

    if constexpr (Reset) {
        reset(pos, data);
    }
}

} // namespace adhoc

#endif // ADHOC_BACKPROPAGATOR2LOSSY_HPP
