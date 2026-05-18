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

#ifndef ADHOC_TAPE_CONVERTER_HPP
#define ADHOC_TAPE_CONVERTER_HPP

#include "passive_id.hpp"
#include "tape_data.hpp"

namespace adhoc {

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
    std::vector<std::uint8_t> on_which_buffer;
    std::vector<std::size_t> pos;
    std::vector<LossyOpCode> lossy_op;
    std::vector<bool> invert_mult;
    std::vector<bool> use_op;
    std::size_t buffer_size{ 0 };
};

struct buffer_t {
    std::vector<double> values;
    std::size_t size{ 0 };
    std::vector<std::size_t> free_positions;
};

template<bool Reset>
auto
convert_to_lossy_tape(std::size_t to,
                      TapeData const& data,
                      std::vector<std::size_t>& node_location_on_buffer,
                      std::vector<std::size_t> const& checkpoints,
                      std::vector<buffer_t>& buffers) -> LossyTape
{
    LossyTape result;

    std::size_t from = data.next_id;
    const auto& ops = data.ops;
    const auto& ids = data.ids;

    std::size_t id_idx = ids.size();

    auto& buffer_free_positions = buffers.back().free_positions;

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

    auto mul_inplace = [&](std::size_t const pos) {
        this_lossy_op.push_back(LossyOpCode::MUL_INPLACE);
        this_pos.push_back(pos);
    };

    auto mul_add = [&](std::size_t const out_pos, std::size_t const in_pos, std::uint8_t const which) {
        this_lossy_op.push_back(LossyOpCode::MUL_ADD);
        this_pos.push_back(out_pos);
        this_pos.push_back(in_pos);
        this_on_which_buffer.push_back(which);
    };

    auto mul_set = [&](std::size_t const out_pos, std::size_t const in_pos) {
        this_lossy_op.push_back(LossyOpCode::MUL_SET);
        this_pos.push_back(out_pos);
        this_pos.push_back(in_pos);
    };

    // output is: is_new, is_current, buffer_id, position
    auto get_loc = [&](std::size_t id) -> std::tuple<bool, std::uint8_t> {
        auto it = std::upper_bound(checkpoints.begin(), checkpoints.end(), id);
        auto buffer_id = static_cast<std::uint8_t>(std::distance(checkpoints.cbegin(), it) - 1);
        return { it == checkpoints.end(), buffer_id };
    };

    auto update_loc = [&](std::size_t& arg_pos, std::uint8_t buffer_id) {
        if (arg_pos == passive_id<std::size_t>) {
            auto& arg_buffer = buffers[buffer_id];
            if (arg_buffer.free_positions.empty()) {
                arg_pos = arg_buffer.size;
                ++arg_buffer.size;
            }
            else {
                arg_pos = arg_buffer.free_positions.back();
                arg_buffer.free_positions.pop_back();
            }
        }
    };

    for (std::size_t op_idx = from; op_idx-- > to;) {
        OpCode const& op = ops[op_idx];
        bool const use_this_op = node_location_on_buffer[op_idx] != passive_id<std::size_t>;

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

                    auto const arg_pos_data = get_loc(arg_id);

                    bool const arg_is_new = (arg_pos == passive_id<std::size_t>);
                    bool const arg_is_in_current_buffer = std::get<0>(arg_pos_data);
                    bool const arg_inplace = arg_is_new && arg_is_in_current_buffer;

                    if constexpr (Reset) {
                        if (arg_inplace) {
                            // res id should now be lhs id, avoiding a copy and a potential buffer increase
                            std::swap(arg_pos, res_pos);
                        }
                        else {
                            update_loc(arg_pos, std::get<1>(arg_pos_data));
                            add(res_pos, arg_pos, std::get<1>(arg_pos_data));
                            buffer_free_positions.push_back(res_pos);
                            res_pos = passive_id<std::size_t>;
                        }
                    }
                    else {
                        update_loc(arg_pos, std::get<1>(arg_pos_data));

                        // we don't reset so we don't free res_pos,
                        // it will be potentially used in the next operations as an active node location
                        if (arg_inplace) {
                            copy(res_pos, arg_pos);
                        }
                        else {
                            add(res_pos, arg_pos, std::get<1>(arg_pos_data));
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

                    std::size_t& res_pos = node_location_on_buffer[res_id];
                    std::size_t& lhs_pos = node_location_on_buffer[lhs_id];
                    std::size_t& rhs_pos = node_location_on_buffer[rhs_id];

                    auto const lhs_pos_data = get_loc(lhs_id);
                    auto const rhs_pos_data = get_loc(rhs_id);

                    bool const lhs_is_new = (lhs_pos == passive_id<std::size_t>);
                    bool const rhs_is_new = (rhs_pos == passive_id<std::size_t>);
                    bool const lhs_inplace = lhs_is_new && std::get<0>(lhs_pos_data);
                    bool const rhs_inplace_pre = rhs_is_new && std::get<0>(rhs_pos_data);
                    bool const rhs_inplace = !lhs_inplace && rhs_inplace_pre;

                    if (!lhs_inplace) {
                        update_loc(lhs_pos, std::get<1>(lhs_pos_data));
                        add(res_pos, lhs_pos, std::get<1>(lhs_pos_data));
                    }

                    if (!rhs_inplace) {
                        update_loc(rhs_pos, std::get<1>(rhs_pos_data));
                        if (rhs_inplace_pre) {
                            copy(res_pos, rhs_pos);
                        }
                        else {
                            add(res_pos, rhs_pos, std::get<1>(rhs_pos_data));
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

                    auto const lhs_pos_data = get_loc(lhs_id);
                    auto const rhs_pos_data = get_loc(rhs_id);

                    bool const lhs_is_new = (lhs_pos == passive_id<std::size_t>);
                    bool const rhs_is_new = (rhs_pos == passive_id<std::size_t>);
                    bool const lhs_inplace = lhs_is_new && std::get<0>(lhs_pos_data);
                    bool const rhs_inplace_pre = rhs_is_new && std::get<0>(rhs_pos_data);
                    bool const rhs_inplace = !lhs_inplace && rhs_inplace_pre;

                    if (!lhs_inplace) {
                        update_loc(lhs_pos, std::get<1>(lhs_pos_data));
                        add(res_pos, lhs_pos, std::get<1>(lhs_pos_data));
                    }

                    if (!rhs_inplace) {
                        update_loc(rhs_pos, std::get<1>(rhs_pos_data));
                        if (rhs_inplace_pre) {
                            copy_minus(res_pos, rhs_pos);
                        }
                        else {
                            sub(res_pos, rhs_pos, std::get<1>(rhs_pos_data));
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
                        buffer_free_positions.push_back(res_pos);
                        res_pos = passive_id<std::size_t>;
                    }
                }
                break;
            }
            case OpCode::MUL: {
                id_idx -= 3;
                result.use_op.push_back(use_this_op);
                if (use_this_op) {
                    std::size_t const lhs_id = ids[id_idx];
                    std::size_t const rhs_id = ids[id_idx + 1];
                    std::size_t const res_id = ids[id_idx + 2];

                    std::size_t& res_pos = node_location_on_buffer[res_id];
                    std::size_t& lhs_pos = node_location_on_buffer[lhs_id];
                    std::size_t& rhs_pos = node_location_on_buffer[rhs_id];

                    auto const lhs_pos_data = get_loc(lhs_id);
                    auto const rhs_pos_data = get_loc(rhs_id);

                    bool const lhs_is_new = (lhs_pos == passive_id<std::size_t>);
                    bool const rhs_is_new = (rhs_pos == passive_id<std::size_t>);
                    bool const lhs_inplace = lhs_is_new && std::get<0>(lhs_pos_data);
                    bool const rhs_inplace_pre = rhs_is_new && std::get<0>(rhs_pos_data);
                    bool const rhs_inplace = !lhs_inplace && rhs_inplace_pre;

                    // this is the only case when rhs_multiplier is used BEFORE lhs_multiplier.
                    result.invert_mult.push_back(lhs_inplace);

                    if (!lhs_inplace) {
                        update_loc(lhs_pos, std::get<1>(lhs_pos_data));
                        mul_add(res_pos, lhs_pos, std::get<1>(lhs_pos_data));
                    }

                    if (!rhs_inplace) {
                        update_loc(rhs_pos, std::get<1>(rhs_pos_data));
                        if (rhs_inplace_pre) {
                            mul_set(res_pos, rhs_pos);
                        }
                        else {
                            mul_add(res_pos, rhs_pos, std::get<1>(rhs_pos_data));
                        }
                    }

                    if (lhs_inplace) {
                        mul_inplace(res_pos);
                        // res id should now be lhs id, avoiding a copy and a potential buffer increase
                        std::swap(lhs_pos, res_pos);
                    }
                    else if (rhs_inplace) {
                        mul_inplace(res_pos);
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

                    auto const arg_pos_data = get_loc(arg_id);
                    bool const arg_is_new = (arg_pos == passive_id<std::size_t>);
                    bool const arg_inplace = arg_is_new && std::get<0>(arg_pos_data);

                    if (arg_inplace) {
                        // res id should now be lhs id, avoiding a copy and a potential buffer increase
                        std::swap(arg_pos, res_pos);
                    }
                    else {
                        update_loc(arg_pos, std::get<1>(arg_pos_data));
                        add(res_pos, arg_pos, std::get<1>(arg_pos_data));
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

                    auto const arg_pos_data = get_loc(arg_id);
                    bool const arg_is_new = (arg_pos == passive_id<std::size_t>);
                    bool const arg_inplace = arg_is_new && std::get<0>(arg_pos_data);

                    if (arg_inplace) {
                        // this is a subtraction, so we need to negate the value in the buffer
                        minus_inplace(res_pos);
                        // res id should now be lhs id, avoiding a copy and a potential buffer increase
                        std::swap(arg_pos, res_pos);
                    }
                    else {
                        update_loc(arg_pos, std::get<1>(arg_pos_data));
                        sub(res_pos, arg_pos, std::get<1>(arg_pos_data));
                        buffer_free_positions.push_back(res_pos);
                        res_pos = passive_id<std::size_t>;
                    }
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
                result.use_op.push_back(use_this_op);
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];

                    std::size_t& res_pos = node_location_on_buffer[res_id];
                    std::size_t& arg_pos = node_location_on_buffer[arg_id];

                    auto const arg_pos_data = get_loc(arg_id);
                    bool const arg_is_new = (arg_pos == passive_id<std::size_t>);
                    bool const arg_inplace = arg_is_new && std::get<0>(arg_pos_data);

                    if (arg_inplace) {
                        mul_inplace(res_pos);
                        // res id should now be arg id, avoiding a copy and a potential buffer increase
                        std::swap(arg_pos, res_pos);
                    }
                    else {
                        update_loc(arg_pos, std::get<1>(arg_pos_data));
                        mul_add(res_pos, arg_pos, std::get<1>(arg_pos_data));
                        buffer_free_positions.push_back(res_pos);
                        res_pos = passive_id<std::size_t>;
                    }
                }
                break;
            }
        }
    }
    result.buffer_size = buffers.back().size;

    return result;
}

} // namespace adhoc

#endif // ADHOC_TAPE_CONVERTER_HPP
