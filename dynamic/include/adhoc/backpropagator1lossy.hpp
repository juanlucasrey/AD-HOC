#ifndef ADHOC_BACKPROPAGATOR1LOSSY_HPP
#define ADHOC_BACKPROPAGATOR1LOSSY_HPP

#include "passive_id.hpp"
#include "tape_data.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

namespace adhoc {

template<class Float, bool Vectorised = false>
class BackPropagatorLossy {
  private:
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

    std::size_t m_num_lanes{ 1 };

    std::vector<std::size_t> node_location_on_buffer;

    // lossy tape
    std::vector<std::uint8_t> on_which_buffer;
    std::vector<std::size_t> pos;
    std::vector<LossyOpCode> lossy_op;
    std::vector<double> values;

    std::size_t from_prev{ 0 };
    std::size_t to_prev{ 0 };

    struct buffer_t {
        std::vector<double> values;
        std::size_t size{ 0 };
        std::vector<std::size_t> free_positions;
    };

    std::vector<std::size_t> checkpoints{ 0 };
    std::vector<buffer_t> buffers{ buffer_t{} };

  public:
    explicit BackPropagatorLossy() = default;

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
        size += 3 * sizeof(std::size_t); // m_num_lanes, from_prev, to_prev
        size += sizeof(std::size_t) * (capacity ? node_location_on_buffer.capacity() : node_location_on_buffer.size());
        size += sizeof(std::uint8_t) * (capacity ? on_which_buffer.capacity() : on_which_buffer.size());
        size += sizeof(std::size_t) * (capacity ? pos.capacity() : pos.size());
        size += sizeof(LossyOpCode) * (capacity ? lossy_op.capacity() : lossy_op.size());
        size += sizeof(double) * (capacity ? values.capacity() : values.size());
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

template<class Float, bool Vectorised>
template<bool Reset, bool ResetInPlace, bool Log>
void
BackPropagatorLossy<Float, Vectorised>::backpropagate_to(std::size_t to, TapeData& data)
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
    if (lossy_op.empty() || from_prev != from || to_prev != to) {
        from_prev = from;
        to_prev = to;

        this->on_which_buffer.clear();
        this->pos.clear();
        this->lossy_op.clear();
        this->values.clear();

        this->node_location_on_buffer.resize(ops.size(), passive_id<std::size_t>);

        auto& buffer_free_positions = buffers.back().free_positions;

        // output is: is_new, is_current, buffer_id, position
        auto const& checkpoints_c = this->checkpoints;
        auto get_loc = [checkpoints_c](std::size_t id) -> std::tuple<bool, std::uint8_t> {
            auto it = std::upper_bound(checkpoints_c.begin(), checkpoints_c.end(), id);
            auto buffed_id = static_cast<std::uint8_t>(std::distance(checkpoints_c.cbegin(), it) - 1);
            return { it == checkpoints_c.end(), buffed_id };
        };

        auto update_loc = [&pos = this->pos, &on_which_buffer = this->on_which_buffer, &buffers = this->buffers](
                            std::size_t res_pos, std::size_t& arg_pos, std::uint8_t buffer_id) {
            pos.push_back(res_pos);

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

            pos.push_back(arg_pos);
            on_which_buffer.push_back(buffer_id);
        };

        auto copy_add =
          [update_loc, &lossy_op = this->lossy_op](std::size_t res_pos, std::size_t& arg_pos, std::uint8_t buffer_id) {
              bool arg_is_new = (arg_pos == passive_id<std::size_t>);
              update_loc(res_pos, arg_pos, buffer_id);

              if (arg_is_new) {
                  // this is a new value, we NEED to override
                  lossy_op.push_back(LossyOpCode::COPY);
              }
              else {
                  lossy_op.push_back(LossyOpCode::ADD);
              }
          };

        auto copy_sub =
          [update_loc, &lossy_op = this->lossy_op](std::size_t res_pos, std::size_t& arg_pos, std::uint8_t buffer_id) {
              bool arg_is_new = (arg_pos == passive_id<std::size_t>);
              update_loc(res_pos, arg_pos, buffer_id);

              if (arg_is_new) {
                  // this is a new value, we NEED to override
                  lossy_op.push_back(LossyOpCode::COPY_MINUS);
              }
              else {
                  lossy_op.push_back(LossyOpCode::SUB);
              }
          };

        auto copy_mul = [update_loc, &lossy_op = this->lossy_op, &values = this->values](
                          std::size_t res_pos, std::size_t& arg_pos, std::uint8_t buffer_id, double multiplier) {
            bool arg_is_new = (arg_pos == passive_id<std::size_t>);
            update_loc(res_pos, arg_pos, buffer_id);

            if (arg_is_new) {
                // this is a new value, we NEED to override
                lossy_op.push_back(LossyOpCode::MUL_SET);
            }
            else {
                lossy_op.push_back(LossyOpCode::MUL_ADD);
            }
            values.push_back(multiplier);
        };

        auto update_univariate =
          [this, copy_mul, get_loc, &node_location_on_buffer = this->node_location_on_buffer, &buffer_free_positions](
            std::size_t const arg_id, std::size_t const res_id, double const der_local_1) {
              auto const res_pos = node_location_on_buffer[res_id];
              node_location_on_buffer[res_id] = passive_id<std::size_t>;
              std::size_t& arg_pos = node_location_on_buffer[arg_id];

              auto const arg_pos_data = get_loc(arg_id);
              bool arg_is_new = (arg_pos == passive_id<std::size_t>);
              bool arg_inplace = arg_is_new && std::get<0>(arg_pos_data);

              if (arg_inplace) {
                  // res id should now be arg id, avoiding a copy and a potential buffer increase
                  arg_pos = res_pos;
                  this->pos.push_back(res_pos);
                  this->values.push_back(der_local_1);
                  this->lossy_op.push_back(LossyOpCode::MUL_INPLACE);
              }
              else {
                  copy_mul(res_pos, arg_pos, std::get<1>(arg_pos_data), der_local_1);
                  buffer_free_positions.push_back(res_pos);
              }
          };

        for (std::size_t op_idx = from; op_idx-- > to;) {
            OpCode const& op = ops[op_idx];
            bool const use_this_op = this->node_location_on_buffer[op_idx] != passive_id<std::size_t>;

            switch (op) {
                case OpCode::REG_INPUT: {
                    id_idx -= 1;
                    break;
                }
                case OpCode::REG_OUTPUT: {
                    id_idx -= 2;
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];

                        auto const res_pos = node_location_on_buffer[res_id];
                        std::size_t& arg_pos = node_location_on_buffer[arg_id];

                        auto const arg_pos_data = get_loc(arg_id);
                        copy_add(res_pos, arg_pos, std::get<1>(arg_pos_data));
                    }
                    break;
                }
                case OpCode::ADD: {
                    id_idx -= 3;
                    if (use_this_op) {
                        std::size_t const lhs_id = ids[id_idx];
                        std::size_t const rhs_id = ids[id_idx + 1];
                        std::size_t const res_id = ids[id_idx + 2];

                        auto const res_pos = node_location_on_buffer[res_id];
                        node_location_on_buffer[res_id] = passive_id<std::size_t>;
                        std::size_t& lhs_pos = node_location_on_buffer[lhs_id];
                        std::size_t& rhs_pos = node_location_on_buffer[rhs_id];

                        auto const lhs_pos_data = get_loc(lhs_id);
                        auto const rhs_pos_data = get_loc(rhs_id);

                        bool lhs_is_new = (lhs_pos == passive_id<std::size_t>);
                        bool rhs_is_new = (rhs_pos == passive_id<std::size_t>);
                        bool lhs_inplace = lhs_is_new && std::get<0>(lhs_pos_data);
                        bool rhs_inplace = !lhs_inplace && rhs_is_new && std::get<0>(rhs_pos_data);

                        if (!lhs_inplace) {
                            copy_add(res_pos, lhs_pos, std::get<1>(lhs_pos_data));
                        }

                        if (!rhs_inplace) {
                            copy_add(res_pos, rhs_pos, std::get<1>(rhs_pos_data));
                        }

                        if (lhs_inplace) {
                            // res id should now be lhs id, avoiding a copy and a potential buffer increase
                            lhs_pos = res_pos;
                        }
                        else if (rhs_inplace) {
                            // res id should now be rhs id, avoiding a copy and a potential buffer increase
                            rhs_pos = res_pos;
                        }
                        else {
                            // don't forget to free res_id from the buffer!
                            buffer_free_positions.push_back(res_pos);
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

                        auto const res_pos = node_location_on_buffer[res_id];
                        node_location_on_buffer[res_id] = passive_id<std::size_t>;
                        std::size_t& lhs_pos = node_location_on_buffer[lhs_id];
                        std::size_t& rhs_pos = node_location_on_buffer[rhs_id];

                        auto const lhs_pos_data = get_loc(lhs_id);
                        auto const rhs_pos_data = get_loc(rhs_id);

                        bool lhs_is_new = (lhs_pos == passive_id<std::size_t>);
                        bool rhs_is_new = (rhs_pos == passive_id<std::size_t>);
                        bool lhs_inplace = lhs_is_new && std::get<0>(lhs_pos_data);
                        bool rhs_inplace = !lhs_inplace && rhs_is_new && std::get<0>(rhs_pos_data);

                        if (!lhs_inplace) {
                            copy_add(res_pos, lhs_pos, std::get<1>(lhs_pos_data));
                        }

                        if (!rhs_inplace) {
                            copy_sub(res_pos, rhs_pos, std::get<1>(rhs_pos_data));
                        }

                        if (lhs_inplace) {
                            // res id should now be lhs id, avoiding a copy and a potential buffer increase
                            lhs_pos = res_pos;
                        }
                        else if (rhs_inplace) {
                            // res id should now be rhs id, avoiding a copy and a potential buffer increase
                            rhs_pos = res_pos;
                            // however this is a subtraction, so we need to negate the value in the buffer
                            // MINUS_INPLACE is always in buffer
                            this->lossy_op.push_back(LossyOpCode::MINUS_INPLACE);
                            this->pos.push_back(rhs_pos);
                        }
                        else {
                            // don't forget to free res_id from the buffer!
                            buffer_free_positions.push_back(res_pos);
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

                        auto const res_pos = node_location_on_buffer[res_id];
                        node_location_on_buffer[res_id] = passive_id<std::size_t>;
                        std::size_t& lhs_pos = node_location_on_buffer[lhs_id];
                        std::size_t& rhs_pos = node_location_on_buffer[rhs_id];

                        auto const lhs_pos_data = get_loc(lhs_id);
                        auto const rhs_pos_data = get_loc(rhs_id);

                        bool lhs_is_new = (lhs_pos == passive_id<std::size_t>);
                        bool rhs_is_new = (rhs_pos == passive_id<std::size_t>);
                        bool lhs_inplace = lhs_is_new && std::get<0>(lhs_pos_data);
                        bool rhs_inplace = !lhs_inplace && rhs_is_new && std::get<0>(rhs_pos_data);

                        if (!lhs_inplace) {
                            copy_mul(res_pos, lhs_pos, std::get<1>(lhs_pos_data), rhs_val);
                        }

                        if (!rhs_inplace) {
                            copy_mul(res_pos, rhs_pos, std::get<1>(rhs_pos_data), lhs_val);
                        }

                        if (lhs_inplace) {
                            // res id should now be lhs id, avoiding a copy and a potential buffer increase
                            lhs_pos = res_pos;
                            this->pos.push_back(lhs_pos);
                            this->values.push_back(rhs_val);
                            this->lossy_op.push_back(LossyOpCode::MUL_INPLACE);
                        }
                        else if (rhs_inplace) {
                            // res id should now be rhs id, avoiding a copy and a potential buffer increase
                            rhs_pos = res_pos;
                            this->pos.push_back(rhs_pos);
                            this->values.push_back(lhs_val);
                            this->lossy_op.push_back(LossyOpCode::MUL_INPLACE);
                        }
                        else {
                            // don't forget to free res_id from the buffer!
                            buffer_free_positions.push_back(res_pos);
                        }
                    }
                    break;
                }
                case OpCode::ADD_C: {
                    id_idx -= 2;
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];

                        auto const res_pos = node_location_on_buffer[res_id];
                        node_location_on_buffer[res_id] = passive_id<std::size_t>;
                        std::size_t& arg_pos = node_location_on_buffer[arg_id];

                        auto const arg_pos_data = get_loc(arg_id);
                        bool arg_is_new = (arg_pos == passive_id<std::size_t>);
                        bool arg_inplace = arg_is_new && std::get<0>(arg_pos_data);

                        if (arg_inplace) {
                            // res id should now be lhs id, avoiding a copy and a potential buffer increase
                            arg_pos = res_pos;
                        }
                        else {
                            copy_add(res_pos, arg_pos, std::get<1>(arg_pos_data));
                            buffer_free_positions.push_back(res_pos);
                        }
                    }
                    break;
                }
                case OpCode::SUB_C: {
                    id_idx -= 2;
                    if (use_this_op) {
                        std::size_t const arg_id = ids[id_idx];
                        std::size_t const res_id = ids[id_idx + 1];

                        auto const res_pos = node_location_on_buffer[res_id];
                        node_location_on_buffer[res_id] = passive_id<std::size_t>;
                        std::size_t& arg_pos = node_location_on_buffer[arg_id];

                        auto const arg_pos_data = get_loc(arg_id);
                        bool arg_is_new = (arg_pos == passive_id<std::size_t>);
                        bool arg_inplace = arg_is_new && std::get<0>(arg_pos_data);

                        if (arg_inplace) {
                            // res id should now be lhs id, avoiding a copy and a potential buffer increase
                            arg_pos = res_pos;
                            // however this is a subtraction, so we need to negate the value in the buffer
                            // MINUS_INPLACE is always in buffer
                            this->lossy_op.push_back(LossyOpCode::MINUS_INPLACE);
                            this->pos.push_back(arg_pos);
                        }
                        else {
                            copy_sub(res_pos, arg_pos, std::get<1>(arg_pos_data));
                            buffer_free_positions.push_back(res_pos);
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
                    val_idx -= 3;
                    id_idx -= 2;
                    if (use_this_op) {
                        double const one_over_in = 1. / vals[val_idx];
                        double const der_local_1 = vals[val_idx + 2] * vals[val_idx + 1] * one_over_in;
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
        data.ops.resize(to);
        data.vals.resize(val_idx);
        data.ids.resize(id_idx);
        data.next_id = to;

        this->node_location_on_buffer.resize(to);

        from_prev = 0;
        to_prev = 0;
    }

    for (auto& b : buffers) {
        b.values.resize(b.size * this->m_num_lanes, 0.);
    }
    auto& buffer = buffers.back().values;

    std::size_t val_idx_l = 0;
    std::size_t id_idx_l = 0;
    std::size_t on_which_buffer_idx = 0;

    for (std::size_t lossy_op_idx = 0; lossy_op_idx < this->lossy_op.size(); ++lossy_op_idx) {
        auto const& op = this->lossy_op[lossy_op_idx];

        switch (op) {
            case LossyOpCode::COPY: {
                std::uint8_t const which = this->on_which_buffer[on_which_buffer_idx++];
                std::size_t const out_pos = this->pos[id_idx_l++];
                std::size_t const in_pos = this->pos[id_idx_l++];
                if constexpr (Vectorised) {
                    const double* src = &buffer[out_pos * this->m_num_lanes];
                    double* dest = &buffers[which].values[in_pos * this->m_num_lanes];
#pragma omp simd
                    for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                        dest[i] = src[i];
                    }
                }
                else {
                    buffers[which].values[in_pos] = buffer[out_pos];
                }
                break;
            }
            case LossyOpCode::COPY_MINUS: {
                std::uint8_t const which = this->on_which_buffer[on_which_buffer_idx++];
                std::size_t const out_pos = this->pos[id_idx_l++];
                std::size_t const in_pos = this->pos[id_idx_l++];
                if constexpr (Vectorised) {
                    const double* src = &buffer[out_pos * this->m_num_lanes];
                    double* dest = &buffers[which].values[in_pos * this->m_num_lanes];
#pragma omp simd
                    for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                        dest[i] = -src[i];
                    }
                }
                else {
                    buffers[which].values[in_pos] = -buffer[out_pos];
                }
                break;
            }
            case LossyOpCode::ADD: {
                std::uint8_t const which = this->on_which_buffer[on_which_buffer_idx++];
                std::size_t const out_pos = this->pos[id_idx_l++];
                std::size_t const in_pos = this->pos[id_idx_l++];
                if constexpr (Vectorised) {
                    const double* src = &buffer[out_pos * this->m_num_lanes];
                    double* dest = &buffers[which].values[in_pos * this->m_num_lanes];
#pragma omp simd
                    for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                        dest[i] += src[i];
                    }
                }
                else {
                    buffers[which].values[in_pos] += buffer[out_pos];
                }
                break;
            }
            case LossyOpCode::SUB: {
                std::uint8_t const which = this->on_which_buffer[on_which_buffer_idx++];
                std::size_t const out_pos = this->pos[id_idx_l++];
                std::size_t const in_pos = this->pos[id_idx_l++];
                if constexpr (Vectorised) {
                    const double* src = &buffer[out_pos * this->m_num_lanes];
                    double* dest = &buffers[which].values[in_pos * this->m_num_lanes];
#pragma omp simd
                    for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                        dest[i] -= src[i];
                    }
                }
                else {
                    buffers[which].values[in_pos] -= buffer[out_pos];
                }
                break;
            }
            case LossyOpCode::MINUS_INPLACE: {
                std::size_t const pos = this->pos[id_idx_l++];
                if constexpr (Vectorised) {
                    double* dest = &buffer[pos * this->m_num_lanes];
#pragma omp simd
                    for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                        dest[i] = -dest[i];
                    }
                }
                else {
                    buffer[pos] = -buffer[pos];
                }
                break;
            }
            case LossyOpCode::MUL_INPLACE: {
                std::size_t const pos = this->pos[id_idx_l++];
                double const multiplier = this->values[val_idx_l++];
                if constexpr (Vectorised) {
                    double* dest = &buffer[pos * this->m_num_lanes];
#pragma omp simd
                    for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                        dest[i] *= multiplier;
                    }
                }
                else {
                    buffer[pos] *= multiplier;
                }
                break;
            }
            case LossyOpCode::MUL_ADD: {
                std::uint8_t const which = this->on_which_buffer[on_which_buffer_idx++];
                std::size_t const out_pos = this->pos[id_idx_l++];
                std::size_t const in_pos = this->pos[id_idx_l++];
                double const multiplier = this->values[val_idx_l++];
                if constexpr (Vectorised) {
                    const double* src = &buffer[out_pos * this->m_num_lanes];
                    double* dest = &buffers[which].values[in_pos * this->m_num_lanes];
#pragma omp simd
                    for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                        dest[i] += src[i] * multiplier;
                    }
                }
                else {
                    buffers[which].values[in_pos] += buffer[out_pos] * multiplier;
                }
                break;
            }
            case LossyOpCode::MUL_SET: {
                std::uint8_t const which = this->on_which_buffer[on_which_buffer_idx++];
                std::size_t const out_pos = this->pos[id_idx_l++];
                std::size_t const in_pos = this->pos[id_idx_l++];
                double const multiplier = this->values[val_idx_l++];
                if constexpr (Vectorised) {
                    const double* src = &buffer[out_pos * this->m_num_lanes];
                    double* dest = &buffers[which].values[in_pos * this->m_num_lanes];
#pragma omp simd
                    for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                        dest[i] = src[i] * multiplier;
                    }
                }
                else {
                    buffers[which].values[in_pos] = buffer[out_pos] * multiplier;
                }
                break;
            }
        }
    }

    if constexpr (Reset) {
        if (!checkpoints.empty() && to == checkpoints.back()) {
            // buffers.pop_back();
            buffers.back() = buffer_t{};
        }
    }
}

} // namespace adhoc

#endif // ADHOC_BACKPROPAGATOR1LOSSY_HPP
