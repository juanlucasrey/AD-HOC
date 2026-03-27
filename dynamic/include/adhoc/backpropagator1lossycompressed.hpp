#ifndef ADHOC_BACKPROPAGATOR1LOSSYCOMPRESSED_HPP
#define ADHOC_BACKPROPAGATOR1LOSSYCOMPRESSED_HPP

#include "opcode_lossy.hpp"
#include "passive_id.hpp"
#include "tape_data.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

namespace adhoc {

// BackPropagatorLossy - creates and uses a LossyTape for backpropagation
template<class Float>
class BackPropagatorLossyCompressed {
  private:
    std::size_t m_num_lanes{ 1 };
    std::vector<std::size_t> input_ids;
    std::vector<std::size_t> output_ids;
    std::vector<double> input_derivs;
    std::vector<double> output_derivs;

    std::vector<bool> use_op;

    // lossy tape
    std::vector<bool> from_interface;
    std::vector<std::size_t> pos;
    std::vector<LossyOpCode> lossy_op;
    std::vector<double> values;

    std::size_t from_prev{ 0 };
    std::size_t to_prev{ 0 };
    std::size_t buffer_size_prev{ 0 };

    std::vector<std::size_t> buffer_active_id;
    std::vector<std::size_t> buffer_active_free_positions;
    std::vector<std::size_t> node_location_on_buffer_active;
    std::vector<double> buffer_active;

  public:
    explicit BackPropagatorLossyCompressed() = default;

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
    void reserve_input(std::size_t count_registered)
    {
        this->input_ids.reserve(this->input_ids.size() + count_registered);
        this->input_derivs.reserve(this->input_derivs.size() + (this->m_num_lanes * count_registered));
    }
    void reserve_output(std::size_t count_registered)
    {
        this->output_ids.reserve(this->output_ids.size() + count_registered);
        this->output_derivs.reserve(this->output_derivs.size() + (this->m_num_lanes * count_registered));
    }
    void register_variable(std::size_t var_id)
    {
        this->input_ids.push_back(var_id);
        this->input_derivs.resize(this->input_derivs.size() + this->m_num_lanes, 0.0);
    }
    void register_output_variable(std::size_t var_id, std::size_t /* ops_size */)
    {
        this->output_ids.push_back(var_id);
        this->output_derivs.resize(this->output_derivs.size() + this->m_num_lanes, 0.0);
    }

    void set_derivative(std::size_t var_id, double deriv, std::size_t /* ops_size */, std::size_t /* lane */ = 0)
    {
        const auto& reg = this->output_ids;
        auto it = std::lower_bound(reg.begin(), reg.end(), var_id);
        if (it == reg.cend() || *it != var_id) {
            // variable not registered!
            // NOT ALLOWED
            throw;
        }

        this->output_derivs[std::distance(reg.begin(), it)] = deriv;
    }

    void set_derivative(std::size_t /* var_id1 */, std::size_t /* var_id2 */, double /* deriv */)
    {
        // not supported for lossy tape (first order only)
        throw;
    }

    auto get_derivative(std::size_t var_id, std::size_t /* lane */) const -> double
    {
        const auto& reg = this->input_ids;
        auto it = std::lower_bound(reg.begin(), reg.end(), var_id);
        if (it == reg.cend() || *it != var_id) {
            // variable not registered!
            // NOT ALLOWED
            throw;
            return 0.;
        }

        return this->input_derivs[std::distance(reg.begin(), it)];
    }

    auto get_derivative(std::size_t /* var_id1 */, std::size_t /* var_id2 */, std::size_t /* lane */) const -> double
    {
        // order is too high for lossy tape
        return 0.;
    }

    void clear() {}

    void zero_adjoints()
    {
        std::fill(this->input_derivs.begin(), this->input_derivs.end(), 0.0);
        std::fill(this->output_derivs.begin(), this->output_derivs.end(), 0.0);
        clear();
    }

    auto size_of(bool /* capacity */ = false) const -> std::size_t { return 0; }

    template<bool Reset = false, bool ResetInPlace = false, bool Log = false>
    void backpropagate_to(std::size_t to, TapeData& data);
};

template<class Float>
template<bool Reset, bool ResetInPlace, bool Log>
void
BackPropagatorLossyCompressed<Float>::backpropagate_to(std::size_t to, TapeData& data)
{
    std::size_t from = data.next_id;
    if (this->m_num_lanes > 1) {
        // This backpropagator is not designed for multiple lanes
        throw;
    }

    const auto& reg_in = this->input_ids;
    const auto& reg_out = this->output_ids;
    auto& der_in = this->input_derivs;
    auto& der_out = this->output_derivs;
    std::size_t val_idx = 0;
    std::size_t id_idx = 0;
    std::vector<std::size_t> buffer_id;
    std::size_t input_decrease_counter = 0;
    std::size_t output_decrease_counter = 0;
    if (lossy_op.empty() || from_prev != from || to_prev != to) {
        from_prev = from;
        to_prev = to;
        from_interface.clear();
        pos.clear();
        lossy_op.clear();
        values.clear();
        const auto& ops = data.ops;
        const auto& vals = data.vals;
        const auto& ids = data.ids;

        // useful to avoid searching for most cases
        // registered inputs are in increasing order so the maximum is the last element
        const auto max_input_id = reg_in.back();

        this->use_op.resize(ops.size());

        // this boolean tells the location of the previous multiplication value
        std::vector<std::size_t> mulinplace_loc(from - to, passive_id<std::size_t>);
        std::vector<std::size_t> buffer_free_positions;
        std::vector<std::size_t> node_location_on_buffer(from, passive_id<std::size_t>);

        // loop over active buffer to insert buffer ids that are in between from and to:
        for (std::size_t i = 0; i < buffer_active_id.size(); ++i) {
            std::size_t const active_id = buffer_active_id[i];
            if (active_id != passive_id<std::size_t>) {
                if (active_id < from && active_id >= to) {
                    buffer_id.push_back(active_id);
                    node_location_on_buffer[active_id] = buffer_id.size() - 1;
                }
            }
        }

        // output is: in_input, in_buffer, position
        auto const& node_location_on_buffer_c = node_location_on_buffer;
        auto get_loc =
          [&reg_in, max_input_id, &node_location_on_buffer_c](std::size_t id) -> std::tuple<bool, bool, std::size_t> {
            if (id <= max_input_id) {
                auto it_in = std::lower_bound(reg_in.begin(), reg_in.end(), id);
                if (*it_in == id) {
                    return { true, false, std::distance(reg_in.cbegin(), it_in) };
                }
            }

            auto const position = node_location_on_buffer_c[id];
            return { false, position != passive_id<std::size_t>, position };
        };

        auto update_univariate =
          [this, &mulinplace_loc, to, &buffer_id, get_loc, &node_location_on_buffer, &buffer_free_positions](
            double const der_local_1, std::size_t const arg_id, std::size_t const res_id) {
              use_op[arg_id] = true;

              auto const position = node_location_on_buffer[res_id];
              if (position == passive_id<std::size_t>) {
                  throw;
              }

              auto const arg_pos_data = get_loc(arg_id);
              bool arg_inplace = !std::get<0>(arg_pos_data) && !std::get<1>(arg_pos_data);

              if (!arg_inplace) {
                  this->pos.push_back(position);
                  this->values.push_back(der_local_1);
                  this->from_interface.push_back(std::get<0>(arg_pos_data));
                  this->pos.push_back(std::get<2>(arg_pos_data));
                  this->lossy_op.push_back(LossyOpCode::MUL_ADD);
                  // adding invalidates a chain of multiplications
                  if (arg_id >= to) {
                      mulinplace_loc[arg_id - to] = passive_id<std::size_t>;
                  }
                  buffer_id[position] = passive_id<std::size_t>;
                  buffer_free_positions.push_back(position);
              }
              else {
                  std::size_t previous_multinplace_loc = mulinplace_loc[res_id - to];
                  if (previous_multinplace_loc != passive_id<std::size_t>) {
                      this->values[previous_multinplace_loc] *= der_local_1;
                      if (arg_id >= to) {
                          mulinplace_loc[arg_id - to] = previous_multinplace_loc;
                      }
                  }
                  else {
                      this->pos.push_back(position);
                      this->values.push_back(der_local_1);
                      this->lossy_op.push_back(LossyOpCode::MUL_INPLACE);
                      if (arg_id >= to) {
                          mulinplace_loc[arg_id - to] = this->values.size() - 1;
                      }
                  }
                  // res id should now be lhs id, avoiding a copy and a potential buffer increase
                  buffer_id[position] = arg_id;
                  node_location_on_buffer[arg_id] = position;
              }
          };

        val_idx = vals.size();
        id_idx = ids.size();

        for (std::size_t op_idx = from; op_idx-- > to;) {
            OpCode const& op = ops[op_idx];

            // OpCode::REG_OUTPUT is the only opcode that is always used in backpropagation
            if (op == OpCode::REG_OUTPUT) {
                id_idx -= 2;
                std::size_t const arg_id = ids[id_idx];
                std::size_t const res_id = ids[id_idx + 1];

                use_op[arg_id] = true;
                use_op[res_id] = true;
                std::size_t out_loc = reg_out.size() - 1 - output_decrease_counter;
                if (res_id != reg_out[out_loc]) {
                    throw;
                }

                this->pos.push_back(out_loc);
                output_decrease_counter++;

                auto const arg_pos_data = get_loc(arg_id);
                std::size_t res_position = std::get<2>(arg_pos_data);

                if (!std::get<0>(arg_pos_data) && !std::get<1>(arg_pos_data)) {
                    if (buffer_free_positions.empty()) {
                        buffer_id.push_back(arg_id);
                        res_position = buffer_id.size() - 1;
                    }
                    else {
                        res_position = buffer_free_positions.back();
                        buffer_free_positions.pop_back();
                        buffer_id[res_position] = arg_id;
                    }
                    node_location_on_buffer[arg_id] = res_position;
                }

                this->pos.push_back(res_position);
                this->from_interface.push_back(std::get<0>(arg_pos_data));
                this->lossy_op.push_back(LossyOpCode::ADD_FROM_SEED);

                if constexpr (Reset && ResetInPlace) {
                    data.ids.resize(data.ids.size() - 2);
                    use_op.resize(use_op.size() - 2);
                }
                continue;
            }

            if (!static_cast<bool>(use_op[op_idx])) {
                // skip operation
                switch (op) {
                    case OpCode::REG_INPUT:
                        input_decrease_counter++;
                        id_idx -= 1;
                        break;
                    case OpCode::REG_OUTPUT:
                        // should not reach here
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
                    case OpCode::INV:
                    case OpCode::EXP:
                    case OpCode::LOG:
                    case OpCode::ERF:
                    case OpCode::ERFC:
                        val_idx -= 1;
                        id_idx -= 2;
                        break;
                    case OpCode::SQRT:
                    case OpCode::COS:
                        val_idx -= 2;
                        id_idx -= 2;
                        break;
                    case OpCode::POW_C:
                        val_idx -= 3;
                        id_idx -= 2;
                        break;
                }

                continue;
            }

            switch (op) {
                case OpCode::REG_INPUT: {
                    input_decrease_counter++;
                    id_idx -= 1;
                    break;
                }
                case OpCode::REG_OUTPUT: {
                    // should not reach here
                    id_idx -= 2;
                    break;
                }
                case OpCode::ADD: {
                    id_idx -= 3;
                    std::size_t const lhs_id = ids[id_idx];
                    std::size_t const rhs_id = ids[id_idx + 1];
                    std::size_t const res_id = ids[id_idx + 2];
                    use_op[lhs_id] = true;
                    use_op[rhs_id] = true;

                    auto const res_buffer_pos = node_location_on_buffer[res_id];
                    if (res_buffer_pos == passive_id<std::size_t>) {
                        throw;
                    }

                    auto const lhs_pos_data = get_loc(lhs_id);
                    auto const rhs_pos_data = get_loc(rhs_id);

                    bool lhs_is_new = !std::get<0>(lhs_pos_data) && !std::get<1>(lhs_pos_data);
                    bool rhs_is_new = !std::get<0>(rhs_pos_data) && !std::get<1>(rhs_pos_data);
                    bool lhs_inplace = lhs_is_new;
                    bool rhs_inplace = !lhs_inplace && rhs_is_new;

                    if (!lhs_inplace) {
                        this->pos.push_back(res_buffer_pos);

                        this->from_interface.push_back(std::get<0>(lhs_pos_data));
                        this->pos.push_back(std::get<2>(lhs_pos_data));
                        this->lossy_op.push_back(LossyOpCode::ADD);
                        // adding invalidates a chain of multiplications
                        if (lhs_id >= to) {
                            mulinplace_loc[lhs_id - to] = passive_id<std::size_t>;
                        }
                    }

                    if (!rhs_inplace) {
                        this->pos.push_back(res_buffer_pos);

                        if (rhs_is_new) {
                            std::size_t free_pos = 0;
                            if (buffer_free_positions.empty()) {
                                buffer_id.push_back(rhs_id);
                                free_pos = buffer_id.size() - 1;
                            }
                            else {
                                free_pos = buffer_free_positions.back();
                                buffer_free_positions.pop_back();
                                buffer_id[free_pos] = rhs_id;
                            }

                            node_location_on_buffer[rhs_id] = free_pos;
                            this->pos.push_back(free_pos);
                            this->lossy_op.push_back(LossyOpCode::COPY);
                        }
                        else {
                            this->from_interface.push_back(std::get<0>(rhs_pos_data));
                            this->pos.push_back(std::get<2>(rhs_pos_data));
                            this->lossy_op.push_back(LossyOpCode::ADD);
                            // adding invalidates a chain of multiplications
                            if (rhs_id >= to) {
                                mulinplace_loc[rhs_id - to] = passive_id<std::size_t>;
                            }
                        }
                    }

                    if (lhs_inplace) {
                        // res id should now be lhs id, avoiding a copy and a potential buffer increase
                        buffer_id[res_buffer_pos] = lhs_id;
                        node_location_on_buffer[lhs_id] = res_buffer_pos;
                    }
                    else if (rhs_inplace) {
                        // res id should now be rhs id, avoiding a copy and a potential buffer increase
                        buffer_id[res_buffer_pos] = rhs_id;
                        node_location_on_buffer[rhs_id] = res_buffer_pos;
                    }
                    else {
                        // don't forget to free res_id from the buffer!
                        buffer_id[res_buffer_pos] = passive_id<std::size_t>;
                        buffer_free_positions.push_back(res_buffer_pos);
                    }
                    break;
                }
                case OpCode::SUB: {
                    id_idx -= 3;
                    std::size_t const lhs_id = ids[id_idx];
                    std::size_t const rhs_id = ids[id_idx + 1];
                    std::size_t const res_id = ids[id_idx + 2];
                    use_op[lhs_id] = true;
                    use_op[rhs_id] = true;

                    auto const res_buffer_pos = node_location_on_buffer[res_id];
                    if (res_buffer_pos == passive_id<std::size_t>) {
                        throw;
                    }

                    auto const lhs_pos_data = get_loc(lhs_id);
                    auto const rhs_pos_data = get_loc(rhs_id);

                    bool lhs_is_new = !std::get<0>(lhs_pos_data) && !std::get<1>(lhs_pos_data);
                    bool rhs_is_new = !std::get<0>(rhs_pos_data) && !std::get<1>(rhs_pos_data);
                    bool lhs_inplace = lhs_is_new;
                    bool rhs_inplace = !lhs_inplace && rhs_is_new;

                    if (!lhs_inplace) {
                        this->pos.push_back(res_buffer_pos);

                        this->from_interface.push_back(std::get<0>(lhs_pos_data));
                        this->pos.push_back(std::get<2>(lhs_pos_data));
                        this->lossy_op.push_back(LossyOpCode::ADD);
                        // adding invalidates a chain of multiplications
                        if (lhs_id >= to) {
                            mulinplace_loc[lhs_id - to] = passive_id<std::size_t>;
                        }
                    }

                    if (!rhs_inplace) {
                        this->pos.push_back(res_buffer_pos);

                        if (rhs_is_new) {
                            std::size_t free_pos = 0;
                            if (buffer_free_positions.empty()) {
                                buffer_id.push_back(rhs_id);
                                free_pos = buffer_id.size() - 1;
                            }
                            else {
                                free_pos = buffer_free_positions.back();
                                buffer_free_positions.pop_back();
                                buffer_id[free_pos] = rhs_id;
                            }

                            node_location_on_buffer[rhs_id] = free_pos;
                            this->pos.push_back(free_pos);
                            this->lossy_op.push_back(LossyOpCode::COPY_MINUS);
                        }
                        else {
                            this->from_interface.push_back(std::get<0>(rhs_pos_data));
                            this->pos.push_back(std::get<2>(rhs_pos_data));
                            this->lossy_op.push_back(LossyOpCode::SUB);
                            // adding invalidates a chain of multiplications
                            if (rhs_id >= to) {
                                mulinplace_loc[rhs_id - to] = passive_id<std::size_t>;
                            }
                        }
                    }

                    if (lhs_inplace) {
                        // res id should now be lhs id, avoiding a copy and a potential buffer increase
                        buffer_id[res_buffer_pos] = lhs_id;
                        node_location_on_buffer[lhs_id] = res_buffer_pos;
                    }
                    else if (rhs_inplace) {
                        this->pos.push_back(res_buffer_pos);
                        this->values.push_back(-1.0);
                        this->lossy_op.push_back(LossyOpCode::MUL_INPLACE);
                        if (rhs_id >= to) {
                            mulinplace_loc[rhs_id - to] = this->values.size() - 1;
                        }

                        // res id should now be rhs id, avoiding a copy and a potential buffer increase
                        buffer_id[res_buffer_pos] = rhs_id;
                        node_location_on_buffer[rhs_id] = res_buffer_pos;
                    }
                    else {
                        // don't forget to free res_id from the buffer!
                        buffer_id[res_buffer_pos] = passive_id<std::size_t>;
                        buffer_free_positions.push_back(res_buffer_pos);
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
                    use_op[lhs_id] = true;
                    use_op[rhs_id] = true;

                    auto const res_buffer_pos = node_location_on_buffer[res_id];
                    if (res_buffer_pos == passive_id<std::size_t>) {
                        throw;
                    }

                    auto const lhs_pos_data = get_loc(lhs_id);
                    auto const rhs_pos_data = get_loc(rhs_id);

                    bool lhs_is_new = !std::get<0>(lhs_pos_data) && !std::get<1>(lhs_pos_data);
                    bool rhs_is_new = !std::get<0>(rhs_pos_data) && !std::get<1>(rhs_pos_data);
                    bool lhs_inplace = lhs_is_new;
                    bool rhs_inplace = !lhs_inplace && rhs_is_new;

                    if (!lhs_inplace) {
                        this->pos.push_back(res_buffer_pos);

                        this->from_interface.push_back(std::get<0>(lhs_pos_data));
                        this->pos.push_back(std::get<2>(lhs_pos_data));
                        this->values.push_back(rhs_val);
                        this->lossy_op.push_back(LossyOpCode::MUL_ADD);
                        // adding invalidates a chain of multiplications
                        if (lhs_id >= to) {
                            mulinplace_loc[lhs_id - to] = passive_id<std::size_t>;
                        }
                    }

                    if (!rhs_inplace) {
                        this->pos.push_back(res_buffer_pos);

                        this->values.push_back(lhs_val);
                        if (rhs_is_new) {
                            std::size_t free_pos = 0;
                            if (buffer_free_positions.empty()) {
                                buffer_id.push_back(rhs_id);
                                free_pos = buffer_id.size() - 1;
                            }
                            else {
                                free_pos = buffer_free_positions.back();
                                buffer_free_positions.pop_back();
                                buffer_id[free_pos] = rhs_id;
                            }

                            node_location_on_buffer[rhs_id] = free_pos;
                            this->pos.push_back(free_pos);
                            this->lossy_op.push_back(LossyOpCode::MUL_SET);
                        }
                        else {
                            this->from_interface.push_back(std::get<0>(rhs_pos_data));
                            this->pos.push_back(std::get<2>(rhs_pos_data));
                            this->lossy_op.push_back(LossyOpCode::MUL_ADD);
                            // adding invalidates a chain of multiplications
                            if (rhs_id >= to) {
                                mulinplace_loc[rhs_id - to] = passive_id<std::size_t>;
                            }
                        }
                    }

                    if (lhs_inplace) {
                        // res id should now be lhs id, avoiding a copy and a potential buffer increase
                        buffer_id[res_buffer_pos] = lhs_id;
                        node_location_on_buffer[lhs_id] = res_buffer_pos;
                        this->pos.push_back(res_buffer_pos);
                        this->values.push_back(rhs_val);
                        this->lossy_op.push_back(LossyOpCode::MUL_INPLACE);
                        if (lhs_id >= to) {
                            mulinplace_loc[lhs_id - to] = this->values.size() - 1;
                        }
                    }
                    else if (rhs_inplace) {
                        // res id should now be rhs id, avoiding a copy and a potential buffer increase
                        buffer_id[res_buffer_pos] = rhs_id;
                        node_location_on_buffer[rhs_id] = res_buffer_pos;
                        this->pos.push_back(res_buffer_pos);
                        this->values.push_back(lhs_val);
                        this->lossy_op.push_back(LossyOpCode::MUL_INPLACE);
                        if (rhs_id >= to) {
                            mulinplace_loc[rhs_id - to] = this->values.size() - 1;
                        }
                    }
                    else {
                        // don't forget to free res_id from the buffer!
                        buffer_id[res_buffer_pos] = passive_id<std::size_t>;
                        buffer_free_positions.push_back(res_buffer_pos);
                    }
                    break;
                }
                case OpCode::ADD_C: {
                    id_idx -= 2;
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    use_op[arg_id] = true;

                    auto const res_buffer_pos = node_location_on_buffer[res_id];
                    if (res_buffer_pos == passive_id<std::size_t>) {
                        throw;
                    }

                    auto const arg_pos_data = get_loc(arg_id);
                    bool arg_inplace = !std::get<0>(arg_pos_data) && !std::get<1>(arg_pos_data);

                    if (!arg_inplace) {
                        this->pos.push_back(res_buffer_pos);

                        this->from_interface.push_back(std::get<0>(arg_pos_data));
                        this->pos.push_back(std::get<2>(arg_pos_data));
                        this->lossy_op.push_back(LossyOpCode::ADD);
                        buffer_id[res_buffer_pos] = passive_id<std::size_t>;
                        buffer_free_positions.push_back(res_buffer_pos);
                        // adding invalidates a chain of multiplications
                        if (arg_id >= to) {
                            mulinplace_loc[arg_id - to] = passive_id<std::size_t>;
                        }
                    }
                    else {
                        // res id should now be lhs id, avoiding a copy and a potential buffer increase
                        buffer_id[res_buffer_pos] = arg_id;
                        node_location_on_buffer[arg_id] = res_buffer_pos;

                        std::size_t previous_multinplace_loc = mulinplace_loc[res_id - to];
                        if (previous_multinplace_loc != passive_id<std::size_t>) {
                            if (arg_id >= to) {
                                mulinplace_loc[arg_id - to] = previous_multinplace_loc;
                            }
                        }
                    }
                    break;
                }
                case OpCode::SUB_C: {
                    id_idx -= 2;
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    use_op[arg_id] = true;

                    auto const res_buffer_pos = node_location_on_buffer[res_id];
                    if (res_buffer_pos == passive_id<std::size_t>) {
                        throw;
                    }

                    auto const arg_pos_data = get_loc(arg_id);

                    bool arg_inplace = !std::get<0>(arg_pos_data) && !std::get<1>(arg_pos_data);

                    if (!arg_inplace) {
                        this->pos.push_back(res_buffer_pos);
                        this->from_interface.push_back(std::get<0>(arg_pos_data));
                        this->pos.push_back(std::get<2>(arg_pos_data));
                        this->lossy_op.push_back(LossyOpCode::SUB);
                        buffer_id[res_buffer_pos] = passive_id<std::size_t>;
                        buffer_free_positions.push_back(res_buffer_pos);
                        // adding invalidates a chain of multiplications
                        if (arg_id >= to) {
                            mulinplace_loc[arg_id - to] = passive_id<std::size_t>;
                        }
                    }
                    else {
                        std::size_t previous_multinplace_loc = mulinplace_loc[res_id - to];
                        if (previous_multinplace_loc != passive_id<std::size_t>) {
                            this->values[previous_multinplace_loc] *= -1.0;
                            if (arg_id >= to) {
                                mulinplace_loc[arg_id - to] = previous_multinplace_loc;
                            }
                        }
                        else {
                            this->pos.push_back(res_buffer_pos);
                            this->values.push_back(-1.0);
                            this->lossy_op.push_back(LossyOpCode::MUL_INPLACE);
                            if (arg_id >= to) {
                                mulinplace_loc[arg_id - to] = this->values.size() - 1;
                            }
                        }

                        // res id should now be lhs id, avoiding a copy and a potential buffer increase
                        buffer_id[res_buffer_pos] = arg_id;
                        node_location_on_buffer[arg_id] = res_buffer_pos;
                    }
                    break;
                }
                case OpCode::MUL_C: {
                    val_idx -= 1;
                    id_idx -= 2;
                    double const der_local_1 = vals[val_idx];
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    update_univariate(der_local_1, arg_id, res_id);
                    break;
                }
                case OpCode::NORM: {
                    val_idx -= 1;
                    id_idx -= 2;
                    double const der_local_1 = 2.0 * vals[val_idx];
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    update_univariate(der_local_1, arg_id, res_id);
                    break;
                }
                case OpCode::INV: {
                    val_idx -= 1;
                    id_idx -= 2;
                    double const der_local_1 = -vals[val_idx] * vals[val_idx];
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    update_univariate(der_local_1, arg_id, res_id);
                    break;
                }
                case OpCode::ABS: {
                    val_idx -= 1;
                    id_idx -= 2;
                    // this could be improved but for now it will do
                    double const der_local_1 = std::copysign(1.0, vals[val_idx]);
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    update_univariate(der_local_1, arg_id, res_id);
                    break;
                }
                case OpCode::EXP: {
                    val_idx -= 1;
                    id_idx -= 2;
                    double const der_local_1 = vals[val_idx];
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    update_univariate(der_local_1, arg_id, res_id);
                    break;
                }
                case OpCode::LOG: {
                    val_idx -= 1;
                    id_idx -= 2;
                    double const der_local_1 = 1.0 / vals[val_idx];
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    update_univariate(der_local_1, arg_id, res_id);
                    break;
                }
                case OpCode::ERF: {
                    val_idx -= 1;
                    id_idx -= 2;
                    constexpr double minus_two_over_root_pi = 2. * std::numbers::inv_sqrtpi_v<double>;
                    double const der_local_1 = std::exp(-vals[val_idx] * vals[val_idx]) * minus_two_over_root_pi;
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    update_univariate(der_local_1, arg_id, res_id);
                    break;
                }
                case OpCode::ERFC: {
                    val_idx -= 1;
                    id_idx -= 2;
                    constexpr double minus_two_over_root_pi = -2. * std::numbers::inv_sqrtpi_v<double>;
                    double const der_local_1 = std::exp(-vals[val_idx] * vals[val_idx]) * minus_two_over_root_pi;
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    update_univariate(der_local_1, arg_id, res_id);
                    break;
                }
                case OpCode::COS: {
                    val_idx -= 2;
                    id_idx -= 2;
                    double const der_local_1 = -std::sin(vals[val_idx]);
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    update_univariate(der_local_1, arg_id, res_id);
                    break;
                }
                case OpCode::SQRT: {
                    val_idx -= 2;
                    id_idx -= 2;
                    double const one_over_in = 1. / vals[val_idx];
                    double const der_local_1 = 0.5 * vals[val_idx + 1] * one_over_in;
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    update_univariate(der_local_1, arg_id, res_id);
                    break;
                }
                case OpCode::POW_C: {
                    val_idx -= 3;
                    id_idx -= 2;
                    double const one_over_in = 1. / vals[val_idx];
                    double const rhs_arg = vals[val_idx + 2];
                    double const der_local_1 = rhs_arg * vals[val_idx + 1] * one_over_in;
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    update_univariate(der_local_1, arg_id, res_id);
                    break;
                }
            }
        }

        buffer_size_prev = buffer_id.size();
    }

    if constexpr (Reset) {
        data.ops.resize(to);
        data.vals.resize(val_idx);
        data.ids.resize(id_idx);
        data.next_id = to;

        use_op.resize(to);

        from_prev = 0;
        to_prev = 0;
        node_location_on_buffer_active.resize(to, passive_id<std::size_t>);
    }
    else {
        node_location_on_buffer_active.resize(from, passive_id<std::size_t>);
    }

    std::vector<double> buffer(buffer_size_prev, 0.0);
    // loop to insert buffer values that are in between from and to:
    std::size_t buffer_val_idx = 0;
    for (std::size_t i = 0; i < buffer_active_id.size(); ++i) {
        std::size_t const active_id = buffer_active_id[i];
        if (active_id != passive_id<std::size_t>) {
            if (active_id < from && active_id >= to) {
                buffer[buffer_val_idx] = buffer_active[i];
                ++buffer_val_idx;
                if constexpr (Reset) {
                    node_location_on_buffer_active[active_id] = passive_id<std::size_t>;
                    buffer_active_id[i] = passive_id<std::size_t>;
                    buffer_active_free_positions.push_back(i);
                }
            }
        }
    }

    val_idx = 0;
    id_idx = 0;
    std::size_t interface_idx = 0;

    for (std::size_t lossy_op_idx = 0; lossy_op_idx < this->lossy_op.size(); ++lossy_op_idx) {
        auto const& op = this->lossy_op[lossy_op_idx];

        switch (op) {
            case LossyOpCode::ADD_FROM_SEED: {
                bool const in_interface = this->from_interface[interface_idx++];
                std::size_t const out_pos = this->pos[id_idx++];
                std::size_t const in_pos = this->pos[id_idx++];

                if (in_interface) {
                    // maybe this should throw. why copy form output to intput directly?
                    der_in[in_pos] += der_out[out_pos];
                }
                else {
                    buffer[in_pos] += der_out[out_pos];
                }
                break;
            }
            case LossyOpCode::COPY: {
                std::size_t const out_pos = this->pos[id_idx++];
                std::size_t const in_pos = this->pos[id_idx++];
                // should never happen on the interface!!
                buffer[in_pos] = buffer[out_pos];
                break;
            }
            case LossyOpCode::COPY_MINUS: {
                std::size_t const out_pos = this->pos[id_idx++];
                std::size_t const in_pos = this->pos[id_idx++];
                // should never happen on the interface!!
                buffer[in_pos] = -buffer[out_pos];
                break;
            }
            case LossyOpCode::ADD: {
                bool const in_interface = this->from_interface[interface_idx++];
                std::size_t const out_pos = this->pos[id_idx++];
                std::size_t const in_pos = this->pos[id_idx++];

                // to separate with different operators
                if (in_interface) {
                    der_in[in_pos] += buffer[out_pos];
                }
                else {
                    buffer[in_pos] += buffer[out_pos];
                }
                break;
            }
            case LossyOpCode::SUB: {
                bool const in_interface = this->from_interface[interface_idx++];
                std::size_t const out_pos = this->pos[id_idx++];
                std::size_t const in_pos = this->pos[id_idx++];

                // to separate with different operators
                if (in_interface) {
                    der_in[in_pos] -= buffer[out_pos];
                }
                else {
                    buffer[in_pos] -= buffer[out_pos];
                }
                break;
            }
            case LossyOpCode::MINUS_INPLACE: {
                // std::size_t const pos = this->pos[id_idx++];
                // buffer[pos] = -buffer[pos];
                // disabled
                throw;
                break;
            }
            case LossyOpCode::MUL_INPLACE: {
                std::size_t const pos = this->pos[id_idx++];
                double const multiplier = this->values[val_idx++];
                buffer[pos] *= multiplier;
                break;
            }
            case LossyOpCode::MUL_ADD: {
                bool const in_interface = this->from_interface[interface_idx++];
                std::size_t const out_pos = this->pos[id_idx++];
                std::size_t const in_pos = this->pos[id_idx++];
                double const multiplier = this->values[val_idx++];

                // to separate with different operators
                if (in_interface) {
                    der_in[in_pos] += multiplier * buffer[out_pos];
                }
                else {
                    buffer[in_pos] += multiplier * buffer[out_pos];
                }
                break;
            }
            case LossyOpCode::MUL_SET: {
                std::size_t const out_pos = this->pos[id_idx++];
                std::size_t const in_pos = this->pos[id_idx++];
                double const multiplier = this->values[val_idx++];
                // should never happen on the interface!!
                buffer[in_pos] = multiplier * buffer[out_pos];
                break;
            }
        }
    }

    // Store non-passive buffer entries into global buffer_active
    for (std::size_t idx = 0; idx < buffer_id.size(); ++idx) {
        std::size_t current_id = buffer_id[idx];
        if (current_id != passive_id<std::size_t>) {
            std::size_t& loc = node_location_on_buffer_active[current_id];
            if (loc != passive_id<std::size_t>) {
                buffer_active[loc] += buffer[idx];
            }
            else {
                std::size_t free_pos = 0;
                if (buffer_active_free_positions.empty()) {
                    buffer_active_id.push_back(current_id);
                    buffer_active.resize(buffer_active.size() + 1, 0.0);
                    free_pos = buffer_active_id.size() - 1;
                }
                else {
                    free_pos = buffer_active_free_positions.back();
                    buffer_active_free_positions.pop_back();
                    buffer_active_id[free_pos] = current_id;
                }
                buffer_active[free_pos] = buffer[idx];
                loc = free_pos;
            }
        }
    }

    if constexpr (Reset) {
        der_in.resize(der_in.size() - input_decrease_counter);
        this->input_ids.resize(this->input_ids.size() - input_decrease_counter);

        der_out.resize(der_out.size() - output_decrease_counter);
        this->output_ids.resize(this->output_ids.size() - output_decrease_counter);

        buffer_size_prev = 0;
    }
}

} // namespace adhoc

#endif // ADHOC_BACKPROPAGATOR1LOSSYCOMPRESSED_HPP
