#ifndef ADHOC_BACKPROPAGATOR1_HPP
#define ADHOC_BACKPROPAGATOR1_HPP

#include "tape_data.hpp"

#include <cmath>
#include <iostream>
#include <numbers>
#include <vector>

namespace adhoc {

template<class Float, bool Vectorised = false>
class BackPropagator {
  private:
    std::size_t m_num_lanes{ 1 };
    std::vector<double> derivatives;
    std::vector<char> use_op;

  public:
    explicit BackPropagator() = default;

    void set_checkpoint(std::size_t /* ops_size */) {}
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
    void reserve_input(std::size_t /* count_registered */) {}
    void reserve_output(std::size_t /* count_registered */) {}
    void register_variable(std::size_t /* var_id */) {}
    void register_output_variable(std::size_t var_id, std::size_t ops_size)
    {
        this->use_op.resize(ops_size);
        this->use_op[var_id] = true;
    }

    void set_derivative(std::size_t var_id, double deriv, std::size_t ops_size, std::size_t lane = 0)
    {
        if (lane < this->m_num_lanes) {
            this->derivatives.resize(ops_size * this->m_num_lanes);
            this->derivatives[var_id * this->m_num_lanes + lane] = deriv;
        }
        else {
            throw;
        }
    }

    void set_derivative(std::size_t /* var_id1 */, std::size_t /* var_id2 */, double /* deriv */) { throw; }

    auto get_derivative(std::size_t var_id, std::size_t lane) const -> double
    {
        if (lane < this->m_num_lanes) {
            return this->derivatives[var_id * this->m_num_lanes + lane];
        }

        throw;
        return 0.;
    }

    auto get_derivative(std::size_t /* var_id1 */, std::size_t /* var_id2 */, std::size_t /* lane */) const -> double
    {
        return 0.;
    }

    void clear()
    {
        this->derivatives.clear();
        this->use_op.clear();
    }

    void zero_adjoints() { std::fill(this->derivatives.begin(), this->derivatives.end(), 0.0); }

    auto size_of(bool capacity = false) const -> std::size_t
    {
        std::size_t size = 0;
        size += sizeof(std::size_t); // m_num_lanes
        size += sizeof(double) * (capacity ? this->derivatives.capacity() : this->derivatives.size());
        size += sizeof(char) * (capacity ? this->use_op.capacity() : this->use_op.size());
        return size;
    }

    template<bool Reset = false, bool ResetInPlace = false, bool Log = false>
    void backpropagate_to(std::size_t to, TapeData& data);
};

template<class Float, bool Vectorised>
template<bool Reset, bool ResetInPlace, bool Log>
void
BackPropagator<Float, Vectorised>::backpropagate_to(std::size_t to, TapeData& data)
{
    std::size_t from = data.next_id;
    const auto& ops = data.ops;
    const auto& vals = data.vals;
    const auto& ids = data.ids;

    this->derivatives.resize(ops.size() * this->m_num_lanes);
    this->use_op.resize(ops.size());

    auto add_derivative = [&](std::size_t arg_id, std::size_t res_id) {
        use_op[arg_id] = true;
        if constexpr (Vectorised) {
            const double* src = &this->derivatives[res_id * this->m_num_lanes];
            double* dest = &this->derivatives[arg_id * this->m_num_lanes];
#pragma omp simd
            for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                dest[i] += src[i];
            }
        }
        else {
            this->derivatives[arg_id] += this->derivatives[res_id];
        }
    };

    auto sub_derivative = [&](std::size_t arg_id, std::size_t res_id) {
        use_op[arg_id] = true;
        if constexpr (Vectorised) {
            const double* src = &this->derivatives[res_id * this->m_num_lanes];
            double* dest = &this->derivatives[arg_id * this->m_num_lanes];
#pragma omp simd
            for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                dest[i] -= src[i];
            }
        }
        else {
            this->derivatives[arg_id] -= this->derivatives[res_id];
        }
    };

    auto mul_add_derivative = [&](std::size_t arg_id, std::size_t res_id, double multiplier) {
        use_op[arg_id] = true;
        if constexpr (Vectorised) {
            const double* src = &this->derivatives[res_id * this->m_num_lanes];
            double* dest = &this->derivatives[arg_id * this->m_num_lanes];
#pragma omp simd
            for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                dest[i] += src[i] * multiplier;
            }
        }
        else {
            this->derivatives[arg_id] += this->derivatives[res_id] * multiplier;
        }
    };

    std::size_t val_idx = vals.size();
    std::size_t id_idx = ids.size();
    for (std::size_t op_idx = from; op_idx-- > to;) {
        OpCode const& op = ops[op_idx];
        bool const use_this_op = use_op[op_idx];

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
                    add_derivative(arg_id, res_id);
                }
                break;
            }
            case OpCode::ADD: {
                id_idx -= 3;
                if (use_this_op) {
                    std::size_t const lhs_id = ids[id_idx];
                    std::size_t const rhs_id = ids[id_idx + 1];
                    std::size_t const res_id = ids[id_idx + 2];
                    add_derivative(lhs_id, res_id);
                    add_derivative(rhs_id, res_id);
                }
                break;
            }
            case OpCode::SUB: {
                id_idx -= 3;
                if (use_this_op) {
                    std::size_t const lhs_id = ids[id_idx];
                    std::size_t const rhs_id = ids[id_idx + 1];
                    std::size_t const res_id = ids[id_idx + 2];
                    add_derivative(lhs_id, res_id);
                    sub_derivative(rhs_id, res_id);
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
                    mul_add_derivative(lhs_id, res_id, rhs_val);
                    mul_add_derivative(rhs_id, res_id, lhs_val);
                }
                break;
            }
            case OpCode::ADD_C: {
                id_idx -= 2;
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    add_derivative(arg_id, res_id);
                }
                break;
            }
            case OpCode::SUB_C: {
                id_idx -= 2;
                if (use_this_op) {
                    std::size_t const arg_id = ids[id_idx];
                    std::size_t const res_id = ids[id_idx + 1];
                    sub_derivative(arg_id, res_id);
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
                    mul_add_derivative(arg_id, res_id, der_local_1);
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
                    mul_add_derivative(arg_id, res_id, der_local_1);
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
                    mul_add_derivative(arg_id, res_id, der_local_1);
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
                    mul_add_derivative(arg_id, res_id, der_local_1);
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
                    mul_add_derivative(arg_id, res_id, der_local_1);
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
                    mul_add_derivative(arg_id, res_id, der_local_1);
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
                    mul_add_derivative(arg_id, res_id, der_local_1);
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
                    mul_add_derivative(arg_id, res_id, der_local_1);
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
                    mul_add_derivative(arg_id, res_id, der_local_1);
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
                    mul_add_derivative(arg_id, res_id, der_local_1);
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
                    mul_add_derivative(arg_id, res_id, der_local_1);
                }
                break;
            }
        }
        if constexpr (Reset && ResetInPlace) {
            this->derivatives.resize(this->derivatives.size() - this->m_num_lanes);
            this->use_op.resize(this->use_op.size() - 1);
        }
    }

    if constexpr (Log) {
        // calculate the percentage of operations saved using the vector use_op
        std::size_t avoided_ops = 0;
        for (std::size_t i = 0; i < use_op.size(); ++i) {
            avoided_ops += !use_op[i];
        }

        double percentage = static_cast<double>(avoided_ops) / static_cast<double>(use_op.size());
        std::cout.precision(std::numeric_limits<double>::max_digits10);
        std::cout << "backpropagation saved operations: " << percentage << std::endl;
    }

    if constexpr (Reset && !ResetInPlace) {
        this->derivatives.resize(to * this->m_num_lanes);
        this->use_op.resize(to);
    }

    if constexpr (Reset) {
        data.ops.resize(to);
        data.vals.resize(val_idx);
        data.ids.resize(id_idx);
        data.next_id = to;
    }
}

} // namespace adhoc

#endif // ADHOC_BACKPROPAGATOR1_HPP
