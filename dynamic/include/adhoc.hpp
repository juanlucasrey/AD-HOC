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

#ifndef ADHOC_HPP
#define ADHOC_HPP

#include <cmath>
#include <complex>
#include <cstddef>
#include <memory>
#include <vector>

#include "adhoc/passive_id.hpp"
#include "adhoc/tape_data.hpp"

namespace adhoc {

template<class Float>
class Tape;

template<class Float>
class adhoc_type;

template<typename T>
class opcode {
  public:
    using type = adhoc_type<T>;
    using tape_t = Tape<T>;
    inline static thread_local tape_t* global_tape = nullptr;
};

template<class Float>
class adhoc_type {
  private:
    friend Tape<Float>;

    double value{ 0. };
    mutable std::size_t id{ passive_id<std::size_t> };

    auto sub_c(double lhs) const -> adhoc_type;
    auto div_c(double lhs) const -> adhoc_type;

  public:
    using mode_t = opcode<Float>;

    adhoc_type() = default;

    adhoc_type(Float val)
      : value(val)
    {
    }
    adhoc_type(const adhoc_type& other);
    adhoc_type(const adhoc_type&& other);

    // Get value
    auto get_value() const -> double { return value; }
    auto get_value() -> double& { return value; }
    auto is_passive() const -> bool { return id == passive_id<std::size_t>; }
    auto is_active() const -> bool { return id != passive_id<std::size_t>; }

    // Declare operators - implementations come after Tape definition
    auto operator-() const -> adhoc_type;

    auto operator+(const adhoc_type& other) const -> adhoc_type;
    auto operator-(const adhoc_type& other) const -> adhoc_type;
    auto operator*(const adhoc_type& other) const -> adhoc_type;
    auto operator/(const adhoc_type& other) const -> adhoc_type;

    auto operator+(double other) const -> adhoc_type;
    auto operator-(double other) const -> adhoc_type;
    auto operator*(double other) const -> adhoc_type;
    auto operator/(double other) const -> adhoc_type;

    auto operator=(const adhoc_type& other) -> adhoc_type&;

    auto operator+=(const adhoc_type& other) -> adhoc_type&;
    auto operator-=(const adhoc_type& other) -> adhoc_type&;
    auto operator*=(const adhoc_type& other) -> adhoc_type&;
    auto operator/=(const adhoc_type& other) -> adhoc_type&;

    auto operator+=(double other) -> adhoc_type&;
    auto operator-=(double other) -> adhoc_type&;
    auto operator*=(double other) -> adhoc_type&;
    auto operator/=(double other) -> adhoc_type&;

    template<class Float2>
    friend auto exp(const adhoc_type<Float2>& x) -> adhoc_type<Float2>;

    template<class Float2>
    friend auto expm1(const adhoc_type<Float2>& x) -> adhoc_type<Float2>;

    template<class Float2>
    friend auto log(const adhoc_type<Float2>& x) -> adhoc_type<Float2>;

    template<class Float2>
    friend auto erf(const adhoc_type<Float2>& x) -> adhoc_type<Float2>;

    template<class Float2>
    friend auto erfc(const adhoc_type<Float2>& x) -> adhoc_type<Float2>;

    template<class Float2>
    friend auto cos(const adhoc_type<Float2>& x) -> adhoc_type<Float2>;

    template<class Float2>
    friend auto norm(const adhoc_type<Float2>& x) -> adhoc_type<Float2>;

    template<class Float2>
    friend auto inv(const adhoc_type<Float2>& x) -> adhoc_type<Float2>;

    template<class Float2>
    friend auto abs(const adhoc_type<Float2>& x) -> adhoc_type<Float2>;

    template<class Float2>
    friend auto sqrt(const adhoc_type<Float2>& x) -> adhoc_type<Float2>;

    // template<class Float2>
    // friend auto pow(const adhoc_type<Float2>& lhs, const adhoc_type<Float2>& rhs) -> adhoc_type<Float2>;

    // template<class T, class Float2>
    // inline auto pow(T /* lhs */, const adhoc_type<Float2>& /* rhs */) -> adhoc_type<Float2>;

    template<class Float2, class T>
    friend auto pow(const adhoc_type<Float2>& /* lhs */, T /* rhs */) -> adhoc_type<Float2>;

    template<class Float2>
    friend auto operator-(double lhs, const adhoc_type<Float2>& rhs) -> adhoc_type<Float2>;

    template<class Float2>
    friend auto operator/(double lhs, const adhoc_type<Float2>& rhs) -> adhoc_type<Float2>;
};

// Forward declare functions - implementations come after Tape definition
template<class Float>
auto
exp(const adhoc_type<Float>& x) -> adhoc_type<Float>;

template<class Float>
auto
expm1(const adhoc_type<Float>& x) -> adhoc_type<Float>;

template<class Float>
auto
log(const adhoc_type<Float>& x) -> adhoc_type<Float>;

template<class Float>
auto
erf(const adhoc_type<Float>& x) -> adhoc_type<Float>;

template<class Float>
auto
erfc(const adhoc_type<Float>& x) -> adhoc_type<Float>;

template<class Float>
auto
cos(const adhoc_type<Float>& x) -> adhoc_type<Float>;

template<class Float>
auto
norm(const adhoc_type<Float>& x) -> adhoc_type<Float>;

template<class Float>
auto
inv(const adhoc_type<Float>& x) -> adhoc_type<Float>;

template<class Float>
auto
sqrt(const adhoc_type<Float>& x) -> adhoc_type<Float>;

template<class Float, class T>
auto
pow(const adhoc_type<Float>& lhs, T rhs) -> adhoc_type<Float>;

enum class Method {
    FirstOrderSimple,
    FirstOrderSimd8,
    SecondOrderSimple,
    SecondOrderSimd8_stdmap,
    SecondOrderSimd8_stdunorderedmap,
    SecondOrderSimd8_ankerl,
    SecondOrderSimd8_boost,
    FirstOrderLossy,
    FirstOrderVLossy,
    FirstOrderLossyCompressed,
    FirstOrderVLossyCompressed
};

template<class Float>
class Tape {
  private:
    struct Impl;

    std::unique_ptr<Impl> impl;
    TapeData data;

    friend adhoc_type<Float>;

    template<class Float2>
    friend auto exp(const adhoc_type<Float2>& x) -> adhoc_type<Float2>;

    template<class Float2>
    friend auto expm1(const adhoc_type<Float2>& x) -> adhoc_type<Float2>;

    template<class Float2>
    friend auto log(const adhoc_type<Float2>& x) -> adhoc_type<Float2>;

    template<class Float2>
    friend auto erf(const adhoc_type<Float2>& x) -> adhoc_type<Float2>;

    template<class Float2>
    friend auto erfc(const adhoc_type<Float2>& x) -> adhoc_type<Float2>;

    template<class Float2>
    friend auto cos(const adhoc_type<Float2>& x) -> adhoc_type<Float2>;

    template<class Float2>
    friend auto norm(const adhoc_type<Float2>& x) -> adhoc_type<Float2>;

    template<class Float2>
    friend auto inv(const adhoc_type<Float2>& x) -> adhoc_type<Float2>;

    template<class Float2>
    friend auto sqrt(const adhoc_type<Float2>& x) -> adhoc_type<Float2>;

    template<class Float2, class T>
    friend auto pow(const adhoc_type<Float2>& /* lhs */, T /* rhs */) -> adhoc_type<Float2>;

    template<class Float2>
    friend auto abs(const adhoc_type<Float2>& x) -> adhoc_type<Float2>;

    void record_unary(OpCode op, std::size_t arg_id, std::size_t res_id)
    {
        data.ops.push_back(op);
        data.ids.push_back(arg_id);
        data.ids.push_back(res_id);
    }

    void record_value(double const value) { data.vals.push_back(value); }

    void record_binary(OpCode op, std::size_t lhs_id, std::size_t rhs_id, std::size_t res_id)
    {
        data.ops.push_back(op);
        data.ids.push_back(lhs_id);
        data.ids.push_back(rhs_id);
        data.ids.push_back(res_id);
    }

    void reserve_input(std::size_t n_registered);
    void reserve_output(std::size_t n_registered);

    auto generate_id() -> std::size_t { return data.next_id++; }

  public:
    void register_variable(adhoc_type<Float> const& var);
    void register_variable(adhoc_type<Float>& var);

    template<typename Container>
    void register_variable(Container&& vars)
    {
        std::size_t const count = std::distance(std::begin(vars), std::end(vars));
        this->reserve_input(count);
        for (auto&& var : vars) {
            register_variable(var);
        }
    }

    void register_output_variable(adhoc_type<Float> const& var);
    void register_output_variable(adhoc_type<Float>& var);

    template<typename Container>
    void register_output_variable(Container&& vars)
    {
        std::size_t const count = std::distance(std::begin(vars), std::end(vars));
        this->reserve_output(count);
        for (auto&& var : vars) {
            register_output_variable(var);
        }
    }

    // Implicit Function Theorem terminology: given f(x, y(x))
    //  -> x is a regular input (to be registered previously as an input)
    //  -> y(x) is an implicit variable (to be registered after the start of the
    //  implicit function)
    //  -> f() is a constraint variable (to be registered before the end of the
    //  implicit function)
    // void start_implicit_function() { ops.push_back(OpCode::IFT_START); }

    // void register_implicit_variable(adhoc_type<Float> &value) {
    //     value.tape = this;
    //     if (value.is_passive()) {
    //         value.id = generate_id();
    //         record_register(OpCode::REG_INPUT_IFT, value);
    //     }
    // }

    // void register_constraint_variable(adhoc_type<Float> &var) {
    //     if (var.is_active()) {
    //         record_register(OpCode::REG_OUTPUT_IFT, var);
    //     }
    // }

    // void end_implicit_function() { ops.push_back(OpCode::IFT_END); }

    Tape();
    ~Tape();

    // only for lossy tapes for now
    void set_lanes(std::size_t num_lanes = 1);
    auto get_lanes() const -> std::size_t;
    void set_method(Method m);
    auto get_method() const -> Method;
    auto get_order() const -> std::size_t;

    void clear()
    {
        data.ops.clear();
        data.vals.clear();
        data.ids.clear();
    }

    void set_checkpoint();
    void backpropagate();
    void backpropagate_to(std::size_t to);

    template<bool ResetInPlace = false, bool Log = false>
    void backpropagate_and_reset_to(std::size_t to);
    void set_derivative(adhoc_type<Float> const& var, double deriv, std::size_t lane = 0);
    auto get_derivative(adhoc_type<Float> const& var, std::size_t lane = 0) const -> double;

    auto get_derivative(adhoc_type<Float> const& var1, adhoc_type<Float> const& var2, std::size_t lane = 0) const
      -> double;
    void zero_adjoints();

    using position_t = std::size_t;

    auto get_position() const -> position_t { return data.next_id; }

    void print() const;

    auto size_of(bool capacity = false) const -> std::size_t;
};

template<class Float>
inline adhoc_type<Float>::adhoc_type(const adhoc_type& other)
  : value(other.value)
  , id(other.id)
{
}

template<class Float>
inline adhoc_type<Float>::adhoc_type(const adhoc_type&& other)
{
    *this = other;
}

template<class Float>
inline auto
adhoc_type<Float>::operator-() const -> adhoc_type
{
    adhoc_type result(-this->value);
    if (this->is_active()) {
        result.id = opcode<Float>::global_tape->generate_id();
        opcode<Float>::global_tape->record_unary(OpCode::SUB_C, this->id, result.id);
    }
    return result;
}

template<class Float>
inline auto
adhoc_type<Float>::operator+(const adhoc_type& other) const -> adhoc_type
{
    if (this->is_passive() && other.is_active()) {
        return other + this->value;
    }

    if (other.is_passive()) {
        return (*this) + other.value;
    }

    adhoc_type result(this->value + other.value);
    result.id = opcode<Float>::global_tape->generate_id();
    opcode<Float>::global_tape->record_binary(OpCode::ADD, this->id, other.id, result.id);
    return result;
}

template<class Float>
inline auto
adhoc_type<Float>::operator+(double other) const -> adhoc_type
{
    adhoc_type result(this->value + other);
    if (this->is_active()) {
        result.id = opcode<Float>::global_tape->generate_id();
        opcode<Float>::global_tape->record_unary(OpCode::ADD_C, this->id, result.id);
    }
    return result;
}

template<class Float>
inline auto
adhoc_type<Float>::operator-(const adhoc_type& other) const -> adhoc_type
{
    if (this->is_passive() && other.is_active()) {
        return this->value - other;
    }

    if (other.is_passive()) {
        return (*this) - other.value;
    }

    adhoc_type result(this->value - other.value);
    result.id = opcode<Float>::global_tape->generate_id();
    opcode<Float>::global_tape->record_binary(OpCode::SUB, this->id, other.id, result.id);
    return result;
}

template<class Float>
inline auto
adhoc_type<Float>::operator-(double other) const -> adhoc_type
{
    adhoc_type result(this->value - other);
    if (this->is_active()) {
        result.id = opcode<Float>::global_tape->generate_id();
        opcode<Float>::global_tape->record_unary(OpCode::ADD_C, this->id, result.id);
    }
    return result;
}

template<class Float>
inline auto
adhoc_type<Float>::operator*(const adhoc_type& other) const -> adhoc_type
{
    if (this->is_passive() && other.is_active()) {
        return other * this->value;
    }

    if (other.is_passive()) {
        return (*this) * other.value;
    }

    adhoc_type result(this->value * other.value);
    result.id = opcode<Float>::global_tape->generate_id();
    if (this->id == other.id) {
        opcode<Float>::global_tape->record_unary(OpCode::NORM, this->id, result.id);
        opcode<Float>::global_tape->record_value(this->value);
    }
    else {
        opcode<Float>::global_tape->record_binary(OpCode::MUL, this->id, other.id, result.id);
        opcode<Float>::global_tape->record_value(this->value);
        opcode<Float>::global_tape->record_value(other.value);
    }

    return result;
}

template<class Float>
inline auto
adhoc_type<Float>::operator*(double other) const -> adhoc_type
{
    adhoc_type result(this->value * other);
    if (this->is_active()) {
        result.id = opcode<Float>::global_tape->generate_id();
        opcode<Float>::global_tape->record_unary(OpCode::MUL_C, this->id, result.id);
        opcode<Float>::global_tape->record_value(other);
    }
    return result;
}

template<class Float>
inline auto
adhoc_type<Float>::operator/(const adhoc_type& other) const -> adhoc_type
{

    if (this->is_passive()) {
        return other.div_c(this->value);
    }

    if (other.is_passive()) {
        return *this / other.value;
    }

    adhoc_type result(this->value / other.value);

    adhoc_type intermediary_result(1.0 / other.value);
    intermediary_result.id = opcode<Float>::global_tape->generate_id();
    result.id = opcode<Float>::global_tape->generate_id();
    opcode<Float>::global_tape->record_unary(OpCode::INV, other.id, intermediary_result.id);
    opcode<Float>::global_tape->record_value(intermediary_result.value);
    opcode<Float>::global_tape->record_binary(OpCode::MUL, this->id, intermediary_result.id, result.id);
    opcode<Float>::global_tape->record_value(this->value);
    opcode<Float>::global_tape->record_value(intermediary_result.value);
    return result;
}

template<class Float>
inline auto
adhoc_type<Float>::operator/(double other) const -> adhoc_type
{
    adhoc_type result(this->value / other);
    if (this->is_active()) {
        result.id = opcode<Float>::global_tape->generate_id();
        opcode<Float>::global_tape->record_unary(OpCode::MUL_C, this->id, result.id);
        opcode<Float>::global_tape->record_value(1.0 / other);
    }
    return result;
}

template<class Float>
inline auto
adhoc_type<Float>::operator=(const adhoc_type& other) -> adhoc_type& = default;

template<class Float>
inline auto
adhoc_type<Float>::operator+=(const adhoc_type& other) -> adhoc_type&
{
    if (this->is_passive() && other.is_active()) {
        return *this = other + this->value;
    }

    if (other.is_passive()) {
        return *this += other.value;
    }

    return *this = *this + other;
}

template<class Float>
inline auto
adhoc_type<Float>::operator+=(double other) -> adhoc_type&
{
    return *this = *this + other;
}

template<class Float>
inline auto
adhoc_type<Float>::operator-=(const adhoc_type& other) -> adhoc_type&
{
    if (this->is_passive() && other.is_active()) {
        return *this = this->value - other;
    }

    if (other.is_passive()) {
        return (*this) -= other.value;
    }

    return *this = *this - other;
}

template<class Float>
inline auto
adhoc_type<Float>::operator-=(double other) -> adhoc_type&
{
    return *this = *this - other;
}

template<class Float>
inline auto
adhoc_type<Float>::operator*=(const adhoc_type& other) -> adhoc_type&
{
    if (this->is_passive() && other.is_active()) {
        return *this = other * this->value;
    }

    if (other.is_passive()) {
        return (*this) *= other.value;
    }

    return *this = *this * other;
}

template<class Float>
inline auto
adhoc_type<Float>::operator*=(double other) -> adhoc_type&
{
    return *this = *this * other;
}

template<class Float>
inline auto
adhoc_type<Float>::operator/=(const adhoc_type& other) -> adhoc_type&
{
    return *this = *this / other;
}

template<class Float>
inline auto
adhoc_type<Float>::operator/=(double other) -> adhoc_type&
{
    return *this = *this / other;
}

template<class Float>
inline auto
exp(const adhoc_type<Float>& arg) -> adhoc_type<Float>
{
    adhoc_type result(std::exp(arg.value));
    if (arg.is_active()) {
        result.id = opcode<Float>::global_tape->generate_id();
        opcode<Float>::global_tape->record_unary(OpCode::EXP, arg.id, result.id);
        opcode<Float>::global_tape->record_value(result.value);
    }
    return result;
}

template<class Float>
inline auto
expm1(const adhoc_type<Float>& arg) -> adhoc_type<Float>
{
    adhoc_type result(std::expm1(arg.value));

    if (arg.is_active()) {
        adhoc_type intermediary_result(std::exp(arg.value));
        intermediary_result.id = opcode<Float>::global_tape->generate_id();
        result.id = opcode<Float>::global_tape->generate_id();
        opcode<Float>::global_tape->record_unary(OpCode::EXP, arg.id, intermediary_result.id);
        opcode<Float>::global_tape->record_value(intermediary_result.value);
        opcode<Float>::global_tape->record_unary(OpCode::ADD_C, intermediary_result.id, result.id);
    }

    return result;
}

template<class Float>
inline auto
log(const adhoc_type<Float>& arg) -> adhoc_type<Float>
{
    adhoc_type result(std::log(arg.value));
    if (arg.is_active()) {
        result.id = opcode<Float>::global_tape->generate_id();
        opcode<Float>::global_tape->record_unary(OpCode::LOG, arg.id, result.id);
        opcode<Float>::global_tape->record_value(arg.value);
    }
    return result;
}

template<class Float>
inline auto
erf(const adhoc_type<Float>& arg) -> adhoc_type<Float>
{
    adhoc_type result(std::erf(arg.value));
    if (arg.is_active()) {
        result.id = opcode<Float>::global_tape->generate_id();
        opcode<Float>::global_tape->record_unary(OpCode::ERF, arg.id, result.id);
        opcode<Float>::global_tape->record_value(arg.value);
    }
    return result;
}

template<class Float>
inline auto
erfc(const adhoc_type<Float>& arg) -> adhoc_type<Float>
{
    adhoc_type result(std::erfc(arg.value));
    if (arg.is_active()) {
        result.id = opcode<Float>::global_tape->generate_id();
        opcode<Float>::global_tape->record_unary(OpCode::ERFC, arg.id, result.id);
        opcode<Float>::global_tape->record_value(arg.value);
    }
    return result;
}

template<class Float>
inline auto
cos(const adhoc_type<Float>& arg) -> adhoc_type<Float>
{
    adhoc_type result(std::cos(arg.value));
    if (arg.is_active()) {
        result.id = opcode<Float>::global_tape->generate_id();
        opcode<Float>::global_tape->record_unary(OpCode::COS, arg.id, result.id);
        opcode<Float>::global_tape->record_value(arg.value);
        opcode<Float>::global_tape->record_value(result.value);
    }
    return result;
}

template<class Float>
inline auto
norm(const adhoc_type<Float>& arg) -> adhoc_type<Float>
{
    adhoc_type result(std::norm(arg.value));
    if (arg.is_active()) {
        result.id = opcode<Float>::global_tape->generate_id();
        opcode<Float>::global_tape->record_unary(OpCode::NORM, arg.id, result.id);
        opcode<Float>::global_tape->record_value(arg.value);
    }
    return result;
}

template<class Float>
inline auto
inv(const adhoc_type<Float>& arg) -> adhoc_type<Float>
{
    adhoc_type result(1.0 / arg.value);
    if (arg.is_active()) {
        result.id = opcode<Float>::global_tape->generate_id();
        opcode<Float>::global_tape->record_unary(OpCode::INV, arg.id, result.id);
        opcode<Float>::global_tape->record_value(result.value);
    }
    return result;
}

template<class Float>
inline auto
abs(const adhoc_type<Float>& arg) -> adhoc_type<Float>
{
    adhoc_type result(std::abs(arg.value));
    if (arg.is_active()) {
        result.id = opcode<Float>::global_tape->generate_id();
        opcode<Float>::global_tape->record_unary(OpCode::ABS, arg.id, result.id);
        opcode<Float>::global_tape->record_value(arg.value);
    }
    return result;
}

template<class Float>
inline auto
sqrt(const adhoc_type<Float>& arg) -> adhoc_type<Float>
{
    adhoc_type result(std::sqrt(arg.value));
    if (arg.is_active()) {
        result.id = opcode<Float>::global_tape->generate_id();
        opcode<Float>::global_tape->record_unary(OpCode::SQRT, arg.id, result.id);
        opcode<Float>::global_tape->record_value(arg.value);
        opcode<Float>::global_tape->record_value(result.value);
    }
    return result;
}

template<class Float, class T>
inline auto
pow(const adhoc_type<Float>& lhs, T rhs) -> adhoc_type<Float>
{
    adhoc_type result(std::pow(lhs.get_value(), rhs));
    if (lhs.is_active()) {
        result.id = opcode<Float>::global_tape->generate_id();
        opcode<Float>::global_tape->record_unary(OpCode::POW_C, lhs.id, result.id);
        opcode<Float>::global_tape->record_value(lhs.value);
        opcode<Float>::global_tape->record_value(result.value);
        opcode<Float>::global_tape->record_value(static_cast<Float>(rhs));
    }
    return result;
}

template<class Float>
inline auto
pow(const adhoc_type<Float>& lhs, const adhoc_type<Float>& rhs) -> adhoc_type<Float>
{
    if (lhs.is_passive() && rhs.is_passive()) {
        return adhoc_type<Float>{ std::pow(lhs.get_value(), rhs.get_value()) };
    }
    // TODO: derivative
    throw;
    return adhoc_type<Float>{ std::pow(lhs.get_value(), rhs.get_value()) };
}

template<class T, class Float>
inline auto
pow(T lhs, const adhoc_type<Float>& rhs) -> adhoc_type<Float>
{
    // TODO: derivative
    throw;
    return adhoc_type<Float>{ std::pow(lhs, rhs.get_value()) };
}

template<class Float>
inline auto
atan2(const adhoc_type<Float>& lhs, const adhoc_type<Float>& rhs) -> adhoc_type<Float>
{
    // TODO: derivative
    throw;
    return adhoc_type<Float>{ std::atan2(lhs.get_value(), rhs.get_value()) };
}

template<class T, class Float>
inline auto
atan2(T lhs, const adhoc_type<Float>& rhs) -> adhoc_type<Float>
{
    // TODO: derivative
    throw;
    return adhoc_type<Float>{ std::atan2(lhs, rhs.get_value()) };
}

template<class Float, class T>
inline auto
atan2(const adhoc_type<Float>& lhs, T rhs) -> adhoc_type<Float>
{
    // TODO: derivative
    throw;
    return adhoc_type<Float>{ std::atan2(lhs.get_value(), rhs) };
}

template<class Float>
inline auto
max(const adhoc_type<Float>& lhs, const adhoc_type<Float>& rhs) -> adhoc_type<Float>
{
    if (lhs.get_value() >= rhs.get_value()) {
        return lhs;
    }

    return rhs;
}

template<class T, class Float>
inline auto
max(T lhs, const adhoc_type<Float>& rhs) -> adhoc_type<Float>
{
    // TODO: derivative
    throw;
    return adhoc_type<Float>{ std::max(lhs, rhs.get_value()) };
}

template<class Float, class T>
inline auto
max(const adhoc_type<Float>& lhs, T rhs) -> adhoc_type<Float>
{
    // TODO: derivative
    throw;
    return adhoc_type<Float>{ std::max(lhs.get_value(), rhs) };
}

template<class Float>
inline auto
min(const adhoc_type<Float>& lhs, const adhoc_type<Float>& rhs) -> adhoc_type<Float>
{
    if (lhs.get_value() <= rhs.get_value()) {
        return lhs;
    }

    return rhs;
}

template<class T, class Float>
inline auto
min(T lhs, const adhoc_type<Float>& rhs) -> adhoc_type<Float>
{
    // TODO: derivative
    throw;
    return adhoc_type<Float>{ std::min(lhs, rhs.get_value()) };
}

template<class Float, class T>
inline auto
min(const adhoc_type<Float>& lhs, T rhs) -> adhoc_type<Float>
{
    // TODO: derivative
    throw;
    return adhoc_type<Float>{ std::min(lhs.get_value(), rhs) };
}

template<class Float>
inline auto
fabs(const adhoc_type<Float>& arg) -> adhoc_type<Float>
{
    // TODO: derivative
    throw;
    return adhoc_type<Float>{ std::abs(arg.get_value()) };
}

template<class Float>
inline auto
floor(const adhoc_type<Float>& arg) -> adhoc_type<Float>
{
    // TODO: derivative
    return adhoc_type<Float>{ std::floor(arg.get_value()) };
}

template<class Float>
inline auto
sin(const adhoc_type<Float>& arg) -> adhoc_type<Float>
{
    // TODO: derivative
    throw;
    return adhoc_type<Float>{ std::sin(arg.get_value()) };
}

template<class Float>
inline auto
cosh(const adhoc_type<Float>& arg) -> adhoc_type<Float>
{
    // TODO: derivative
    throw;
    return adhoc_type<Float>{ std::cosh(arg.get_value()) };
}

template<class Float>
inline auto
sinh(const adhoc_type<Float>& arg) -> adhoc_type<Float>
{
    // TODO: derivative
    throw;
    return adhoc_type<Float>{ std::sinh(arg.get_value()) };
}

template<class Float>
inline auto
asin(const adhoc_type<Float>& arg) -> adhoc_type<Float>
{
    // TODO: derivative
    throw;
    return adhoc_type<Float>{ std::asin(arg.get_value()) };
}

template<class Float>
inline auto
isfinite(const adhoc_type<Float>& arg) -> bool
{
    return std::isfinite(arg.get_value());
}

template<class Float>
auto
adhoc_type<Float>::sub_c(double other) const -> adhoc_type
{
    adhoc_type result(other - this->value);
    if (this->is_active()) {
        result.id = opcode<Float>::global_tape->generate_id();
        opcode<Float>::global_tape->record_unary(OpCode::SUB_C, this->id, result.id);
    }
    return result;
}

template<class Float>
auto
adhoc_type<Float>::div_c(double other) const -> adhoc_type
{
    adhoc_type result(other / this->value);

    if (this->is_active()) {
        adhoc_type intermediary_result(1.0 / this->value);
        intermediary_result.id = opcode<Float>::global_tape->generate_id();
        result.id = opcode<Float>::global_tape->generate_id();
        opcode<Float>::global_tape->record_unary(OpCode::INV, this->id, intermediary_result.id);
        opcode<Float>::global_tape->record_value(intermediary_result.value);
        opcode<Float>::global_tape->record_value(other);
        opcode<Float>::global_tape->record_unary(OpCode::MUL_C, intermediary_result.id, result.id);
    }

    return result;
}

template<class Float>
inline auto
operator+(double lhs, const adhoc_type<Float>& rhs) -> adhoc_type<Float>
{
    return rhs + lhs;
}

template<class Float>
inline auto
operator-(double lhs, const adhoc_type<Float>& rhs) -> adhoc_type<Float>
{
    return rhs.sub_c(lhs);
}

template<class Float>
inline auto
operator*(double lhs, const adhoc_type<Float>& rhs) -> adhoc_type<Float>
{
    return rhs * lhs;
}

template<class Float>
inline auto
operator/(double lhs, const adhoc_type<Float>& rhs) -> adhoc_type<Float>
{
    return rhs.div_c(lhs);
}

template<class T, class Float>
inline auto
operator<(T lhs, const adhoc_type<Float>& rhs) -> bool
{
    return lhs < rhs.get_value();
}

template<class Float, class T>
inline auto
operator<(const adhoc_type<Float>& lhs, T rhs) -> bool
{
    return lhs.get_value() < rhs;
}

template<class Float1, class Float2>
inline auto
operator<(const adhoc_type<Float1>& lhs, const adhoc_type<Float2>& rhs) -> bool
{
    return lhs.get_value() < rhs.get_value();
}

template<class T, class Float>
inline auto
operator>(T lhs, const adhoc_type<Float>& rhs) -> bool
{
    return lhs > rhs.get_value();
}

template<class Float, class T>
inline auto
operator>(const adhoc_type<Float>& lhs, T rhs) -> bool
{
    return lhs.get_value() > rhs;
}

template<class Float1, class Float2>
inline auto
operator>(const adhoc_type<Float1>& lhs, const adhoc_type<Float2>& rhs) -> bool
{
    return lhs.get_value() > rhs.get_value();
}

template<class T, class Float>
inline auto
operator<=(T lhs, const adhoc_type<Float>& rhs) -> bool
{
    return lhs <= rhs.get_value();
}

template<class Float, class T>
inline auto
operator<=(const adhoc_type<Float>& lhs, T rhs) -> bool
{
    return lhs.get_value() <= rhs;
}

template<class Float1, class Float2>
inline auto
operator<=(const adhoc_type<Float1>& lhs, const adhoc_type<Float2>& rhs) -> bool
{
    return lhs.get_value() <= rhs.get_value();
}

template<class T, class Float>
inline auto
operator>=(T lhs, const adhoc_type<Float>& rhs) -> bool
{
    return lhs >= rhs.get_value();
}

template<class Float, class T>
inline auto
operator>=(const adhoc_type<Float>& lhs, T rhs) -> bool
{
    return lhs.get_value() >= rhs;
}

template<class Float1, class Float2>
inline auto
operator>=(const adhoc_type<Float1>& lhs, const adhoc_type<Float2>& rhs) -> bool
{
    return lhs.get_value() >= rhs.get_value();
}

template<class T, class Float>
inline auto
operator==(T lhs, const adhoc_type<Float>& rhs) -> bool
{
    return lhs == rhs.get_value();
}

template<class Float, class T>
inline auto
operator==(const adhoc_type<Float>& lhs, T rhs) -> bool
{
    return lhs.get_value() == rhs;
}

template<class Float1, class Float2>
inline auto
operator==(const adhoc_type<Float1>& lhs, const adhoc_type<Float2>& rhs) -> bool
{
    return lhs.get_value() == rhs.get_value();
}

template<class Float>
static inline auto
operator<<(std::ostream& out, const adhoc_type<Float>& x) -> std::ostream&
{
    out << x.get_value();
    return out;
}

// smart pointer that manages the lifetime of the static tape. It has an internal counter to track how many instances
// are using the tape, and only deletes the tape when the last instance is destroyed. This allows for multiple instances
// of adhoc to share the same tape without worrying about ownership issues.
template<typename mode_t>
struct smart_tape_ptr_t {
  private:
    using tape_t = typename mode_t::tape_t;
    inline static thread_local std::size_t ref_count = 0;

  public:
    smart_tape_ptr_t()
    {
        if (ref_count == 0) {
            mode_t::global_tape = new tape_t;
        }
        ++ref_count;
    }

    smart_tape_ptr_t(const smart_tape_ptr_t& /*other*/) { ++ref_count; }

    smart_tape_ptr_t(smart_tape_ptr_t&& /*other*/) noexcept { ++ref_count; }

    auto operator=(const smart_tape_ptr_t& other) -> smart_tape_ptr_t&
    {
        if (this != &other) {
            // No change in ref_count since we're already tracking this instance
        }
        return *this;
    }

    auto operator=(smart_tape_ptr_t&& /*other*/) noexcept->smart_tape_ptr_t&
    {
        // No change in ref_count since we're already tracking this instance
        return *this;
    }

    ~smart_tape_ptr_t()
    {
        --ref_count;
        if (ref_count == 0) {
            delete mode_t::global_tape;
            mode_t::global_tape = nullptr;
        }
    }

    auto get() const -> tape_t* { return mode_t::global_tape; }

    auto operator->() const -> tape_t* { return mode_t::global_tape; }

    auto operator*() const -> tape_t& { return *mode_t::global_tape; }

    explicit operator bool() const { return mode_t::global_tape != nullptr; }

    static auto use_count() -> std::size_t { return ref_count; }
};

template<typename T>
auto
size_of(const Tape<T>& arg, bool capacity = false) -> std::size_t
{
    return arg.size_of(capacity);
}

template<typename T>
auto
size_of(const Tape<T>* arg, bool capacity = false) -> std::size_t
{
    return arg->size_of(capacity);
}

template<typename T>
auto
size_of(const adhoc::smart_tape_ptr_t<adhoc::opcode<T> >& arg, bool capacity = false) -> std::size_t
{
    return arg->size_of(capacity);
}

template<class Float>
inline auto
isnan(adhoc::adhoc_type<Float> const& x) -> bool
{
    return std::isnan(x.get_value());
}

template<class Float>
inline auto
lround(adhoc::adhoc_type<Float> const& x) -> long
{
    return std::lround(x.get_value());
}

template<class Float>
inline auto
passive_value(const adhoc::adhoc_type<Float>& x) -> double
{
    return x.get_value();
}

template<class Float>
inline auto
passive_value(adhoc::adhoc_type<Float>& x) -> double&
{
    return x.get_value();
}

template<class Float>
inline auto
passive_value(adhoc::adhoc_type<Float>&& x) -> double
{
    return x.get_value();
}

template<class Float>
inline auto
passive_value(std::vector<adhoc::adhoc_type<Float> >& x) -> std::vector<Float>
{
    std::vector<Float> result;
    result.reserve(x.size());
    for (auto const& item : x) {
        result.push_back(item.get_value());
    }
    return result;
}

template<class Float>
inline auto
passive_value(std::vector<adhoc::adhoc_type<Float> > const& x) -> std::vector<Float>
{
    std::vector<Float> result;
    result.reserve(x.size());
    for (auto const& item : x) {
        result.push_back(item.get_value());
    }
    return result;
}

} // namespace adhoc

#ifdef ADHOC_HEADER_ONLY
#include "adhoc.cpp"
#endif

#endif // ADHOC_HPP
