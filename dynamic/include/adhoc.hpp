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
#include <vector>

#include <concepts>

#include "adhoc/passive_id.hpp"
#include "adhoc/tape_data.hpp"

namespace adhoc {

template<class type>
class Tape;

template<class Float, class TapeDataType>
class adhoc_type;

template<class mode_t>
class smart_tape_ptr_t;

template<class FloatPrimal,
         class FloatTape = FloatPrimal,
         EnumVectorType enumvectype = EnumVectorType::Simple,
         IdxVectorType idxvectype = IdxVectorType::Simple>
class opcode {
  public:
    using tape_data_t = TapeData<FloatTape, enumvectype, idxvectype>;
    using type = adhoc_type<FloatPrimal, tape_data_t>;
    using tape_t = Tape<type>;
    inline static thread_local tape_t* global_tape = nullptr;

    // windows doesn't like it when these are private, even though they are only used in the friend classes
#ifndef _MSC_VER
  private:
#endif
    friend type;
    friend class smart_tape_ptr_t<opcode>;
    inline static thread_local tape_data_t* global_tape_data = nullptr;
};

template<class Float, class TapeDataType>
class adhoc_type {
  private:
    friend Tape<adhoc_type<Float, TapeDataType> >;
    using mode_t =
      opcode<Float, typename TapeDataType::Float, TapeDataType::tape_enumvector_t, TapeDataType::tape_idxvector_t>;

    Float value{ 0. };
    mutable std::size_t id{ passive_id<std::size_t> };

  public:
    using tape_data_t = TapeDataType;
    adhoc_type() = default;

    adhoc_type(Float val)
      : value(val)
    {
    }

    // Get value
    auto get_value() const -> Float { return value; }
    auto get_value() -> Float& { return value; }
    auto is_passive() const -> bool { return id == passive_id<std::size_t>; }
    auto is_active() const -> bool { return id != passive_id<std::size_t>; }

    auto operator-() const -> adhoc_type { return 0.0 - *this; }

    auto operator+(const adhoc_type& rhs) const -> adhoc_type
    {
        if (this->is_passive() && rhs.is_active()) {
            return rhs + this->value;
        }

        if (rhs.is_passive()) {
            return (*this) + rhs.value;
        }

        adhoc_type result(this->value + rhs.value);
        result.id = mode_t::global_tape_data->generate_id();
        if (this->id == rhs.id) {
            mode_t::global_tape_data->record_unary(OpCode::MUL_C, this->id, result.id);
            mode_t::global_tape_data->record_value(2.0);
        }
        else {
            mode_t::global_tape_data->record_binary(OpCode::ADD, this->id, rhs.id, result.id);
        }
        return result;
    }

    auto operator-(const adhoc_type& rhs) const -> adhoc_type
    {
        if (this->is_passive() && rhs.is_active()) {
            return this->value - rhs;
        }

        if (rhs.is_passive()) {
            return (*this) - rhs.value;
        }

        adhoc_type result(this->value - rhs.value);
        result.id = mode_t::global_tape_data->generate_id();
        if (this->id == rhs.id) {
            mode_t::global_tape_data->record_unary(OpCode::MUL_C, this->id, result.id);
            mode_t::global_tape_data->record_value(0.0);
        }
        else {
            mode_t::global_tape_data->record_binary(OpCode::SUB, this->id, rhs.id, result.id);
        }
        return result;
    }

    auto operator*(const adhoc_type& rhs) const -> adhoc_type
    {
        if (this->is_passive() && rhs.is_active()) {
            return rhs * this->value;
        }

        if (rhs.is_passive()) {
            return (*this) * rhs.value;
        }

        adhoc_type result(this->value * rhs.value);
        result.id = mode_t::global_tape_data->generate_id();
        if (this->id == rhs.id) {
            mode_t::global_tape_data->record_unary(OpCode::NORM, this->id, result.id);
            mode_t::global_tape_data->record_value(this->value);
        }
        else {
            mode_t::global_tape_data->record_binary(OpCode::MUL, this->id, rhs.id, result.id);
            mode_t::global_tape_data->record_value(this->value);
            mode_t::global_tape_data->record_value(rhs.value);
        }

        return result;
    }

    auto operator/(const adhoc_type& rhs) const -> adhoc_type
    {

        if (this->is_passive()) {
            return this->value / rhs;
        }

        if (rhs.is_passive()) {
            return *this / rhs.value;
        }

        adhoc_type result(this->value / rhs.value);

        adhoc_type intermediary_result(1.0 / rhs.value);
        intermediary_result.id = mode_t::global_tape_data->generate_id();
        result.id = mode_t::global_tape_data->generate_id();
        mode_t::global_tape_data->record_unary(OpCode::INV, rhs.id, intermediary_result.id);
        mode_t::global_tape_data->record_value(intermediary_result.value);
        mode_t::global_tape_data->record_binary(OpCode::MUL, this->id, intermediary_result.id, result.id);
        mode_t::global_tape_data->record_value(this->value);
        mode_t::global_tape_data->record_value(intermediary_result.value);
        return result;
    }

    template<class T>
    auto operator+(T rhs) const -> adhoc_type
    {
        adhoc_type result(this->value + rhs);
        if (this->is_active()) {
            result.id = mode_t::global_tape_data->generate_id();
            mode_t::global_tape_data->record_unary(OpCode::ADD_C, this->id, result.id);
        }
        return result;
    }

    template<class T>
    auto operator-(T rhs) const -> adhoc_type
    {
        adhoc_type result(this->value - rhs);
        if (this->is_active()) {
            result.id = mode_t::global_tape_data->generate_id();
            mode_t::global_tape_data->record_unary(OpCode::ADD_C, this->id, result.id);
        }
        return result;
    }

    template<class T>
    auto operator*(T rhs) const -> adhoc_type
    {
        adhoc_type result(this->value * rhs);
        if (this->is_active()) {
            result.id = mode_t::global_tape_data->generate_id();
            mode_t::global_tape_data->record_unary(OpCode::MUL_C, this->id, result.id);
            mode_t::global_tape_data->record_value(rhs);
        }
        return result;
    }

    template<class T>
    auto operator/(T rhs) const -> adhoc_type
    {
        adhoc_type result(this->value / rhs);
        if (this->is_active()) {
            result.id = mode_t::global_tape_data->generate_id();
            mode_t::global_tape_data->record_unary(OpCode::MUL_C, this->id, result.id);
            mode_t::global_tape_data->record_value(1.0 / rhs);
        }
        return result;
    }

    auto operator+=(const adhoc_type& arg) -> adhoc_type&
    {
        if (this->is_passive() && arg.is_active()) {
            return *this = arg + this->value;
        }

        if (arg.is_passive()) {
            return *this += arg.value;
        }

        return *this = *this + arg;
    }

    auto operator-=(const adhoc_type& arg) -> adhoc_type&
    {
        if (this->is_passive() && arg.is_active()) {
            return *this = this->value - arg;
        }

        if (arg.is_passive()) {
            return (*this) -= arg.value;
        }

        return *this = *this - arg;
    }

    auto operator*=(const adhoc_type& arg) -> adhoc_type&
    {
        if (this->is_passive() && arg.is_active()) {
            return *this = arg * this->value;
        }

        if (arg.is_passive()) {
            return (*this) *= arg.value;
        }

        return *this = *this * arg;
    }

    auto operator/=(const adhoc_type& arg) -> adhoc_type& { return *this = *this / arg; }
    template<class T>
    auto operator+=(T arg) -> adhoc_type&
    {
        return *this = *this + arg;
    }
    template<class T>
    auto operator-=(T arg) -> adhoc_type&
    {
        return *this = *this - arg;
    }
    template<class T>
    auto operator*=(T arg) -> adhoc_type&
    {
        return *this = *this * arg;
    }
    template<class T>
    auto operator/=(T arg) -> adhoc_type&
    {
        return *this = *this / arg;
    }

    // these functions need to be defined here, since they are templated
    // (templated friend functions are a dark corner of the standard)
    friend auto exp(const adhoc_type& arg) -> adhoc_type
    {
        adhoc_type result(std::exp(arg.value));
        if (arg.is_active()) {
            result.id = mode_t::global_tape_data->generate_id();
            mode_t::global_tape_data->record_unary(OpCode::EXP, arg.id, result.id);
            mode_t::global_tape_data->record_value(result.value);
        }
        return result;
    }

    friend auto expm1(const adhoc_type& arg) -> adhoc_type
    {
        adhoc_type result(std::expm1(arg.value));

        if (arg.is_active()) {
            adhoc_type intermediary_result(std::exp(arg.value));
            intermediary_result.id = mode_t::global_tape_data->generate_id();
            result.id = mode_t::global_tape_data->generate_id();
            mode_t::global_tape_data->record_unary(OpCode::EXP, arg.id, intermediary_result.id);
            mode_t::global_tape_data->record_value(intermediary_result.value);
            mode_t::global_tape_data->record_unary(OpCode::ADD_C, intermediary_result.id, result.id);
        }

        return result;
    }

    friend auto log(const adhoc_type& arg) -> adhoc_type
    {
        adhoc_type result(std::log(arg.value));
        if (arg.is_active()) {
            result.id = mode_t::global_tape_data->generate_id();
            mode_t::global_tape_data->record_unary(OpCode::LOG, arg.id, result.id);
            mode_t::global_tape_data->record_value(arg.value);
        }
        return result;
    }

    friend auto erf(const adhoc_type& arg) -> adhoc_type
    {
        adhoc_type result(std::erf(arg.value));
        if (arg.is_active()) {
            result.id = mode_t::global_tape_data->generate_id();
            mode_t::global_tape_data->record_unary(OpCode::ERF, arg.id, result.id);
            mode_t::global_tape_data->record_value(arg.value);
        }
        return result;
    }

    friend auto erfc(const adhoc_type& arg) -> adhoc_type
    {
        adhoc_type result(std::erfc(arg.value));
        if (arg.is_active()) {
            result.id = mode_t::global_tape_data->generate_id();
            mode_t::global_tape_data->record_unary(OpCode::ERFC, arg.id, result.id);
            mode_t::global_tape_data->record_value(arg.value);
        }
        return result;
    }

    friend auto cos(const adhoc_type& arg) -> adhoc_type
    {
        adhoc_type result(std::cos(arg.value));
        if (arg.is_active()) {
            result.id = mode_t::global_tape_data->generate_id();
            mode_t::global_tape_data->record_unary(OpCode::COS, arg.id, result.id);
            mode_t::global_tape_data->record_value(arg.value);
            mode_t::global_tape_data->record_value(result.value);
        }
        return result;
    }

    friend auto norm(const adhoc_type& arg) -> adhoc_type
    {
        adhoc_type result(std::norm(arg.value));
        if (arg.is_active()) {
            result.id = mode_t::global_tape_data->generate_id();
            mode_t::global_tape_data->record_unary(OpCode::NORM, arg.id, result.id);
            mode_t::global_tape_data->record_value(arg.value);
        }
        return result;
    }

    friend auto abs(const adhoc_type& arg) -> adhoc_type
    {
        adhoc_type result(std::abs(arg.value));
        if (arg.is_active()) {
            result.id = mode_t::global_tape_data->generate_id();
            mode_t::global_tape_data->record_unary(OpCode::ABS, arg.id, result.id);
            mode_t::global_tape_data->record_value(arg.value);
        }
        return result;
    }

    friend auto sqrt(const adhoc_type& arg) -> adhoc_type
    {
        adhoc_type result(std::sqrt(arg.value));
        if (arg.is_active()) {
            result.id = mode_t::global_tape_data->generate_id();
            mode_t::global_tape_data->record_unary(OpCode::SQRT, arg.id, result.id);
            mode_t::global_tape_data->record_value(arg.value);
            mode_t::global_tape_data->record_value(result.value);
        }
        return result;
    }

    friend auto pow(const adhoc_type& lhs, Float rhs) -> adhoc_type
    {
        adhoc_type result(std::pow(lhs.value, rhs));
        if (lhs.is_active()) {
            result.id = mode_t::global_tape_data->generate_id();
            mode_t::global_tape_data->record_unary(OpCode::POW_C, lhs.id, result.id);
            mode_t::global_tape_data->record_value(lhs.value);
            mode_t::global_tape_data->record_value(rhs);
        }
        return result;
    }

    template<class T>
    friend auto operator+(T lhs, const adhoc_type& rhs) -> adhoc_type
    {
        return rhs + lhs;
    }

    template<class T>
    friend auto operator-(T lhs, const adhoc_type& rhs) -> adhoc_type
    {
        adhoc_type result(lhs - rhs.value);
        if (rhs.is_active()) {
            result.id = mode_t::global_tape_data->generate_id();
            mode_t::global_tape_data->record_unary(OpCode::SUB_C, rhs.id, result.id);
        }
        return result;
    }

    template<std::floating_point T>
    friend auto operator*(T lhs, const adhoc_type& rhs) -> adhoc_type
    {
        return rhs * lhs;
    }

    template<class T>
    friend auto operator/(T lhs, const adhoc_type& rhs) -> adhoc_type
    {
        adhoc_type result(lhs / rhs.value);

        if (rhs.is_active()) {
            adhoc_type intermediary_result(1.0 / rhs.value);
            intermediary_result.id = mode_t::global_tape_data->generate_id();
            result.id = mode_t::global_tape_data->generate_id();
            mode_t::global_tape_data->record_unary(OpCode::INV, rhs.id, intermediary_result.id);
            mode_t::global_tape_data->record_value(intermediary_result.value);
            mode_t::global_tape_data->record_value(lhs);
            mode_t::global_tape_data->record_unary(OpCode::MUL_C, intermediary_result.id, result.id);
        }

        return result;
    }

    template<class T>
    friend auto operator<(T lhs, const adhoc_type& rhs) -> bool
    {
        return lhs < rhs.value;
    }

    template<class T>
    friend auto operator<(const adhoc_type& lhs, T rhs) -> bool
    {
        return lhs.value < rhs;
    }

    friend auto operator<(const adhoc_type& lhs, const adhoc_type& rhs) -> bool { return lhs.value < rhs.value; }

    template<class T>
    friend auto operator>(T lhs, const adhoc_type& rhs) -> bool
    {
        return lhs > rhs.value;
    }

    template<class T>
    friend auto operator>(const adhoc_type& lhs, T rhs) -> bool
    {
        return lhs.value > rhs;
    }

    friend auto operator>(const adhoc_type& lhs, const adhoc_type& rhs) -> bool { return lhs.value > rhs.value; }

    template<class T>
    friend auto operator<=(T lhs, const adhoc_type& rhs) -> bool
    {
        return lhs <= rhs.value;
    }

    template<class T>
    friend auto operator<=(const adhoc_type& lhs, T rhs) -> bool
    {
        return lhs.value <= rhs;
    }

    friend auto operator<=(const adhoc_type& lhs, const adhoc_type& rhs) -> bool { return lhs.value <= rhs.value; }

    template<class T>
    friend auto operator>=(T lhs, const adhoc_type& rhs) -> bool
    {
        return lhs >= rhs.value;
    }

    template<class T>
    friend auto operator>=(const adhoc_type& lhs, T rhs) -> bool
    {
        return lhs.value >= rhs;
    }

    friend auto operator>=(const adhoc_type& lhs, const adhoc_type& rhs) -> bool { return lhs.value >= rhs.value; }

    template<class T>
    friend auto operator==(T lhs, const adhoc_type& rhs) -> bool
    {
        return lhs == rhs.value;
    }

    template<class T>
    friend auto operator==(const adhoc_type& lhs, T rhs) -> bool
    {
        return lhs.value == rhs;
    }

    friend auto operator==(const adhoc_type& lhs, const adhoc_type& rhs) -> bool { return lhs.value == rhs.value; }

    friend auto operator<<(std::ostream& out, const adhoc_type& arg) -> std::ostream&
    {
        out << arg.value;
        return out;
    }

    friend auto isfinite(const adhoc_type& arg) -> bool { return std::isfinite(arg.value); }
    friend auto isnan(adhoc_type const& arg) -> bool { return std::isnan(arg.value); }
    friend auto lround(adhoc_type const& arg) -> long { return std::lround(arg.value); }

    // unimplemented functions
    friend auto pow(const adhoc_type& lhs, const adhoc_type& rhs) -> adhoc_type
    {
        if (lhs.is_passive() && rhs.is_passive()) {
            return adhoc_type{ std::pow(lhs.value, rhs.value) };
        }
        // TODO: derivative
        throw;
        return adhoc_type{ std::pow(lhs.value, rhs.value) };
    }

    template<class T>
    friend auto pow(T lhs, const adhoc_type& rhs) -> adhoc_type
    {
        // TODO: derivative
        throw;
        return adhoc_type{ std::pow(lhs, rhs.value) };
    }

    friend auto atan2(const adhoc_type& lhs, const adhoc_type& rhs) -> adhoc_type
    {
        // TODO: derivative
        throw;
        return adhoc_type{ std::atan2(lhs.value, rhs.value) };
    }

    template<class T>
    friend auto atan2(T lhs, const adhoc_type& rhs) -> adhoc_type
    {
        // TODO: derivative
        throw;
        return adhoc_type{ std::atan2(lhs, rhs.value) };
    }

    template<class T>
    friend auto atan2(const adhoc_type& lhs, T rhs) -> adhoc_type
    {
        // TODO: derivative
        throw;
        return adhoc_type{ std::atan2(lhs.value, rhs) };
    }

    friend auto max(const adhoc_type& lhs, const adhoc_type& rhs) -> adhoc_type
    {
        if (lhs.value >= rhs.value) {
            return lhs;
        }

        return rhs;
    }

    template<class T>
    friend auto max(T lhs, const adhoc_type& rhs) -> adhoc_type
    {
        // TODO: derivative
        throw;
        return adhoc_type{ std::max(lhs, rhs.value) };
    }

    template<class T>
    friend auto max(const adhoc_type& lhs, T rhs) -> adhoc_type
    {
        // TODO: derivative
        throw;
        return adhoc_type{ std::max(lhs.value, rhs) };
    }

    friend auto min(const adhoc_type& lhs, const adhoc_type& rhs) -> adhoc_type
    {
        if (lhs.value <= rhs.value) {
            return lhs;
        }

        return rhs;
    }

    template<class T>
    friend auto min(T lhs, const adhoc_type& rhs) -> adhoc_type
    {
        // TODO: derivative
        throw;
        return adhoc_type{ std::min(lhs, rhs.value) };
    }

    template<class T>
    friend auto min(const adhoc_type& lhs, T rhs) -> adhoc_type
    {
        // TODO: derivative
        throw;
        return adhoc_type{ std::min(lhs.value, rhs) };
    }

    friend auto fabs(const adhoc_type& arg) -> adhoc_type
    {
        // TODO: derivative
        throw;
        return adhoc_type{ std::abs(arg.value) };
    }

    friend auto floor(const adhoc_type& arg) -> adhoc_type
    {
        // TODO: derivative
        return adhoc_type{ std::floor(arg.value) };
    }

    friend auto sin(const adhoc_type& arg) -> adhoc_type
    {
        // TODO: derivative
        throw;
        return adhoc_type{ std::sin(arg.value) };
    }

    friend auto cosh(const adhoc_type& arg) -> adhoc_type
    {
        // TODO: derivative
        throw;
        return adhoc_type{ std::cosh(arg.value) };
    }

    friend auto sinh(const adhoc_type& arg) -> adhoc_type
    {
        // TODO: derivative
        throw;
        return adhoc_type{ std::sinh(arg.value) };
    }

    friend auto asin(const adhoc_type& arg) -> adhoc_type
    {
        // TODO: derivative
        throw;
        return adhoc_type{ std::asin(arg.value) };
    }
};

template<class Float, class TapeDataType>
inline auto
passive_value(const adhoc_type<Float, TapeDataType>& arg) -> Float
{
    return arg.get_value();
}

template<class Float, class TapeDataType>
inline auto
passive_value(adhoc_type<Float, TapeDataType>& arg) -> Float&
{
    return arg.get_value();
}

template<class Float, class TapeDataType>
inline auto
passive_value(adhoc_type<Float, TapeDataType>&& arg) -> Float
{
    return arg.get_value();
}

template<class Float, class TapeDataType>
inline auto
passive_value(std::vector<adhoc_type<Float, TapeDataType> >& arg) -> std::vector<Float>
{
    std::vector<Float> result;
    result.reserve(arg.size());
    for (auto const& item : arg) {
        result.push_back(item.get_value());
    }
    return result;
}

template<class Float, class TapeDataType>
inline auto
passive_value(std::vector<adhoc_type<Float, TapeDataType> > const& arg) -> std::vector<Float>
{
    std::vector<Float> result;
    result.reserve(arg.size());
    for (auto const& item : arg) {
        result.push_back(item.get_value());
    }
    return result;
}

} // namespace adhoc

#endif // ADHOC_HPP
