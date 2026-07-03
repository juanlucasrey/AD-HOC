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

    double value{ 0. };
    mutable std::size_t id{ passive_id<std::size_t> };

  public:
    using tape_data_t = TapeDataType;
    adhoc_type() = default;

    adhoc_type(Float val)
      : value(val)
    {
    }

    adhoc_type(const adhoc_type& other)
      : value(other.value)
      , id(other.id)
    {
    }

    adhoc_type(const adhoc_type&& other) noexcept { *this = other; }

    // Get value
    auto get_value() const -> double { return value; }
    auto get_value() -> double& { return value; }
    auto is_passive() const -> bool { return id == passive_id<std::size_t>; }
    auto is_active() const -> bool { return id != passive_id<std::size_t>; }

    auto operator-() const -> adhoc_type { return 0.0 - *this; }

    auto operator+(const adhoc_type& other) const -> adhoc_type
    {
        if (this->is_passive() && other.is_active()) {
            return other + this->value;
        }

        if (other.is_passive()) {
            return (*this) + other.value;
        }

        adhoc_type result(this->value + other.value);
        result.id = mode_t::global_tape_data->generate_id();
        if (this->id == other.id) {
            mode_t::global_tape_data->record_unary(OpCode::MUL_C, this->id, result.id);
            mode_t::global_tape_data->record_value(2.0);
        }
        else {
            mode_t::global_tape_data->record_binary(OpCode::ADD, this->id, other.id, result.id);
        }
        return result;
    }

    auto operator-(const adhoc_type& other) const -> adhoc_type
    {
        if (this->is_passive() && other.is_active()) {
            return this->value - other;
        }

        if (other.is_passive()) {
            return (*this) - other.value;
        }

        adhoc_type result(this->value - other.value);
        result.id = mode_t::global_tape_data->generate_id();
        if (this->id == other.id) {
            mode_t::global_tape_data->record_unary(OpCode::MUL_C, this->id, result.id);
            mode_t::global_tape_data->record_value(0.0);
        }
        else {
            mode_t::global_tape_data->record_binary(OpCode::SUB, this->id, other.id, result.id);
        }
        return result;
    }

    auto operator*(const adhoc_type& other) const -> adhoc_type
    {
        if (this->is_passive() && other.is_active()) {
            return other * this->value;
        }

        if (other.is_passive()) {
            return (*this) * other.value;
        }

        adhoc_type result(this->value * other.value);
        result.id = mode_t::global_tape_data->generate_id();
        if (this->id == other.id) {
            mode_t::global_tape_data->record_unary(OpCode::NORM, this->id, result.id);
            mode_t::global_tape_data->record_value(this->value);
        }
        else {
            mode_t::global_tape_data->record_binary(OpCode::MUL, this->id, other.id, result.id);
            mode_t::global_tape_data->record_value(this->value);
            mode_t::global_tape_data->record_value(other.value);
        }

        return result;
    }

    auto operator/(const adhoc_type& other) const -> adhoc_type
    {

        if (this->is_passive()) {
            return this->value / other;
        }

        if (other.is_passive()) {
            return *this / other.value;
        }

        adhoc_type result(this->value / other.value);

        adhoc_type intermediary_result(1.0 / other.value);
        intermediary_result.id = mode_t::global_tape_data->generate_id();
        result.id = mode_t::global_tape_data->generate_id();
        mode_t::global_tape_data->record_unary(OpCode::INV, other.id, intermediary_result.id);
        mode_t::global_tape_data->record_value(intermediary_result.value);
        mode_t::global_tape_data->record_binary(OpCode::MUL, this->id, intermediary_result.id, result.id);
        mode_t::global_tape_data->record_value(this->value);
        mode_t::global_tape_data->record_value(intermediary_result.value);
        return result;
    }

    auto operator+(double other) const -> adhoc_type
    {
        adhoc_type result(this->value + other);
        if (this->is_active()) {
            result.id = mode_t::global_tape_data->generate_id();
            mode_t::global_tape_data->record_unary(OpCode::ADD_C, this->id, result.id);
        }
        return result;
    }

    auto operator-(double other) const -> adhoc_type
    {
        adhoc_type result(this->value - other);
        if (this->is_active()) {
            result.id = mode_t::global_tape_data->generate_id();
            mode_t::global_tape_data->record_unary(OpCode::ADD_C, this->id, result.id);
        }
        return result;
    }

    auto operator*(double other) const -> adhoc_type
    {
        adhoc_type result(this->value * other);
        if (this->is_active()) {
            result.id = mode_t::global_tape_data->generate_id();
            mode_t::global_tape_data->record_unary(OpCode::MUL_C, this->id, result.id);
            mode_t::global_tape_data->record_value(other);
        }
        return result;
    }

    auto operator/(double other) const -> adhoc_type
    {
        adhoc_type result(this->value / other);
        if (this->is_active()) {
            result.id = mode_t::global_tape_data->generate_id();
            mode_t::global_tape_data->record_unary(OpCode::MUL_C, this->id, result.id);
            mode_t::global_tape_data->record_value(1.0 / other);
        }
        return result;
    }

    auto operator=(const adhoc_type& other) -> adhoc_type& = default;

    auto operator+=(const adhoc_type& other) -> adhoc_type&
    {
        if (this->is_passive() && other.is_active()) {
            return *this = other + this->value;
        }

        if (other.is_passive()) {
            return *this += other.value;
        }

        return *this = *this + other;
    }

    auto operator-=(const adhoc_type& other) -> adhoc_type&
    {
        if (this->is_passive() && other.is_active()) {
            return *this = this->value - other;
        }

        if (other.is_passive()) {
            return (*this) -= other.value;
        }

        return *this = *this - other;
    }

    auto operator*=(const adhoc_type& other) -> adhoc_type&
    {
        if (this->is_passive() && other.is_active()) {
            return *this = other * this->value;
        }

        if (other.is_passive()) {
            return (*this) *= other.value;
        }

        return *this = *this * other;
    }

    auto operator/=(const adhoc_type& other) -> adhoc_type& { return *this = *this / other; }
    auto operator+=(double other) -> adhoc_type& { return *this = *this + other; }
    auto operator-=(double other) -> adhoc_type& { return *this = *this - other; }
    auto operator*=(double other) -> adhoc_type& { return *this = *this * other; }
    auto operator/=(double other) -> adhoc_type& { return *this = *this / other; }

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

    // friend auto inv(const adhoc_type& arg) -> adhoc_type
    // {
    //     adhoc_type result(1.0 / arg.value);
    //     if (arg.is_active()) {
    //         result.id = mode_t::global_tape_data->generate_id();
    //         mode_t::global_tape_data->record_unary(OpCode::INV, arg.id, result.id);
    //         mode_t::global_tape_data->record_value(result.value);
    //     }
    //     return result;
    // }

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

    // friend auto pow(const adhoc_type<Float>& lhs, const adhoc_type<Float>& rhs) -> adhoc_type<Float>;

    // template<class T>
    // inline auto pow(T /* lhs */, const adhoc_type<Float>& /* rhs */) -> adhoc_type<Float>;

    friend auto pow(const adhoc_type& lhs, Float rhs) -> adhoc_type
    {
        adhoc_type result(std::pow(lhs.get_value(), rhs));
        if (lhs.is_active()) {
            result.id = mode_t::global_tape_data->generate_id();
            mode_t::global_tape_data->record_unary(OpCode::POW_C, lhs.id, result.id);
            mode_t::global_tape_data->record_value(lhs.value);
            mode_t::global_tape_data->record_value(rhs);
        }
        return result;
    }

    friend auto operator-(double lhs, const adhoc_type& rhs) -> adhoc_type
    {
        adhoc_type result(lhs - rhs.value);
        if (rhs.is_active()) {
            result.id = mode_t::global_tape_data->generate_id();
            mode_t::global_tape_data->record_unary(OpCode::SUB_C, rhs.id, result.id);
        }
        return result;
    }

    friend auto operator/(double lhs, const adhoc_type& rhs) -> adhoc_type
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
};

template<class Float, class TapeDataType>
inline auto
operator+(double lhs, const adhoc_type<Float, TapeDataType>& rhs) -> adhoc_type<Float, TapeDataType>
{
    return rhs + lhs;
}

template<class Float, class TapeDataType>
inline auto
operator*(double lhs, const adhoc_type<Float, TapeDataType>& rhs) -> adhoc_type<Float, TapeDataType>
{
    return rhs * lhs;
}

template<class T, class Float, class TapeDataType>
inline auto
operator<(T lhs, const adhoc_type<Float, TapeDataType>& rhs) -> bool
{
    return lhs < rhs.get_value();
}

template<class Float, class TapeDataType, class T>
inline auto
operator<(const adhoc_type<Float, TapeDataType>& lhs, T rhs) -> bool
{
    return lhs.get_value() < rhs;
}

template<class Float1, class TapeDataType1, class Float2, class TapeDataType2>
inline auto
operator<(const adhoc_type<Float1, TapeDataType1>& lhs, const adhoc_type<Float2, TapeDataType2>& rhs) -> bool
{
    return lhs.get_value() < rhs.get_value();
}

template<class T, class Float, class TapeDataType>
inline auto
operator>(T lhs, const adhoc_type<Float, TapeDataType>& rhs) -> bool
{
    return lhs > rhs.get_value();
}

template<class Float, class TapeDataType, class T>
inline auto
operator>(const adhoc_type<Float, TapeDataType>& lhs, T rhs) -> bool
{
    return lhs.get_value() > rhs;
}

template<class Float1, class TapeDataType1, class Float2, class TapeDataType2>
inline auto
operator>(const adhoc_type<Float1, TapeDataType1>& lhs, const adhoc_type<Float2, TapeDataType2>& rhs) -> bool
{
    return lhs.get_value() > rhs.get_value();
}

template<class T, class Float, class TapeDataType>
inline auto
operator<=(T lhs, const adhoc_type<Float, TapeDataType>& rhs) -> bool
{
    return lhs <= rhs.get_value();
}

template<class Float, class TapeDataType, class T>
inline auto
operator<=(const adhoc_type<Float, TapeDataType>& lhs, T rhs) -> bool
{
    return lhs.get_value() <= rhs;
}

template<class Float1, class TapeDataType1, class Float2, class TapeDataType2>
inline auto
operator<=(const adhoc_type<Float1, TapeDataType1>& lhs, const adhoc_type<Float2, TapeDataType2>& rhs) -> bool
{
    return lhs.get_value() <= rhs.get_value();
}

template<class T, class Float, class TapeDataType>
inline auto
operator>=(T lhs, const adhoc_type<Float, TapeDataType>& rhs) -> bool
{
    return lhs >= rhs.get_value();
}

template<class Float, class TapeDataType, class T>
inline auto
operator>=(const adhoc_type<Float, TapeDataType>& lhs, T rhs) -> bool
{
    return lhs.get_value() >= rhs;
}

template<class Float1, class TapeDataType1, class Float2, class TapeDataType2>
inline auto
operator>=(const adhoc_type<Float1, TapeDataType1>& lhs, const adhoc_type<Float2, TapeDataType2>& rhs) -> bool
{
    return lhs.get_value() >= rhs.get_value();
}

template<class T, class Float, class TapeDataType>
inline auto
operator==(T lhs, const adhoc_type<Float, TapeDataType>& rhs) -> bool
{
    return lhs == rhs.get_value();
}

template<class Float, class TapeDataType, class T>
inline auto
operator==(const adhoc_type<Float, TapeDataType>& lhs, T rhs) -> bool
{
    return lhs.get_value() == rhs;
}

template<class Float1, class TapeDataType1, class Float2, class TapeDataType2>
inline auto
operator==(const adhoc_type<Float1, TapeDataType1>& lhs, const adhoc_type<Float2, TapeDataType2>& rhs) -> bool
{
    return lhs.get_value() == rhs.get_value();
}

template<class Float, class TapeDataType>
static inline auto
operator<<(std::ostream& out, const adhoc_type<Float, TapeDataType>& x) -> std::ostream&
{
    out << x.get_value();
    return out;
}

template<class Float, class TapeDataType>
inline auto
isfinite(const adhoc_type<Float, TapeDataType>& arg) -> bool
{
    return std::isfinite(arg.get_value());
}

template<class Float, class TapeDataType>
inline auto
isnan(adhoc::adhoc_type<Float, TapeDataType> const& x) -> bool
{
    return std::isnan(x.get_value());
}

template<class Float, class TapeDataType>
inline auto
lround(adhoc::adhoc_type<Float, TapeDataType> const& x) -> long
{
    return std::lround(x.get_value());
}

template<class Float, class TapeDataType>
inline auto
passive_value(const adhoc::adhoc_type<Float, TapeDataType>& x) -> double
{
    return x.get_value();
}

template<class Float, class TapeDataType>
inline auto
passive_value(adhoc::adhoc_type<Float, TapeDataType>& x) -> double&
{
    return x.get_value();
}

template<class Float, class TapeDataType>
inline auto
passive_value(adhoc::adhoc_type<Float, TapeDataType>&& x) -> double
{
    return x.get_value();
}

template<class Float, class TapeDataType>
inline auto
passive_value(std::vector<adhoc::adhoc_type<Float, TapeDataType> >& x) -> std::vector<Float>
{
    std::vector<Float> result;
    result.reserve(x.size());
    for (auto const& item : x) {
        result.push_back(item.get_value());
    }
    return result;
}

template<class Float, class TapeDataType>
inline auto
passive_value(std::vector<adhoc::adhoc_type<Float, TapeDataType> > const& x) -> std::vector<Float>
{
    std::vector<Float> result;
    result.reserve(x.size());
    for (auto const& item : x) {
        result.push_back(item.get_value());
    }
    return result;
}

// unimplemented functions
template<class Float, class TapeDataType>
inline auto
pow(const adhoc_type<Float, TapeDataType>& lhs, const adhoc_type<Float, TapeDataType>& rhs)
  -> adhoc_type<Float, TapeDataType>
{
    if (lhs.is_passive() && rhs.is_passive()) {
        return adhoc_type<Float, TapeDataType>{ std::pow(lhs.get_value(), rhs.get_value()) };
    }
    // TODO: derivative
    throw;
    return adhoc_type<Float, TapeDataType>{ std::pow(lhs.get_value(), rhs.get_value()) };
}

template<class Float, class TapeDataType>
inline auto
pow(double lhs, const adhoc_type<Float, TapeDataType>& rhs) -> adhoc_type<Float, TapeDataType>
{
    // TODO: derivative
    throw;
    return adhoc_type<Float, TapeDataType>{ std::pow(lhs, rhs.get_value()) };
}

template<class Float, class TapeDataType>
inline auto
atan2(const adhoc_type<Float, TapeDataType>& lhs, const adhoc_type<Float, TapeDataType>& rhs)
  -> adhoc_type<Float, TapeDataType>
{
    // TODO: derivative
    throw;
    return adhoc_type<Float, TapeDataType>{ std::atan2(lhs.get_value(), rhs.get_value()) };
}

template<class T, class Float, class TapeDataType>
inline auto
atan2(T lhs, const adhoc_type<Float, TapeDataType>& rhs) -> adhoc_type<Float, TapeDataType>
{
    // TODO: derivative
    throw;
    return adhoc_type<Float, TapeDataType>{ std::atan2(lhs, rhs.get_value()) };
}

template<class Float, class TapeDataType, class T>
inline auto
atan2(const adhoc_type<Float, TapeDataType>& lhs, T rhs) -> adhoc_type<Float, TapeDataType>
{
    // TODO: derivative
    throw;
    return adhoc_type<Float, TapeDataType>{ std::atan2(lhs.get_value(), rhs) };
}

template<class Float, class TapeDataType>
inline auto
max(const adhoc_type<Float, TapeDataType>& lhs, const adhoc_type<Float, TapeDataType>& rhs)
  -> adhoc_type<Float, TapeDataType>
{
    if (lhs.get_value() >= rhs.get_value()) {
        return lhs;
    }

    return rhs;
}

template<class T, class Float, class TapeDataType>
inline auto
max(T lhs, const adhoc_type<Float, TapeDataType>& rhs) -> adhoc_type<Float, TapeDataType>
{
    // TODO: derivative
    throw;
    return adhoc_type<Float, TapeDataType>{ std::max(lhs, rhs.get_value()) };
}

template<class Float, class TapeDataType, class T>
inline auto
max(const adhoc_type<Float, TapeDataType>& lhs, T rhs) -> adhoc_type<Float, TapeDataType>
{
    // TODO: derivative
    throw;
    return adhoc_type<Float, TapeDataType>{ std::max(lhs.get_value(), rhs) };
}

template<class Float, class TapeDataType>
inline auto
min(const adhoc_type<Float, TapeDataType>& lhs, const adhoc_type<Float, TapeDataType>& rhs)
  -> adhoc_type<Float, TapeDataType>
{
    if (lhs.get_value() <= rhs.get_value()) {
        return lhs;
    }

    return rhs;
}

template<class T, class Float, class TapeDataType>
inline auto
min(T lhs, const adhoc_type<Float, TapeDataType>& rhs) -> adhoc_type<Float, TapeDataType>
{
    // TODO: derivative
    throw;
    return adhoc_type<Float, TapeDataType>{ std::min(lhs, rhs.get_value()) };
}

template<class Float, class TapeDataType, class T>
inline auto
min(const adhoc_type<Float, TapeDataType>& lhs, T rhs) -> adhoc_type<Float, TapeDataType>
{
    // TODO: derivative
    throw;
    return adhoc_type<Float, TapeDataType>{ std::min(lhs.get_value(), rhs) };
}

template<class Float, class TapeDataType>
inline auto
fabs(const adhoc_type<Float, TapeDataType>& arg) -> adhoc_type<Float, TapeDataType>
{
    // TODO: derivative
    throw;
    return adhoc_type<Float, TapeDataType>{ std::abs(arg.get_value()) };
}

template<class Float, class TapeDataType>
inline auto
floor(const adhoc_type<Float, TapeDataType>& arg) -> adhoc_type<Float, TapeDataType>
{
    // TODO: derivative
    return adhoc_type<Float, TapeDataType>{ std::floor(arg.get_value()) };
}

template<class Float, class TapeDataType>
inline auto
sin(const adhoc_type<Float, TapeDataType>& arg) -> adhoc_type<Float, TapeDataType>
{
    // TODO: derivative
    throw;
    return adhoc_type<Float, TapeDataType>{ std::sin(arg.get_value()) };
}

template<class Float, class TapeDataType>
inline auto
cosh(const adhoc_type<Float, TapeDataType>& arg) -> adhoc_type<Float, TapeDataType>
{
    // TODO: derivative
    throw;
    return adhoc_type<Float, TapeDataType>{ std::cosh(arg.get_value()) };
}

template<class Float, class TapeDataType>
inline auto
sinh(const adhoc_type<Float, TapeDataType>& arg) -> adhoc_type<Float, TapeDataType>
{
    // TODO: derivative
    throw;
    return adhoc_type<Float, TapeDataType>{ std::sinh(arg.get_value()) };
}

template<class Float, class TapeDataType>
inline auto
asin(const adhoc_type<Float, TapeDataType>& arg) -> adhoc_type<Float, TapeDataType>
{
    // TODO: derivative
    throw;
    return adhoc_type<Float, TapeDataType>{ std::asin(arg.get_value()) };
}

} // namespace adhoc

#endif // ADHOC_HPP
