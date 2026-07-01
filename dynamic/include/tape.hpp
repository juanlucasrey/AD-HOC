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

#ifndef TAPE_HPP
#define TAPE_HPP

#include "adhoc.hpp"

#include <memory>

namespace adhoc {

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
    FirstOrderVLossyCompressed,
    FirstOrderLossyPathReuse,
    FirstOrderVLossyPathReuse,
    FirstOrderLossyCompressedPathReuse,
    FirstOrderVLossyCompressedPathReuse,
    FirstOrderLossyCompressedPathReuseV,
    FirstOrderVLossyCompressedPathReuseV,
    SecondOrderLossy,
    SecondOrderVLossy,
};

struct PositionImpl;

template<class Float, class TapeDataType>
class Tape {
  private:
    struct Impl;

    std::unique_ptr<Impl> impl;
    TapeDataType& data;

    void reserve_input(std::size_t n_registered);
    void reserve_output(std::size_t n_registered);

  public:
    void register_variable(adhoc_type<Float, TapeDataType> const& var);
    void register_variable(adhoc_type<Float, TapeDataType>& var);

    template<class Container>
    void register_variable(Container&& vars)
    {
        std::size_t const count = std::distance(std::begin(vars), std::end(vars));
        this->reserve_input(count);
        for (auto&& var : vars) {
            register_variable(var);
        }
    }

    void register_output_variable(adhoc_type<Float, TapeDataType> const& var);
    void register_output_variable(adhoc_type<Float, TapeDataType>& var);

    template<class Container>
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

    Tape(TapeDataType& tape_data);
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

    struct position_t {
        position_t();
        ~position_t();
        auto operator=(position_t other) -> position_t&;
        position_t(const position_t& other);

      private:
        friend Tape;
        std::unique_ptr<PositionImpl> impl;
    };

    void set_checkpoint();
    void backpropagate();
    void backpropagate_to(position_t const& to);

    template<bool ResetInPlace = false, bool Log = false>
    void backpropagate_and_reset_to(position_t const& to);
    void set_derivative(adhoc_type<Float, TapeDataType> const& var, double deriv, std::size_t lane = 0);
    auto get_derivative(adhoc_type<Float, TapeDataType> const& var, std::size_t lane = 0) const -> double;

    auto get_derivative(adhoc_type<Float, TapeDataType> const& var1,
                        adhoc_type<Float, TapeDataType> const& var2,
                        std::size_t lane = 0) const -> double;
    void zero_adjoints();

    auto get_position() const -> position_t;

    void print() const;

    auto size_of(bool capacity = false) const -> std::size_t;
};

// smart pointer that manages the lifetime of the static tape. It has an internal counter to track how many instances
// are using the tape, and only deletes the tape when the last instance is destroyed. This allows for multiple instances
// of adhoc to share the same tape without worrying about ownership issues.
template<class mode_t>
class smart_tape_ptr_t {
  private:
    using tape_t = typename mode_t::tape_t;
    inline static thread_local std::size_t ref_count = 0;

  public:
    smart_tape_ptr_t()
    {
        if (ref_count == 0) {
            mode_t::global_tape_data = new typename mode_t::tape_data_t;
            mode_t::global_tape = new tape_t(*mode_t::global_tape_data);
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
            delete mode_t::global_tape_data;
            mode_t::global_tape_data = nullptr;
        }
    }

    auto get() const -> tape_t* { return mode_t::global_tape; }

    auto operator->() const -> tape_t* { return mode_t::global_tape; }

    auto operator*() const -> tape_t& { return *mode_t::global_tape; }

    explicit operator bool() const { return mode_t::global_tape != nullptr; }

    static auto use_count() -> std::size_t { return ref_count; }
};

template<class Float, class TapeDataType>
auto
size_of(const Tape<Float, TapeDataType>& arg, bool capacity = false) -> std::size_t
{
    return arg.size_of(capacity);
}

template<class Float, class TapeDataType>
auto
size_of(const Tape<Float, TapeDataType>* arg, bool capacity = false) -> std::size_t
{
    return arg->size_of(capacity);
}

template<class mode_t>
auto
size_of(const adhoc::smart_tape_ptr_t<mode_t>& arg, bool capacity = false) -> std::size_t
{
    return arg->size_of(capacity);
}

} // namespace adhoc

#ifdef ADHOC_HEADER_ONLY
#include "tape.cpp"
#endif

#endif // TAPE_HPP
