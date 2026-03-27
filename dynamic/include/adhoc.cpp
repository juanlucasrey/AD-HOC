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

#include "adhoc.hpp"

#include "adhoc/backpropagator1.hpp"
#include "adhoc/backpropagator1lossy.hpp"
#include "adhoc/backpropagator1lossycompressed.hpp"
#include "adhoc/backpropagator1vlossycompressed.hpp"
#include "adhoc/backpropagator2.hpp"
#include "adhoc/backpropagator2v.hpp"

#include <cstddef>
#include <iostream>
#include <variant>

namespace adhoc {

namespace {

void
record_register(std::vector<OpCode>& ops, std::vector<std::size_t>& ids, OpCode op, std::size_t arg_id)
{
    ops.push_back(op);
    ids.push_back(arg_id);
}

} // namespace

template<class Float>
struct Tape<Float>::Impl {
    std::variant<BackPropagator<Float>,
                 BackPropagator<Float, true>,
                 BackPropagator2<Float>,
                 BackPropagator2Simd8<Float, MapType::STD_MAP>,
                 BackPropagator2Simd8<Float, MapType::STD_UNORDERED_MAP>,
                 BackPropagator2Simd8<Float, MapType::ANKERL_UNORDERED_DENSE>,
                 BackPropagator2Simd8<Float, MapType::BOOST_UNORDERED_MAP>,
                 BackPropagatorLossy<Float>,
                 BackPropagatorLossy<Float, true>,
                 BackPropagatorLossyCompressed<Float>,
                 BackPropagatorVLossyCompressed<Float> >
      bp = BackPropagator<Float>();
};

template<class Float>
Tape<Float>::Tape()
  : impl(std::make_unique<Impl>())
{
}

template<class Float>
Tape<Float>::~Tape() = default;

template<class Float>
void
Tape<Float>::reserve_input(std::size_t count_registered)
{
    std::visit([count_registered](auto& arg) { arg.reserve_input(count_registered); }, this->impl->bp);
}

template<class Float>
void
Tape<Float>::reserve_output(std::size_t count_registered)
{
    std::visit([count_registered](auto& arg) { arg.reserve_output(count_registered); }, this->impl->bp);
}

template<class Float>
void
Tape<Float>::register_variable(adhoc_type<Float> const& var)
{
    if (var.is_passive()) {
        std::size_t const new_id = generate_id();
        var.id = new_id;
        record_register(data.ops, data.ids, OpCode::REG_INPUT, new_id);
        std::visit([new_id](auto& arg) { arg.register_variable(new_id); }, this->impl->bp);
    }
}

template<class Float>
void
Tape<Float>::register_variable(adhoc_type<Float>& var)
{
    if (var.is_passive()) {
        std::size_t const new_id = generate_id();
        var.id = new_id;
        record_register(data.ops, data.ids, OpCode::REG_INPUT, new_id);
        std::visit([new_id](auto& arg) { arg.register_variable(new_id); }, this->impl->bp);
    }
}

template<class Float>
void
Tape<Float>::register_output_variable(adhoc_type<Float> const& var)
{
    if (var.is_active()) {
        std::size_t const new_id = generate_id();
        data.ids.push_back(var.id);
        var.id = new_id;
        record_register(data.ops, data.ids, OpCode::REG_OUTPUT, new_id);
        std::visit([new_id, ops_size = data.ops.size()](auto& arg) { arg.register_output_variable(new_id, ops_size); },
                   this->impl->bp);
    }
}

template<class Float>
void
Tape<Float>::register_output_variable(adhoc_type<Float>& var)
{
    if (var.is_active()) {
        std::size_t const new_id = generate_id();
        data.ids.push_back(var.id);
        var.id = new_id;
        record_register(data.ops, data.ids, OpCode::REG_OUTPUT, new_id);
        std::visit([new_id, ops_size = data.ops.size()](auto& arg) { arg.register_output_variable(new_id, ops_size); },
                   this->impl->bp);
    }
}

template<class Float>
void
Tape<Float>::set_lanes(std::size_t num_lanes)
{
    std::visit([num_lanes](auto& arg) { arg.set_lanes(num_lanes); }, this->impl->bp);
}

template<class Float>
auto
Tape<Float>::get_lanes() const -> std::size_t
{
    return std::visit([](auto& arg) { return arg.get_lanes(); }, this->impl->bp);
}

template<class Float>
void
Tape<Float>::set_method(Method m)
{
    if (m == Method::FirstOrderSimple) {
        this->impl->bp.template emplace<BackPropagator<double> >();
    }
    else if (m == Method::FirstOrderSimd8) {
        this->impl->bp.template emplace<BackPropagator<double, true> >();
    }
    else if (m == Method::SecondOrderSimple) {
        this->impl->bp.template emplace<BackPropagator2<double> >();
    }
    else if (m == Method::SecondOrderSimd8_stdmap) {
        this->impl->bp.template emplace<BackPropagator2Simd8<double, MapType::STD_MAP> >();
    }
    else if (m == Method::SecondOrderSimd8_stdunorderedmap) {
        this->impl->bp.template emplace<BackPropagator2Simd8<double, MapType::STD_UNORDERED_MAP> >();
    }
    else if (m == Method::SecondOrderSimd8_ankerl) {
        this->impl->bp.template emplace<BackPropagator2Simd8<double, MapType::ANKERL_UNORDERED_DENSE> >();
    }
    else if (m == Method::SecondOrderSimd8_boost) {
        this->impl->bp.template emplace<BackPropagator2Simd8<double, MapType::BOOST_UNORDERED_MAP> >();
    }
    else if (m == Method::FirstOrderLossy) {
        this->impl->bp.template emplace<BackPropagatorLossy<double> >();
    }
    else if (m == Method::FirstOrderVLossy) {
        this->impl->bp.template emplace<BackPropagatorLossy<double, true> >();
    }
    else if (m == Method::FirstOrderLossyCompressed) {
        this->impl->bp.template emplace<BackPropagatorLossyCompressed<double> >();
    }
    else if (m == Method::FirstOrderVLossyCompressed) {
        this->impl->bp.template emplace<BackPropagatorVLossyCompressed<double> >();
    }
}

template<class Float>
auto
Tape<Float>::get_method() const -> Method
{
    if (std::holds_alternative<BackPropagator<Float> >(this->impl->bp)) {
        return Method::FirstOrderSimple;
    }

    if (std::holds_alternative<BackPropagator<Float, true> >(this->impl->bp)) {
        return Method::FirstOrderSimd8;
    }

    if (std::holds_alternative<BackPropagator2<Float> >(this->impl->bp)) {
        return Method::SecondOrderSimple;
    }

    if (std::holds_alternative<BackPropagator2Simd8<Float, MapType::STD_MAP> >(this->impl->bp)) {
        return Method::SecondOrderSimd8_stdmap;
    }

    if (std::holds_alternative<BackPropagator2Simd8<Float, MapType::STD_UNORDERED_MAP> >(this->impl->bp)) {
        return Method::SecondOrderSimd8_stdunorderedmap;
    }

    if (std::holds_alternative<BackPropagator2Simd8<Float, MapType::ANKERL_UNORDERED_DENSE> >(this->impl->bp)) {
        return Method::SecondOrderSimd8_ankerl;
    }

    if (std::holds_alternative<BackPropagator2Simd8<Float, MapType::BOOST_UNORDERED_MAP> >(this->impl->bp)) {
        return Method::SecondOrderSimd8_boost;
    }

    if (std::holds_alternative<BackPropagatorLossy<Float> >(this->impl->bp)) {
        return Method::FirstOrderLossy;
    }

    if (std::holds_alternative<BackPropagatorLossy<Float, true> >(this->impl->bp)) {
        return Method::FirstOrderVLossy;
    }

    if (std::holds_alternative<BackPropagatorLossyCompressed<Float> >(this->impl->bp)) {
        return Method::FirstOrderLossyCompressed;
    }

    if (std::holds_alternative<BackPropagatorVLossyCompressed<Float> >(this->impl->bp)) {
        return Method::FirstOrderVLossyCompressed;
    }

    throw std::runtime_error("Invalid backpropagator type");
}

template<class Float>
auto
Tape<Float>::get_order() const -> std::size_t
{
    if (std::holds_alternative<BackPropagator<Float> >(this->impl->bp)) {
        return 1;
    }

    if (std::holds_alternative<BackPropagator<Float, true> >(this->impl->bp)) {
        return 1;
    }

    if (std::holds_alternative<BackPropagator2<Float> >(this->impl->bp)) {
        return 2;
    }

    if (std::holds_alternative<BackPropagator2Simd8<Float, MapType::STD_MAP> >(this->impl->bp)) {
        return 2;
    }

    if (std::holds_alternative<BackPropagator2Simd8<Float, MapType::STD_UNORDERED_MAP> >(this->impl->bp)) {
        return 2;
    }

    if (std::holds_alternative<BackPropagator2Simd8<Float, MapType::ANKERL_UNORDERED_DENSE> >(this->impl->bp)) {
        return 2;
    }

    if (std::holds_alternative<BackPropagator2Simd8<Float, MapType::BOOST_UNORDERED_MAP> >(this->impl->bp)) {
        return 2;
    }

    if (std::holds_alternative<BackPropagatorLossy<Float> >(this->impl->bp)) {
        return 1;
    }

    if (std::holds_alternative<BackPropagatorLossy<Float, true> >(this->impl->bp)) {
        return 1;
    }

    if (std::holds_alternative<BackPropagatorLossyCompressed<Float> >(this->impl->bp)) {
        return 1;
    }

    if (std::holds_alternative<BackPropagatorVLossyCompressed<Float> >(this->impl->bp)) {
        return 1;
    }

    throw std::runtime_error("Invalid backpropagator type");
    return 0;
}

template<class Float>
void
Tape<Float>::set_checkpoint()
{
    std::visit([ops_size = data.ops.size()](auto& arg) { arg.set_checkpoint(ops_size); }, this->impl->bp);
}

template<class Float>
void
Tape<Float>::backpropagate()
{
    std::visit([&data = this->data](auto& arg) { arg.backpropagate_to(0, data); }, this->impl->bp);
}

template<class Float>
void
Tape<Float>::backpropagate_to(std::size_t to)
{
    std::visit([to, &data = this->data](auto& arg) { arg.backpropagate_to(to, data); }, this->impl->bp);
}

template<class Float>
template<bool ResetInPlace, bool Log>
void
Tape<Float>::backpropagate_and_reset_to(std::size_t to)
{
    std::visit(
      [to, &data = this->data](auto& arg) { arg.template backpropagate_to<true, ResetInPlace, Log>(to, data); },
      this->impl->bp);
}

template<class Float>
void
Tape<Float>::set_derivative(adhoc_type<Float> const& var, double deriv, std::size_t lane)
{
    if (var.is_active()) {
        std::visit([var_id = var.id, deriv, lane, ops_size = data.ops.size()](
                     auto& arg) { arg.set_derivative(var_id, deriv, ops_size, lane); },
                   this->impl->bp);
    }
}

template<class Float>
auto
Tape<Float>::get_derivative(adhoc_type<Float> const& var, std::size_t lane) const -> double
{
    if (var.is_active()) {
        return std::visit([var_id = var.id, lane](auto& arg) { return arg.get_derivative(var_id, lane); },
                          this->impl->bp);
    }
    return 0.;
}

template<class Float>
auto
Tape<Float>::get_derivative(adhoc_type<Float> const& var1, adhoc_type<Float> const& var2, std::size_t lane) const
  -> double
{
    if (var1.is_active() && var2.is_active()) {
        return std::visit([var_id1 = var1.id, var_id2 = var2.id, lane](
                            auto& arg) { return arg.get_derivative(var_id1, var_id2, lane); },
                          this->impl->bp);
    }
    return 0;
}

template<class Float>
void
Tape<Float>::zero_adjoints()
{
    std::visit([](auto& arg) { arg.zero_adjoints(); }, this->impl->bp);
}

template<class Float>
void
Tape<Float>::print() const
{
    std::cout << "Tape contains " << data.ops.size() << " operations:" << std::endl;

    std::size_t id_idx = 0;
    for (std::size_t i = 0; i < data.ops.size(); ++i) {
        const auto& op = data.ops[i];
        std::cout << "  [" << i << "] ";

        switch (op) {
            case OpCode::REG_INPUT:
                std::cout << "REG_INPUT: id " << data.ids[id_idx];
                id_idx += 1;
                break;
            case OpCode::REG_OUTPUT:
                std::cout << "REG_OUTPUT: id " << data.ids[id_idx + 1] << " = " << data.ids[id_idx];
                id_idx += 2;
                break;
            case OpCode::ADD:
                std::cout << "ADD: " << data.ids[id_idx + 2] << " = " << data.ids[id_idx] << " + "
                          << data.ids[id_idx + 1];
                id_idx += 3;
                break;
            case OpCode::SUB:
                std::cout << "SUB: " << data.ids[id_idx + 2] << " = " << data.ids[id_idx] << " - "
                          << data.ids[id_idx + 1];
                id_idx += 3;
                break;
            case OpCode::MUL:
                std::cout << "MUL: " << data.ids[id_idx + 2] << " = " << data.ids[id_idx] << " * "
                          << data.ids[id_idx + 1];
                id_idx += 3;
                break;
            case OpCode::ADD_C:
                std::cout << "ADD_C: " << data.ids[id_idx + 1] << " = " << data.ids[id_idx] << " + constant";
                id_idx += 2;
                break;
            case OpCode::SUB_C:
                std::cout << "SUB_C: " << data.ids[id_idx + 1] << " = constant - " << data.ids[id_idx];
                id_idx += 2;
                break;
            case OpCode::MUL_C:
                std::cout << "MUL_C: " << data.ids[id_idx + 1] << " = " << data.ids[id_idx] << " * constant";
                id_idx += 2;
                break;
            case OpCode::EXP:
                std::cout << "EXP: " << data.ids[id_idx + 1] << " = exp(" << data.ids[id_idx] << ")";
                id_idx += 2;
                break;
            case OpCode::LOG:
                std::cout << "LOG: " << data.ids[id_idx + 1] << " = log(" << data.ids[id_idx] << ")";
                id_idx += 2;
                break;
            case OpCode::ERF:
                std::cout << "ERF: " << data.ids[id_idx + 1] << " = erf(" << data.ids[id_idx] << ")";
                id_idx += 2;
                break;
            case OpCode::ERFC:
                std::cout << "ERFC: " << data.ids[id_idx + 1] << " = erfc(" << data.ids[id_idx] << ")";
                id_idx += 2;
                break;
            case OpCode::COS:
                std::cout << "COS: " << data.ids[id_idx + 1] << " = cos(" << data.ids[id_idx] << ")";
                id_idx += 2;
                break;
            case OpCode::NORM:
                std::cout << "NORM: " << data.ids[id_idx + 1] << " = norm(" << data.ids[id_idx] << ")";
                id_idx += 2;
                break;
            case OpCode::INV:
                std::cout << "INV: " << data.ids[id_idx + 1] << " = inv(" << data.ids[id_idx] << ")";
                id_idx += 2;
                break;
            case OpCode::ABS:
                std::cout << "ABS: " << data.ids[id_idx + 1] << " = abs(" << data.ids[id_idx] << ")";
                id_idx += 2;
                break;
            case OpCode::SQRT:
                std::cout << "SQRT: " << data.ids[id_idx + 1] << " = sqrt(" << data.ids[id_idx] << ")";
                id_idx += 2;
                break;
            case OpCode::POW_C:
                std::cout << "POW_C: " << data.ids[id_idx + 1] << " = pow(" << data.ids[id_idx] << ", constant)";
                id_idx += 2;
                break;
        }
        std::cout << std::endl;
    }
}

template<class Float>
auto
Tape<Float>::size_of(bool capacity) const -> std::size_t
{
    std::size_t size = 0;
    size += sizeof(*impl);
    size += sizeof(data.ops) + sizeof(OpCode) * (capacity ? data.ops.capacity() : data.ops.size());
    size += sizeof(data.ids) + sizeof(std::size_t) * (capacity ? data.ids.capacity() : data.ids.size());
    size += sizeof(data.vals) + sizeof(double) * (capacity ? data.vals.capacity() : data.vals.size());
    size += std::visit([capacity](auto& arg) { return arg.size_of(capacity); }, this->impl->bp);
    return size;
}

// no need to instantiate in header only mode
#ifndef ADHOC_HEADER_ONLY
template void
Tape<double>::backpropagate_and_reset_to<true, true>(std::size_t to);
template void
Tape<double>::backpropagate_and_reset_to<true, false>(std::size_t to);
template void
Tape<double>::backpropagate_and_reset_to<false, true>(std::size_t to);
template void
Tape<double>::backpropagate_and_reset_to<false, false>(std::size_t to);

template class Tape<double>;
#endif

} // namespace adhoc
