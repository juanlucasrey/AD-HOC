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

#include "tape.hpp"

#include "adhoc/backpropagator1.hpp"
#include "adhoc/backpropagator1lossy.hpp"
#include "adhoc/backpropagator1lossycompressed.hpp"
#include "adhoc/backpropagator1lossycompressedpathreuse.hpp"
#include "adhoc/backpropagator1lossycompressedpathreusev.hpp"
#include "adhoc/backpropagator1lossypathreuse.hpp"
#include "adhoc/backpropagator2.hpp"
#include "adhoc/backpropagator2lossy.hpp"
#include "adhoc/position_impl.hpp"
#include "adhoc/tape_data.hpp"
#include "adhoc/vector_size_of.hpp"

#include <cstddef>
#include <iostream>
#include <variant>

namespace adhoc {

namespace {

template<class ContainerOps, class ContainerIds>
void
record_register(ContainerOps& ops, ContainerIds& ids, OpCode op, std::size_t arg_id)
{
    ops.push_back(op);
    ids.push_back(arg_id);
}

} // namespace

template<class Float, class TapeDataType>
Tape<Float, TapeDataType>::position_t::position_t()
  : impl(std::make_unique<PositionImpl>()){};

template<class Float, class TapeDataType>
Tape<Float, TapeDataType>::position_t::~position_t() = default;

template<class Float, class TapeDataType>
Tape<Float, TapeDataType>::position_t&
Tape<Float, TapeDataType>::position_t::operator=(position_t other)
{
    std::swap(this->impl, other.impl);
    return *this;
}

template<class Float, class TapeDataType>
Tape<Float, TapeDataType>::position_t::position_t(const position_t& other)
  : impl(std::make_unique<PositionImpl>(*other.impl))
{
}

template<class Float, class TapeDataType>
struct Tape<Float, TapeDataType>::Impl {
    std::variant<BackPropagator<Float>,
                 BackPropagator<Float, true>,
                 BackPropagator2<Float, MapType::ANKERL_UNORDERED_DENSE>,
                 BackPropagator2<Float, MapType::STD_MAP, true>,
                 BackPropagator2<Float, MapType::STD_UNORDERED_MAP, true>,
                 BackPropagator2<Float, MapType::ANKERL_UNORDERED_DENSE, true>,
                 BackPropagator2<Float, MapType::BOOST_UNORDERED_MAP, true>,
                 BackPropagatorLossy<Float>,
                 BackPropagatorLossy<Float, true>,
                 BackPropagatorLossyCompressed<Float>,
                 BackPropagatorLossyCompressed<Float, true>,
                 BackPropagatorLossyPathReuse<Float>,
                 BackPropagatorLossyPathReuse<Float, true>,
                 BackPropagatorLossyCompressedPathReuse<Float>,
                 BackPropagatorLossyCompressedPathReuse<Float, true>,
                 BackPropagatorLossyCompressedPathReuseV<Float>,
                 BackPropagatorLossyCompressedPathReuseV<Float, true>,
                 BackPropagator2Lossy<double, MapType::ANKERL_UNORDERED_DENSE>,
                 BackPropagator2Lossy<double, MapType::ANKERL_UNORDERED_DENSE, true> >
      bp = BackPropagator<Float>();
};

template<class Float, class TapeDataType>
Tape<Float, TapeDataType>::Tape(TapeDataType& tape_data)
  : impl(std::make_unique<Impl>())
  , data(tape_data)
{
}

template<class Float, class TapeDataType>
Tape<Float, TapeDataType>::~Tape() = default;

template<class Float, class TapeDataType>
void
Tape<Float, TapeDataType>::reserve_input(std::size_t count_registered)
{
    std::visit([count_registered](auto& arg) { arg.reserve_input(count_registered); }, this->impl->bp);
}

template<class Float, class TapeDataType>
void
Tape<Float, TapeDataType>::reserve_output(std::size_t count_registered)
{
    std::visit([count_registered](auto& arg) { arg.reserve_output(count_registered); }, this->impl->bp);
}

template<class Float, class TapeDataType>
void
Tape<Float, TapeDataType>::register_variable(adhoc_type<Float, TapeDataType> const& var)
{
    if (var.is_passive()) {
        std::size_t const new_id = this->data.generate_id();
        var.id = new_id;
        record_register(data.ops, data.ids, OpCode::REG_INPUT, new_id);
        std::visit([new_id](auto& arg) { arg.register_variable(new_id); }, this->impl->bp);
    }
}

template<class Float, class TapeDataType>
void
Tape<Float, TapeDataType>::register_variable(adhoc_type<Float, TapeDataType>& var)
{
    if (var.is_passive()) {
        std::size_t const new_id = this->data.generate_id();
        var.id = new_id;
        record_register(data.ops, data.ids, OpCode::REG_INPUT, new_id);
        std::visit([new_id](auto& arg) { arg.register_variable(new_id); }, this->impl->bp);
    }
}

template<class Float, class TapeDataType>
void
Tape<Float, TapeDataType>::register_output_variable(adhoc_type<Float, TapeDataType> const& var)
{
    if (var.is_active()) {
        std::size_t const new_id = this->data.generate_id();
        data.ids.push_back(var.id);
        var.id = new_id;
        record_register(data.ops, data.ids, OpCode::REG_OUTPUT, new_id);
        std::visit([new_id, ops_size = data.ops.size()](auto& arg) { arg.register_output_variable(new_id, ops_size); },
                   this->impl->bp);
    }
}

template<class Float, class TapeDataType>
void
Tape<Float, TapeDataType>::register_output_variable(adhoc_type<Float, TapeDataType>& var)
{
    if (var.is_active()) {
        std::size_t const new_id = this->data.generate_id();
        data.ids.push_back(var.id);
        var.id = new_id;
        record_register(data.ops, data.ids, OpCode::REG_OUTPUT, new_id);
        std::visit([new_id, ops_size = data.ops.size()](auto& arg) { arg.register_output_variable(new_id, ops_size); },
                   this->impl->bp);
    }
}

template<class Float, class TapeDataType>
void
Tape<Float, TapeDataType>::set_lanes(std::size_t num_lanes)
{
    std::visit([num_lanes](auto& arg) { arg.set_lanes(num_lanes); }, this->impl->bp);
}

template<class Float, class TapeDataType>
auto
Tape<Float, TapeDataType>::get_lanes() const -> std::size_t
{
    return std::visit([](auto& arg) { return arg.get_lanes(); }, this->impl->bp);
}

template<class Float, class TapeDataType>
void
Tape<Float, TapeDataType>::set_method(Method m)
{
    if (m == Method::FirstOrderSimple) {
        this->impl->bp.template emplace<BackPropagator<double> >();
    }
    else if (m == Method::FirstOrderSimd8) {
        this->impl->bp.template emplace<BackPropagator<double, true> >();
    }
    else if (m == Method::SecondOrderSimple) {
        this->impl->bp.template emplace<BackPropagator2<double, MapType::ANKERL_UNORDERED_DENSE> >();
    }
    else if (m == Method::SecondOrderSimd8_stdmap) {
        this->impl->bp.template emplace<BackPropagator2<double, MapType::STD_MAP, true> >();
    }
    else if (m == Method::SecondOrderSimd8_stdunorderedmap) {
        this->impl->bp.template emplace<BackPropagator2<double, MapType::STD_UNORDERED_MAP, true> >();
    }
    else if (m == Method::SecondOrderSimd8_ankerl) {
        this->impl->bp.template emplace<BackPropagator2<double, MapType::ANKERL_UNORDERED_DENSE, true> >();
    }
    else if (m == Method::SecondOrderSimd8_boost) {
        this->impl->bp.template emplace<BackPropagator2<double, MapType::BOOST_UNORDERED_MAP, true> >();
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
        this->impl->bp.template emplace<BackPropagatorLossyCompressed<double, true> >();
    }
    else if (m == Method::FirstOrderLossyPathReuse) {
        this->impl->bp.template emplace<BackPropagatorLossyPathReuse<double> >();
    }
    else if (m == Method::FirstOrderVLossyPathReuse) {
        this->impl->bp.template emplace<BackPropagatorLossyPathReuse<double, true> >();
    }
    else if (m == Method::FirstOrderLossyCompressedPathReuse) {
        this->impl->bp.template emplace<BackPropagatorLossyCompressedPathReuse<double> >();
    }
    else if (m == Method::FirstOrderVLossyCompressedPathReuse) {
        this->impl->bp.template emplace<BackPropagatorLossyCompressedPathReuse<double, true> >();
    }
    else if (m == Method::FirstOrderLossyCompressedPathReuseV) {
        this->impl->bp.template emplace<BackPropagatorLossyCompressedPathReuseV<double> >();
    }
    else if (m == Method::FirstOrderVLossyCompressedPathReuseV) {
        this->impl->bp.template emplace<BackPropagatorLossyCompressedPathReuseV<double, true> >();
    }
    else if (m == Method::SecondOrderLossy) {
        this->impl->bp.template emplace<BackPropagator2Lossy<double, MapType::ANKERL_UNORDERED_DENSE> >();
    }
    else if (m == Method::SecondOrderVLossy) {
        this->impl->bp.template emplace<BackPropagator2Lossy<double, MapType::ANKERL_UNORDERED_DENSE, true> >();
    }
}

template<class Float, class TapeDataType>
auto
Tape<Float, TapeDataType>::get_method() const -> Method
{
    if (std::holds_alternative<BackPropagator<Float> >(this->impl->bp)) {
        return Method::FirstOrderSimple;
    }

    if (std::holds_alternative<BackPropagator<Float, true> >(this->impl->bp)) {
        return Method::FirstOrderSimd8;
    }

    if (std::holds_alternative<BackPropagator2<Float, MapType::ANKERL_UNORDERED_DENSE> >(this->impl->bp)) {
        return Method::SecondOrderSimple;
    }

    if (std::holds_alternative<BackPropagator2<Float, MapType::STD_MAP, true> >(this->impl->bp)) {
        return Method::SecondOrderSimd8_stdmap;
    }

    if (std::holds_alternative<BackPropagator2<Float, MapType::STD_UNORDERED_MAP, true> >(this->impl->bp)) {
        return Method::SecondOrderSimd8_stdunorderedmap;
    }

    if (std::holds_alternative<BackPropagator2<Float, MapType::ANKERL_UNORDERED_DENSE, true> >(this->impl->bp)) {
        return Method::SecondOrderSimd8_ankerl;
    }

    if (std::holds_alternative<BackPropagator2<Float, MapType::BOOST_UNORDERED_MAP, true> >(this->impl->bp)) {
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

    if (std::holds_alternative<BackPropagatorLossyCompressed<Float, true> >(this->impl->bp)) {
        return Method::FirstOrderVLossyCompressed;
    }

    if (std::holds_alternative<BackPropagatorLossyPathReuse<Float> >(this->impl->bp)) {
        return Method::FirstOrderLossyPathReuse;
    }

    if (std::holds_alternative<BackPropagatorLossyPathReuse<Float, true> >(this->impl->bp)) {
        return Method::FirstOrderVLossyPathReuse;
    }

    if (std::holds_alternative<BackPropagatorLossyCompressedPathReuse<Float> >(this->impl->bp)) {
        return Method::FirstOrderLossyCompressedPathReuse;
    }

    if (std::holds_alternative<BackPropagatorLossyCompressedPathReuse<Float, true> >(this->impl->bp)) {
        return Method::FirstOrderVLossyCompressedPathReuse;
    }

    if (std::holds_alternative<BackPropagatorLossyCompressedPathReuseV<Float> >(this->impl->bp)) {
        return Method::FirstOrderLossyCompressedPathReuseV;
    }

    if (std::holds_alternative<BackPropagatorLossyCompressedPathReuseV<Float, true> >(this->impl->bp)) {
        return Method::FirstOrderVLossyCompressedPathReuseV;
    }

    if (std::holds_alternative<BackPropagator2Lossy<double, MapType::ANKERL_UNORDERED_DENSE> >(this->impl->bp)) {
        return Method::SecondOrderLossy;
    }

    if (std::holds_alternative<BackPropagator2Lossy<double, MapType::ANKERL_UNORDERED_DENSE, true> >(this->impl->bp)) {
        return Method::SecondOrderVLossy;
    }

    throw std::runtime_error("Invalid backpropagator type");
}

template<class Float, class TapeDataType>
auto
Tape<Float, TapeDataType>::get_order() const -> std::size_t
{
    if (std::holds_alternative<BackPropagator<Float> >(this->impl->bp)) {
        return 1;
    }

    if (std::holds_alternative<BackPropagator<Float, true> >(this->impl->bp)) {
        return 1;
    }

    if (std::holds_alternative<BackPropagator2<Float, MapType::ANKERL_UNORDERED_DENSE> >(this->impl->bp)) {
        return 2;
    }

    if (std::holds_alternative<BackPropagator2<Float, MapType::STD_MAP, true> >(this->impl->bp)) {
        return 2;
    }

    if (std::holds_alternative<BackPropagator2<Float, MapType::STD_UNORDERED_MAP, true> >(this->impl->bp)) {
        return 2;
    }

    if (std::holds_alternative<BackPropagator2<Float, MapType::ANKERL_UNORDERED_DENSE, true> >(this->impl->bp)) {
        return 2;
    }

    if (std::holds_alternative<BackPropagator2<Float, MapType::BOOST_UNORDERED_MAP, true> >(this->impl->bp)) {
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

    if (std::holds_alternative<BackPropagatorLossyCompressed<Float, true> >(this->impl->bp)) {
        return 1;
    }

    if (std::holds_alternative<BackPropagatorLossyPathReuse<Float> >(this->impl->bp)) {
        return 1;
    }

    if (std::holds_alternative<BackPropagatorLossyPathReuse<Float, true> >(this->impl->bp)) {
        return 1;
    }

    if (std::holds_alternative<BackPropagatorLossyCompressedPathReuse<Float> >(this->impl->bp)) {
        return 1;
    }

    if (std::holds_alternative<BackPropagatorLossyCompressedPathReuse<Float, true> >(this->impl->bp)) {
        return 1;
    }

    if (std::holds_alternative<BackPropagatorLossyCompressedPathReuseV<Float> >(this->impl->bp)) {
        return 1;
    }

    if (std::holds_alternative<BackPropagatorLossyCompressedPathReuseV<Float, true> >(this->impl->bp)) {
        return 1;
    }

    if (std::holds_alternative<BackPropagator2Lossy<double, MapType::ANKERL_UNORDERED_DENSE> >(this->impl->bp)) {
        return 2;
    }

    if (std::holds_alternative<BackPropagator2Lossy<double, MapType::ANKERL_UNORDERED_DENSE, true> >(this->impl->bp)) {
        return 2;
    }

    throw std::runtime_error("Invalid backpropagator type");
    return 0;
}

template<class Float, class TapeDataType>
void
Tape<Float, TapeDataType>::set_checkpoint()
{
    std::visit([ops_size = data.ops.size()](auto& arg) { arg.set_checkpoint(ops_size); }, this->impl->bp);
}

template<class Float, class TapeDataType>
void
Tape<Float, TapeDataType>::backpropagate()
{
    std::visit(
      [&data = this->data](auto& arg) {
          PositionImpl pos0;
          arg.template backpropagate_to<false, false, false>(pos0, data);
      },
      this->impl->bp);
}

template<class Float, class TapeDataType>
void
Tape<Float, TapeDataType>::backpropagate_to(position_t const& pos)
{
    std::visit(
      [pos, &data = this->data](auto& arg) { arg.template backpropagate_to<false, false, false>(*pos.impl, data); },
      this->impl->bp);
}

template<class Float, class TapeDataType>
template<bool ResetInPlace, bool Log>
void
Tape<Float, TapeDataType>::backpropagate_and_reset_to(position_t const& pos)
{
    std::visit(
      [pos, &data = this->data](auto& arg) { arg.template backpropagate_to<true, ResetInPlace, Log>(*pos.impl, data); },
      this->impl->bp);
}

template<class Float, class TapeDataType>
void
Tape<Float, TapeDataType>::set_derivative(adhoc_type<Float, TapeDataType> const& var, double deriv, std::size_t lane)
{
    if (var.is_active()) {
        std::visit([var_id = var.id, deriv, lane, ops_size = data.ops.size()](
                     auto& arg) { arg.set_derivative(var_id, deriv, ops_size, lane); },
                   this->impl->bp);
    }
}

template<class Float, class TapeDataType>
auto
Tape<Float, TapeDataType>::get_derivative(adhoc_type<Float, TapeDataType> const& var, std::size_t lane) const -> double
{
    if (var.is_active()) {
        return std::visit([var_id = var.id, lane](auto& arg) { return arg.get_derivative(var_id, lane); },
                          this->impl->bp);
    }
    return 0.;
}

template<class Float, class TapeDataType>
auto
Tape<Float, TapeDataType>::get_derivative(adhoc_type<Float, TapeDataType> const& var1,
                                          adhoc_type<Float, TapeDataType> const& var2,
                                          std::size_t lane) const -> double
{
    if (var1.is_active() && var2.is_active()) {
        return std::visit([var_id1 = var1.id, var_id2 = var2.id, lane](
                            auto& arg) { return arg.get_derivative(var_id1, var_id2, lane); },
                          this->impl->bp);
    }
    return 0;
}

template<class Float, class TapeDataType>
void
Tape<Float, TapeDataType>::zero_adjoints()
{
    std::visit([](auto& arg) { arg.zero_adjoints(); }, this->impl->bp);
}

template<class Float, class TapeDataType>
auto
Tape<Float, TapeDataType>::get_position() const -> position_t
{
    position_t result;
    result.impl = std::make_unique<PositionImpl>(data.ops.size(), data.ids.size(), data.vals.size());
    return result;
}

template<class Float, class TapeDataType>
void
Tape<Float, TapeDataType>::print() const
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

template<class Float, class TapeDataType>
auto
Tape<Float, TapeDataType>::size_of(bool capacity) const -> std::size_t
{
    std::size_t size = 0;
    size += sizeof(*impl);
    size += sizeof(data.ops) + vector_size_of(data.ops, capacity);
    size += sizeof(data.ids) + vector_size_of(data.ids, capacity);
    size += sizeof(data.vals) + vector_size_of(data.vals, capacity);
    size += std::visit([capacity](auto& arg) { return arg.size_of(capacity); }, this->impl->bp);
    return size;
}

// no need to instantiate in header only mode
#ifndef ADHOC_HEADER_ONLY
template void
Tape<double, TapeData<EnumVectorType::Simple, IdxVectorType::Simple> >::backpropagate_and_reset_to<true, true>(
  position_t const& to);
template void
Tape<double, TapeData<EnumVectorType::Simple, IdxVectorType::Simple> >::backpropagate_and_reset_to<true, false>(
  position_t const& to);
template void
Tape<double, TapeData<EnumVectorType::Simple, IdxVectorType::Simple> >::backpropagate_and_reset_to<false, true>(
  position_t const& to);
template void
Tape<double, TapeData<EnumVectorType::Simple, IdxVectorType::Simple> >::backpropagate_and_reset_to<false, false>(
  position_t const& to);

template class Tape<double, TapeData<EnumVectorType::Simple, IdxVectorType::Simple> >;

#ifndef _MSC_VER
template void
Tape<double, TapeData<EnumVectorType::BitCompression, IdxVectorType::Simple> >::backpropagate_and_reset_to<true, true>(
  position_t const& to);
template void
Tape<double, TapeData<EnumVectorType::BitCompression, IdxVectorType::Simple> >::backpropagate_and_reset_to<true, false>(
  position_t const& to);
template void
Tape<double, TapeData<EnumVectorType::BitCompression, IdxVectorType::Simple> >::backpropagate_and_reset_to<false, true>(
  position_t const& to);
template void
Tape<double, TapeData<EnumVectorType::BitCompression, IdxVectorType::Simple> >::backpropagate_and_reset_to<false,
                                                                                                           false>(
  position_t const& to);

template class Tape<double, TapeData<EnumVectorType::BitCompression, IdxVectorType::Simple> >;

template void
Tape<double, TapeData<EnumVectorType::Valuecompression, IdxVectorType::Simple> >::backpropagate_and_reset_to<true,
                                                                                                             true>(
  position_t const& to);
template void
Tape<double, TapeData<EnumVectorType::Valuecompression, IdxVectorType::Simple> >::backpropagate_and_reset_to<true,
                                                                                                             false>(
  position_t const& to);
template void
Tape<double, TapeData<EnumVectorType::Valuecompression, IdxVectorType::Simple> >::backpropagate_and_reset_to<false,
                                                                                                             true>(
  position_t const& to);
template void
Tape<double, TapeData<EnumVectorType::Valuecompression, IdxVectorType::Simple> >::backpropagate_and_reset_to<false,
                                                                                                             false>(
  position_t const& to);

template class Tape<double, TapeData<EnumVectorType::Valuecompression, IdxVectorType::Simple> >;

//
template void
Tape<double, TapeData<EnumVectorType::Simple, IdxVectorType::BitCompression> >::backpropagate_and_reset_to<true, true>(
  position_t const& to);
template void
Tape<double, TapeData<EnumVectorType::Simple, IdxVectorType::BitCompression> >::backpropagate_and_reset_to<true, false>(
  position_t const& to);
template void
Tape<double, TapeData<EnumVectorType::Simple, IdxVectorType::BitCompression> >::backpropagate_and_reset_to<false, true>(
  position_t const& to);
template void
Tape<double, TapeData<EnumVectorType::Simple, IdxVectorType::BitCompression> >::backpropagate_and_reset_to<false,
                                                                                                           false>(
  position_t const& to);

template class Tape<double, TapeData<EnumVectorType::Simple, IdxVectorType::BitCompression> >;

template void
Tape<double, TapeData<EnumVectorType::BitCompression, IdxVectorType::BitCompression> >::
  backpropagate_and_reset_to<true, true>(position_t const& to);
template void
Tape<double, TapeData<EnumVectorType::BitCompression, IdxVectorType::BitCompression> >::
  backpropagate_and_reset_to<true, false>(position_t const& to);
template void
Tape<double, TapeData<EnumVectorType::BitCompression, IdxVectorType::BitCompression> >::
  backpropagate_and_reset_to<false, true>(position_t const& to);
template void
Tape<double, TapeData<EnumVectorType::BitCompression, IdxVectorType::BitCompression> >::
  backpropagate_and_reset_to<false, false>(position_t const& to);

template class Tape<double, TapeData<EnumVectorType::BitCompression, IdxVectorType::BitCompression> >;
#endif

template void
Tape<double, TapeData<EnumVectorType::Valuecompression, IdxVectorType::BitCompression> >::
  backpropagate_and_reset_to<true, true>(position_t const& to);
template void
Tape<double, TapeData<EnumVectorType::Valuecompression, IdxVectorType::BitCompression> >::
  backpropagate_and_reset_to<true, false>(position_t const& to);
template void
Tape<double, TapeData<EnumVectorType::Valuecompression, IdxVectorType::BitCompression> >::
  backpropagate_and_reset_to<false, true>(position_t const& to);
template void
Tape<double, TapeData<EnumVectorType::Valuecompression, IdxVectorType::BitCompression> >::
  backpropagate_and_reset_to<false, false>(position_t const& to);

template class Tape<double, TapeData<EnumVectorType::Valuecompression, IdxVectorType::BitCompression> >;
#endif

} // namespace adhoc
