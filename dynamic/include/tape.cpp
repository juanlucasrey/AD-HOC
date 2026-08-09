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
#include "adhoc/fwdpropagator1.hpp"
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

template<class type>
Tape<type>::position_t::position_t()
  : impl(std::make_unique<PositionImpl>()){};

template<class type>
Tape<type>::position_t::~position_t() = default;

template<class type>
Tape<type>::position_t&
Tape<type>::position_t::operator=(position_t other)
{
    std::swap(this->impl, other.impl);
    return *this;
}

template<class type>
Tape<type>::position_t::position_t(const position_t& other)
  : impl(std::make_unique<PositionImpl>(*other.impl))
{
}

template<class type>
struct Tape<type>::Impl {
    std::variant<BackPropagator<double>,
                 BackPropagator<double, true>,
                 BackPropagator2<double, MapType::ANKERL_UNORDERED_DENSE>,
                 BackPropagator2<double, MapType::ANKERL_UNORDERED_DENSE, true>,
                 BackPropagator2<double, MapType::STD_MAP>,
                 BackPropagator2<double, MapType::STD_MAP, true>,
                 BackPropagator2<double, MapType::STD_UNORDERED_MAP>,
                 BackPropagator2<double, MapType::STD_UNORDERED_MAP, true>,
                 BackPropagator2<double, MapType::BOOST_UNORDERED_MAP>,
                 BackPropagator2<double, MapType::BOOST_UNORDERED_MAP, true>,
                 BackPropagatorLossy<double>,
                 BackPropagatorLossy<double, true>,
                 BackPropagatorLossyCompressed<double>,
                 BackPropagatorLossyCompressed<double, true>,
                 BackPropagatorLossyPathReuse<double>,
                 BackPropagatorLossyPathReuse<double, true>,
                 BackPropagatorLossyCompressedPathReuse<double>,
                 BackPropagatorLossyCompressedPathReuse<double, true>,
                 BackPropagatorLossyCompressedPathReuseV<double>,
                 BackPropagatorLossyCompressedPathReuseV<double, true>,
                 BackPropagator2Lossy<double, MapType::ANKERL_UNORDERED_DENSE>,
                 BackPropagator2Lossy<double, MapType::ANKERL_UNORDERED_DENSE, true>,
                 FwdPropagator<double>,
                 FwdPropagator<double, true> >
      bp = BackPropagator<double>();

    std::size_t m_n_inputs{ 1 };
    std::size_t m_n_outputs{ 1 };
};

template<class type>
Tape<type>::Tape(tape_data_t& tape_data)
  : impl(std::make_unique<Impl>())
  , data(tape_data)
{
}

template<class type>
Tape<type>::~Tape() = default;

template<class type>
void
Tape<type>::reserve_input(std::size_t count_registered)
{
    std::visit([count_registered](auto& arg) { arg.reserve_input(count_registered); }, this->impl->bp);
}

template<class type>
void
Tape<type>::reserve_output(std::size_t count_registered)
{
    std::visit([count_registered](auto& arg) { arg.reserve_output(count_registered); }, this->impl->bp);
}

template<class type>
void
Tape<type>::register_variable(type const& var)
{
    if (var.is_passive()) {
        std::size_t const new_id = this->data.generate_id();
        var.id = new_id;
        record_register(data.ops, data.ids, OpCode::REG_INPUT, new_id);
        std::visit([new_id](auto& arg) { arg.register_variable(new_id); }, this->impl->bp);
    }
}

template<class type>
void
Tape<type>::register_variable(type& var)
{
    this->register_variable(static_cast<type const&>(var));
}

template<class type>
void
Tape<type>::register_output_variable(type const& var)
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

template<class type>
void
Tape<type>::register_output_variable(type& var)
{
    this->register_output_variable(static_cast<type const&>(var));
}

template<class type>
auto
Tape<type>::get_lanes() const -> std::size_t
{
    return std::visit([](auto& arg) { return arg.get_lanes(); }, this->impl->bp);
}

template<class type>
auto
Tape<type>::get_size() const -> std::size_t
{
    auto const method = this->get_method();
    bool const is_fwd = (method == Method::Fwd);
    return is_fwd ? this->impl->m_n_inputs : this->impl->m_n_outputs;
}

template<class type>
void
Tape<type>::configure(Method m, std::size_t n_inputs, std::size_t n_outputs, std::size_t num_lanes)
{
    this->impl->m_n_inputs = n_inputs;
    this->impl->m_n_outputs = n_outputs;
    bool const is_fwd = (m == Method::Fwd);
    std::size_t const max_lanes = is_fwd ? n_inputs : n_outputs;
    if (num_lanes == 0) {
        num_lanes = max_lanes;
    }

    num_lanes = std::min(num_lanes, max_lanes);

    if (m == Method::Fwd) {
        if (num_lanes == 1) {
            this->impl->bp.template emplace<FwdPropagator<double> >(n_inputs, n_outputs, num_lanes);
        }
        else {
            this->impl->bp.template emplace<FwdPropagator<double, true> >(n_inputs, n_outputs, num_lanes);
        }
    }
    else if (m == Method::Bwd) {
        if (num_lanes == 1) {
            this->impl->bp.template emplace<BackPropagator<double> >();
        }
        else {
            this->impl->bp.template emplace<BackPropagator<double, true> >();
            std::visit([num_lanes](auto& arg) { arg.set_lanes(num_lanes); }, this->impl->bp);
        }
    }
    else if (m == Method::BwdBuffer) {
        if (num_lanes == 1) {
            this->impl->bp.template emplace<BackPropagatorLossy<double> >();
        }
        else {
            this->impl->bp.template emplace<BackPropagatorLossy<double, true> >();
            std::visit([num_lanes](auto& arg) { arg.set_lanes(num_lanes); }, this->impl->bp);
        }
    }
    else if (m == Method::BwdBufferCompressed) {
        if (num_lanes == 1) {
            this->impl->bp.template emplace<BackPropagatorLossyCompressed<double> >();
        }
        else {
            this->impl->bp.template emplace<BackPropagatorLossyCompressed<double, true> >();
            std::visit([num_lanes](auto& arg) { arg.set_lanes(num_lanes); }, this->impl->bp);
        }
    }
    else if (m == Method::BwdBufferPathReuse) {
        if (num_lanes == 1) {
            this->impl->bp.template emplace<BackPropagatorLossyPathReuse<double> >();
        }
        else {
            this->impl->bp.template emplace<BackPropagatorLossyPathReuse<double, true> >();
            std::visit([num_lanes](auto& arg) { arg.set_lanes(num_lanes); }, this->impl->bp);
        }
    }
    else if (m == Method::BwdBufferCompressedPathReuse) {
        if (num_lanes == 1) {
            this->impl->bp.template emplace<BackPropagatorLossyCompressedPathReuse<double> >();
        }
        else {
            this->impl->bp.template emplace<BackPropagatorLossyCompressedPathReuse<double, true> >();
            std::visit([num_lanes](auto& arg) { arg.set_lanes(num_lanes); }, this->impl->bp);
        }
    }
    else if (m == Method::BwdBufferCompressedPathReuseV) {
        if (num_lanes == 1) {
            this->impl->bp.template emplace<BackPropagatorLossyCompressedPathReuseV<double> >();
        }
        else {
            this->impl->bp.template emplace<BackPropagatorLossyCompressedPathReuseV<double, true> >();
            std::visit([num_lanes](auto& arg) { arg.set_lanes(num_lanes); }, this->impl->bp);
        }
    }

    else if (m == Method::SecondOrderSimple || m == Method::SecondOrderSimd8_ankerl) {
        if (num_lanes == 1) {
            this->impl->bp.template emplace<BackPropagator2<double, MapType::ANKERL_UNORDERED_DENSE> >();
        }
        else {
            this->impl->bp.template emplace<BackPropagator2<double, MapType::ANKERL_UNORDERED_DENSE, true> >();
            // only allows 8 lanes for now, but could be extended to allow more
            std::visit([](auto& arg) { arg.set_lanes(8); }, this->impl->bp);
        }
    }
    else if (m == Method::SecondOrderSimd8_stdmap) {
        if (num_lanes == 1) {
            this->impl->bp.template emplace<BackPropagator2<double, MapType::STD_MAP> >();
        }
        else {
            this->impl->bp.template emplace<BackPropagator2<double, MapType::STD_MAP, true> >();
            // only allows 8 lanes for now, but could be extended to allow more
            std::visit([](auto& arg) { arg.set_lanes(8); }, this->impl->bp);
        }
    }
    else if (m == Method::SecondOrderSimd8_stdunorderedmap) {
        if (num_lanes == 1) {
            this->impl->bp.template emplace<BackPropagator2<double, MapType::STD_UNORDERED_MAP> >();
        }
        else {
            this->impl->bp.template emplace<BackPropagator2<double, MapType::STD_UNORDERED_MAP, true> >();
            // only allows 8 lanes for now, but could be extended to allow more
            std::visit([](auto& arg) { arg.set_lanes(8); }, this->impl->bp);
        }
    }
    else if (m == Method::SecondOrderSimd8_boost) {
        if (num_lanes == 1) {
            this->impl->bp.template emplace<BackPropagator2<double, MapType::BOOST_UNORDERED_MAP> >();
        }
        else {
            this->impl->bp.template emplace<BackPropagator2<double, MapType::BOOST_UNORDERED_MAP, true> >();
            // only allows 8 lanes for now, but could be extended to allow more
            std::visit([](auto& arg) { arg.set_lanes(8); }, this->impl->bp);
        }
    }
    else if (m == Method::SecondOrderLossy || m == Method::SecondOrderVLossy) {
        if (num_lanes == 1) {
            this->impl->bp.template emplace<BackPropagator2Lossy<double, MapType::ANKERL_UNORDERED_DENSE> >();
        }
        else {
            this->impl->bp.template emplace<BackPropagator2Lossy<double, MapType::ANKERL_UNORDERED_DENSE, true> >();
            std::visit([num_lanes](auto& arg) { arg.set_lanes(num_lanes); }, this->impl->bp);
        }
    }
}

template<class type>
auto
Tape<type>::get_method() const -> Method
{
    if (std::holds_alternative<BackPropagator<double> >(this->impl->bp)) {
        return Method::Bwd;
    }

    if (std::holds_alternative<BackPropagator<double, true> >(this->impl->bp)) {
        return Method::Bwd;
    }

    if (std::holds_alternative<BackPropagatorLossy<double> >(this->impl->bp)) {
        return Method::BwdBuffer;
    }

    if (std::holds_alternative<BackPropagatorLossy<double, true> >(this->impl->bp)) {
        return Method::BwdBuffer;
    }

    if (std::holds_alternative<BackPropagatorLossyCompressed<double> >(this->impl->bp)) {
        return Method::BwdBufferCompressed;
    }

    if (std::holds_alternative<BackPropagatorLossyCompressed<double, true> >(this->impl->bp)) {
        return Method::BwdBufferCompressed;
    }

    if (std::holds_alternative<BackPropagatorLossyPathReuse<double> >(this->impl->bp)) {
        return Method::BwdBufferPathReuse;
    }

    if (std::holds_alternative<BackPropagatorLossyPathReuse<double, true> >(this->impl->bp)) {
        return Method::BwdBufferPathReuse;
    }

    if (std::holds_alternative<BackPropagatorLossyCompressedPathReuse<double> >(this->impl->bp)) {
        return Method::BwdBufferCompressedPathReuse;
    }

    if (std::holds_alternative<BackPropagatorLossyCompressedPathReuse<double, true> >(this->impl->bp)) {
        return Method::BwdBufferCompressedPathReuse;
    }

    if (std::holds_alternative<BackPropagatorLossyCompressedPathReuseV<double> >(this->impl->bp)) {
        return Method::BwdBufferCompressedPathReuseV;
    }

    if (std::holds_alternative<BackPropagatorLossyCompressedPathReuseV<double, true> >(this->impl->bp)) {
        return Method::BwdBufferCompressedPathReuseV;
    }

    if (std::holds_alternative<BackPropagator2<double, MapType::ANKERL_UNORDERED_DENSE> >(this->impl->bp)) {
        return Method::SecondOrderSimple;
    }

    if (std::holds_alternative<BackPropagator2<double, MapType::STD_MAP, true> >(this->impl->bp)) {
        return Method::SecondOrderSimd8_stdmap;
    }

    if (std::holds_alternative<BackPropagator2<double, MapType::STD_UNORDERED_MAP, true> >(this->impl->bp)) {
        return Method::SecondOrderSimd8_stdunorderedmap;
    }

    if (std::holds_alternative<BackPropagator2<double, MapType::ANKERL_UNORDERED_DENSE, true> >(this->impl->bp)) {
        return Method::SecondOrderSimd8_ankerl;
    }

    if (std::holds_alternative<BackPropagator2<double, MapType::BOOST_UNORDERED_MAP, true> >(this->impl->bp)) {
        return Method::SecondOrderSimd8_boost;
    }

    if (std::holds_alternative<BackPropagator2Lossy<double, MapType::ANKERL_UNORDERED_DENSE> >(this->impl->bp)) {
        return Method::SecondOrderLossy;
    }

    if (std::holds_alternative<BackPropagator2Lossy<double, MapType::ANKERL_UNORDERED_DENSE, true> >(this->impl->bp)) {
        return Method::SecondOrderVLossy;
    }

    if (std::holds_alternative<FwdPropagator<double> >(this->impl->bp)) {
        return Method::Fwd;
    }

    if (std::holds_alternative<FwdPropagator<double, true> >(this->impl->bp)) {
        return Method::Fwd;
    }

    throw std::runtime_error("Invalid backpropagator type");
}

template<class type>
auto
Tape<type>::get_order() const -> std::size_t
{
    if (std::holds_alternative<BackPropagator<double> >(this->impl->bp)) {
        return 1;
    }

    if (std::holds_alternative<BackPropagator<double, true> >(this->impl->bp)) {
        return 1;
    }

    if (std::holds_alternative<BackPropagator2<double, MapType::ANKERL_UNORDERED_DENSE> >(this->impl->bp)) {
        return 2;
    }

    if (std::holds_alternative<BackPropagator2<double, MapType::STD_MAP, true> >(this->impl->bp)) {
        return 2;
    }

    if (std::holds_alternative<BackPropagator2<double, MapType::STD_UNORDERED_MAP, true> >(this->impl->bp)) {
        return 2;
    }

    if (std::holds_alternative<BackPropagator2<double, MapType::ANKERL_UNORDERED_DENSE, true> >(this->impl->bp)) {
        return 2;
    }

    if (std::holds_alternative<BackPropagator2<double, MapType::BOOST_UNORDERED_MAP, true> >(this->impl->bp)) {
        return 2;
    }

    if (std::holds_alternative<BackPropagatorLossy<double> >(this->impl->bp)) {
        return 1;
    }

    if (std::holds_alternative<BackPropagatorLossy<double, true> >(this->impl->bp)) {
        return 1;
    }

    if (std::holds_alternative<BackPropagatorLossyCompressed<double> >(this->impl->bp)) {
        return 1;
    }

    if (std::holds_alternative<BackPropagatorLossyCompressed<double, true> >(this->impl->bp)) {
        return 1;
    }

    if (std::holds_alternative<BackPropagatorLossyPathReuse<double> >(this->impl->bp)) {
        return 1;
    }

    if (std::holds_alternative<BackPropagatorLossyPathReuse<double, true> >(this->impl->bp)) {
        return 1;
    }

    if (std::holds_alternative<BackPropagatorLossyCompressedPathReuse<double> >(this->impl->bp)) {
        return 1;
    }

    if (std::holds_alternative<BackPropagatorLossyCompressedPathReuse<double, true> >(this->impl->bp)) {
        return 1;
    }

    if (std::holds_alternative<BackPropagatorLossyCompressedPathReuseV<double> >(this->impl->bp)) {
        return 1;
    }

    if (std::holds_alternative<BackPropagatorLossyCompressedPathReuseV<double, true> >(this->impl->bp)) {
        return 1;
    }

    if (std::holds_alternative<BackPropagator2Lossy<double, MapType::ANKERL_UNORDERED_DENSE> >(this->impl->bp)) {
        return 2;
    }

    if (std::holds_alternative<BackPropagator2Lossy<double, MapType::ANKERL_UNORDERED_DENSE, true> >(this->impl->bp)) {
        return 2;
    }

    if (std::holds_alternative<FwdPropagator<double> >(this->impl->bp)) {
        return 1;
    }

    if (std::holds_alternative<FwdPropagator<double, true> >(this->impl->bp)) {
        return 1;
    }

    throw std::runtime_error("Invalid backpropagator type");
    return 0;
}

template<class type>
void
Tape<type>::set_checkpoint()
{
    std::visit([ops_size = data.ops.size()](auto& arg) { arg.set_checkpoint(ops_size); }, this->impl->bp);
}

template<class type>
void
Tape<type>::backpropagate()
{
    position_t pos0;
    this->backpropagate_to(pos0);
}

template<class type>
void
Tape<type>::backpropagate_to(position_t const& pos)
{
    std::visit([pos, &data = this->data](auto& arg) { arg.template backpropagate_to<false>(*pos.impl, data); },
               this->impl->bp);
}

template<class type>
void
Tape<type>::backpropagate_and_reset_to(position_t const& pos)
{
    std::visit([pos, &data = this->data](auto& arg) { arg.template backpropagate_to<true>(*pos.impl, data); },
               this->impl->bp);

    this->reset_to(pos);
}

template<class type>
void
Tape<type>::reset_to(position_t const& pos)
{
    std::visit([pos, &data = this->data](auto& arg) { arg.reset(*pos.impl); }, this->impl->bp);
    this->data.reset(pos.impl->op_position, pos.impl->val_position, pos.impl->id_position);
}

template<class type>
void
Tape<type>::set_derivative(type const& var, double deriv, std::size_t lane)
{
    if (var.is_active()) {
        std::visit([var_id = var.id, deriv, lane, ops_size = data.ops.size()](
                     auto& arg) { arg.set_derivative(var_id, deriv, ops_size, lane); },
                   this->impl->bp);
    }
}

template<class type>
auto
Tape<type>::get_derivative(type const& var, std::size_t lane) const -> double
{
    if (var.is_active()) {
        return std::visit([var_id = var.id, lane](auto& arg) { return arg.get_derivative(var_id, lane); },
                          this->impl->bp);
    }
    return 0.;
}

template<class type>
auto
Tape<type>::get_derivative(std::size_t idx_input, std::size_t idx_output) const -> double
{
    if (std::holds_alternative<FwdPropagator<double> >(this->impl->bp)) {
        auto const& arg = std::get<FwdPropagator<double> >(this->impl->bp);
        return arg.get_derivative2(idx_input, idx_output);
    }

    if (std::holds_alternative<FwdPropagator<double, true> >(this->impl->bp)) {
        auto const& arg = std::get<FwdPropagator<double, true> >(this->impl->bp);
        return arg.get_derivative2(idx_input, idx_output);
    }

    return 0;
}

template<class type>
auto
Tape<type>::get_derivative(type const& var1, type const& var2, std::size_t lane) const -> double
{
    if (var1.is_active() && var2.is_active()) {
        return std::visit([var_id1 = var1.id, var_id2 = var2.id, lane](
                            auto& arg) { return arg.get_derivative(var_id1, var_id2, lane); },
                          this->impl->bp);
    }
    return 0;
}

template<class type>
void
Tape<type>::zero_adjoints()
{
    std::visit([](auto& arg) { arg.zero_adjoints(); }, this->impl->bp);
}

template<class type>
auto
Tape<type>::get_position() const -> position_t
{
    position_t result;
    result.impl = std::make_unique<PositionImpl>(data.ops.size(), data.ids.size(), data.vals.size());
    return result;
}

template<class type>
void
Tape<type>::print() const
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

template<class type>
auto
Tape<type>::size_of(bool capacity) const -> std::size_t
{
    std::size_t size = 0;
    size += sizeof(*impl);
    size += sizeof(data.ops) + vector_size_of(data.ops, capacity);
    size += sizeof(data.ids) + vector_size_of(data.ids, capacity);
    size += sizeof(data.vals) + vector_size_of(data.vals, capacity);
    size += std::visit([capacity](auto& arg) { return arg.size_of(capacity); }, this->impl->bp);
    return size;
}

template<class type>
Tape<type>::subrange_t::subrange_t(std::size_t size, std::size_t lanes)
  : m_size(size)
  , m_lanes(lanes)
{
}

template<class type>
auto
Tape<type>::subrange_t::begin() const -> range_t
{
    return range_t{ *this, std::false_type{} };
}

template<class type>
auto
Tape<type>::subrange_t::end() const -> range_t
{
    return range_t{ *this, std::true_type{} };
}

template<class type>
Tape<type>::subrange_t::range_t::range_t() = default;

template<class type>
template<bool IsEnd>
Tape<type>::subrange_t::range_t::range_t(subrange_t const& subrange, std::bool_constant<IsEnd> /* isend */)
  : m_size(subrange.m_size)
  , m_lanes(subrange.m_lanes)
  , m_global_index(IsEnd ? subrange.m_size : 0)
{
}

template<class type>
auto
Tape<type>::subrange_t::range_t::operator++() -> range_t&
{
    this->m_global_index = std::min(this->m_global_index + this->m_lanes, this->m_size);
    return *this;
}

template<class type>
auto
Tape<type>::subrange_t::range_t::operator++(int) -> range_t
{
    range_t temp = *this;
    ++(*this);
    return temp;
}

template<class type>
auto
Tape<type>::subrange_t::range_t::operator==(range_t const& rhs) const -> bool
{
    // m_size and m_lanes never change
    return rhs.m_global_index ==
           this->m_global_index /* && rhs.m_size == this->m_size && rhs.m_lanes == this->m_lanes */;
}

template<class type>
auto
Tape<type>::subrange_t::range_t::operator!=(range_t const& rhs) const -> bool
{
    return !(this->operator==(rhs));
}

template<class type>
auto
Tape<type>::subrange_t::range_t::operator*() const -> range_info_t
{
    return range_info_t{ this->m_size, this->m_lanes, this->m_global_index };
}

template<class type>
auto
Tape<type>::subrange_t::range_t::range_info_t::begin() const -> inner_range_t
{
    return inner_range_t{ *this, std::false_type{} };
}

template<class type>
auto
Tape<type>::subrange_t::range_t::range_info_t::end() const -> inner_range_t
{
    return inner_range_t{ *this, std::true_type{} };
}

template<class type>
Tape<type>::subrange_t::range_t::range_info_t::inner_range_t::inner_range_t() = default;

template<class type>
template<bool IsEnd>
Tape<type>::subrange_t::range_t::range_info_t::inner_range_t::inner_range_t(range_info_t const& range,
                                                                            std::bool_constant<IsEnd> /* isend */)
  : m_sub_index(IsEnd ? std::min(range.global_index + range.m_lanes, range.m_size) - range.global_index : 0)
  , m_global_index(range.global_index + m_sub_index)
{
}

template<class type>
auto
Tape<type>::subrange_t::range_t::range_info_t::inner_range_t::operator++() -> inner_range_t&
{
    ++m_sub_index;
    ++m_global_index;
    return *this;
}

template<class type>
auto
Tape<type>::subrange_t::range_t::range_info_t::inner_range_t::operator++(int) -> inner_range_t
{
    inner_range_t temp = *this;
    ++(*this);
    return temp;
}

template<class type>
auto
Tape<type>::subrange_t::range_t::range_info_t::inner_range_t::operator==(inner_range_t const& rhs) const -> bool
{
    return rhs.m_sub_index == this->m_sub_index && rhs.m_global_index == this->m_global_index;
}

template<class type>
auto
Tape<type>::subrange_t::range_t::range_info_t::inner_range_t::operator!=(inner_range_t const& rhs) const -> bool
{
    return !(this->operator==(rhs));
}

template<class type>
auto
Tape<type>::subrange_t::range_t::range_info_t::inner_range_t::operator*() const -> inner_range_info_t
{
    return inner_range_info_t{ this->m_sub_index, this->m_global_index };
}

// no need to instantiate in header only mode
#ifndef ADHOC_HEADER_ONLY
template class Tape<adhoc_type<double, "main", TapeData<double, EnumVectorType::Simple, IdxVectorType::Simple> > >;
template class Tape<adhoc_type<double, "driver", TapeData<double, EnumVectorType::Simple, IdxVectorType::Simple> > >;
#endif

} // namespace adhoc
