/*
 * AD-HOC, Automatic Differentiation for High Order Calculations
 *
 * This file is part of the AD-HOC distribution
 * (https://github.com/juanlucasrey/AD-HOC).
 * Copyright (c) 2024 Juan Lucas Rey
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

#ifndef ADHOC4_BACKPROPAGATOR_HPP
#define ADHOC4_BACKPROPAGATOR_HPP

#include "adhoc.hpp"
#include "backpropagator_aux.hpp"
#include "backpropagator_buffer_size.hpp"
#include "backpropagator_tools.hpp"
#include "calc_tree.hpp"
#include "dependency.hpp"
#include "differential_operator.hpp"
#include "sort.hpp"
#include "tape_buffer_size.hpp"
#include "utils/tuple.hpp"

#include <array>
#include <cstddef>
#include <tuple>
#include <utility>

namespace adhoc4 {

namespace detail {

template <class... Ids, std::size_t... Orders, class Node>
constexpr auto
create_single_tuple(std::tuple<der::d<Orders, Ids>...> diff_operator,
                    Node /* in */) {
    std::size_t constexpr pos = find<std::tuple<Ids...>, Node>();
    if constexpr (pos < sizeof...(Ids)) {
        return std::make_tuple(std::get<pos>(diff_operator));
    } else {
        return std::tuple<>{};
    }
};

template <class... Ids, std::size_t... Orders, class... CalculationNodes>
constexpr auto
order_differential_operator(std::tuple<der::d<Orders, Ids>...> diff_operator,
                            std::tuple<CalculationNodes...> nodes) {
    return std::apply(
        [diff_operator](auto &&...args) {
            return std::tuple_cat(create_single_tuple(diff_operator, args)...);
        },
        nodes);
}

} // namespace detail

enum class DerivativeType { TaylorCoefficient, Derivative };

template <class... InputsAndOutputsDers> class BackPropagator {
  private:
    using InputDers = decltype(std::tuple_cat(
        std::conditional_t<detail::is_derivative_input(InputsAndOutputsDers{}),
                           std::tuple<InputsAndOutputsDers>,
                           std::tuple<>>{}...));

    using OutputDers = decltype(std::tuple_cat(
        std::conditional_t<!detail::is_derivative_input(InputsAndOutputsDers{}),
                           std::tuple<InputsAndOutputsDers>,
                           std::tuple<>>{}...));

    std::array<double, sizeof...(InputsAndOutputsDers)> m_derivatives{};

  public:
    explicit BackPropagator() = delete;
    explicit BackPropagator(InputsAndOutputsDers... /* out */) {
        static_assert(std::tuple_size_v<InputDers> != 0,
                      "input derivatives empty!");
        static_assert(std::tuple_size_v<OutputDers> != 0,
                      "output derivatives empty!");
        // more checks to be included here later
    }

    template <DerivativeType DT = DerivativeType::TaylorCoefficient,
              class Input>
    auto inline get(Input in) const -> double {
        constexpr auto idx = find<std::tuple<InputsAndOutputsDers...>, Input>();
        static_assert(idx < sizeof...(InputsAndOutputsDers),
                      "derivative not found");

        if constexpr (is_one(in)) {
            return this->m_derivatives[idx];
        } else {
            return this->m_derivatives[idx] * multiplicity(in);
        }
    }

    template <class Input> auto inline set(Input /* in */) -> double & {
        constexpr auto idx = find<std::tuple<InputsAndOutputsDers...>, Input>();
        static_assert(idx < sizeof...(InputsAndOutputsDers),
                      "derivative not found");
        return this->m_derivatives[idx];
    }

    void inline reset() {
        std::fill(std::begin(this->m_derivatives),
                  std::end(this->m_derivatives), 0.);
    }

  public:
    template <class CalcTree> constexpr auto buffer_size() -> std::size_t {
        using PrimalNodesInverse = CalcTree::ValuesTupleInverse;
        constexpr auto primal_nodes_inverted = PrimalNodesInverse{};

        constexpr auto ordered_derivatives =
            std::tuple_cat(std::make_tuple(detail::order_differential_operator(
                InputsAndOutputsDers{}, primal_nodes_inverted))...);

        constexpr auto flags_inputs = std::apply(
            [](auto... der) {
                return std::integer_sequence<bool, detail::is_derivative_input(
                                                       der)...>{};
            },
            ordered_derivatives);

        constexpr auto ordered_derivatives_with_pos =
            detail::add_position(ordered_derivatives);

        constexpr auto inputs_and_outputs_with_pos =
            separate(ordered_derivatives_with_pos, flags_inputs);

        constexpr auto outputs_sorted_with_pos = sort_differential_operators(
            std::get<1>(inputs_and_outputs_with_pos), primal_nodes_inverted);

        constexpr auto inputs_and_outputs =
            separate(ordered_derivatives, flags_inputs);

        return detail::backpropagate_buffer_size<CalcTree>(
            outputs_sorted_with_pos, ordered_derivatives,
            std::get<0>(inputs_and_outputs));
    }

  public:
    template <class CalcTree> inline void backpropagate(CalcTree const &ct) {

        using PrimalNodesInverse = CalcTree::ValuesTupleInverse;
        constexpr auto primal_nodes_inverted = PrimalNodesInverse{};

        constexpr auto ordered_derivatives =
            std::tuple_cat(std::make_tuple(detail::order_differential_operator(
                InputsAndOutputsDers{}, primal_nodes_inverted))...);

        constexpr auto flags_inputs = std::apply(
            [](auto... der) {
                return std::integer_sequence<bool, detail::is_derivative_input(
                                                       der)...>{};
            },
            ordered_derivatives);

        constexpr auto ordered_derivatives_with_pos =
            detail::add_position(ordered_derivatives);

        constexpr auto inputs_and_outputs_with_pos =
            separate(ordered_derivatives_with_pos, flags_inputs);

        constexpr auto outputs_sorted_with_pos = sort_differential_operators(
            std::get<1>(inputs_and_outputs_with_pos), primal_nodes_inverted);

        constexpr auto inputs_and_outputs =
            separate(ordered_derivatives, flags_inputs);

        constexpr auto buffer_size =
            detail::backpropagate_buffer_size<CalcTree>(
                outputs_sorted_with_pos, ordered_derivatives,
                std::get<0>(inputs_and_outputs));

        std::array<double, buffer_size> buffer{};
        constexpr auto buffer_types =
            std::tuple_cat(std::array<detail::available_t, buffer_size>{});

        detail::backpropagate_aux(outputs_sorted_with_pos, ct,
                                  ordered_derivatives, this->m_derivatives,
                                  buffer_types, buffer,
                                  std::get<0>(inputs_and_outputs));
    }

    template <std::size_t BufferSize, class CalcTree>
    inline void backpropagate2(CalcTree const &ct) {

        using PrimalNodesInverse = CalcTree::ValuesTupleInverse;
        constexpr auto primal_nodes_inverted = PrimalNodesInverse{};

        constexpr auto ordered_derivatives =
            std::tuple_cat(std::make_tuple(detail::order_differential_operator(
                InputsAndOutputsDers{}, primal_nodes_inverted))...);

        constexpr auto flags_inputs = std::apply(
            [](auto... der) {
                return std::integer_sequence<bool, detail::is_derivative_input(
                                                       der)...>{};
            },
            ordered_derivatives);

        constexpr auto ordered_derivatives_with_pos =
            detail::add_position(ordered_derivatives);

        constexpr auto flags_outputs = std::apply(
            [](auto... der) {
                return std::integer_sequence<bool, !detail::is_derivative_input(
                                                       der)...>{};
            },
            ordered_derivatives);

        constexpr auto outputs_with_pos =
            filter(ordered_derivatives_with_pos, flags_outputs);

        constexpr auto outputs_sorted_with_pos = sort_differential_operators(
            outputs_with_pos, primal_nodes_inverted);

        constexpr auto inputs = filter(ordered_derivatives, flags_inputs);

        std::array<double, BufferSize> buffer{};
        constexpr auto buffer_types =
            std::tuple_cat(std::array<detail::available_t, BufferSize>{});

        detail::backpropagate_aux(outputs_sorted_with_pos, ct,
                                  ordered_derivatives, this->m_derivatives,
                                  buffer_types, buffer, inputs);
    }
};

} // namespace adhoc4

#endif // ADHOC4_BACKPROPAGATOR_HPP
