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
#include <tuple>

#include "../tests4/type_name.hpp"
#include <iostream>

namespace adhoc4 {

namespace detail {

template <class I> constexpr auto is_derivative_input(I /* in */) -> bool {
    return false;
}

template <StringLiteral... literals, std::size_t... Orders>
constexpr auto
is_derivative_input(std::tuple<der::d<Orders, double_t<literals>>...> /* in */)
    -> bool {
    return true;
}

template <StringLiteral... literals, std::size_t... Orders>
constexpr auto is_derivative_input(
    std::tuple<const der::d<Orders, double_t<literals>>...> /* in */) -> bool {
    return true;
}

// template <class... Ids, std::size_t... Orders, std::size_t... Powers>
// constexpr auto order(std::tuple<der::d<Orders, Ids>...> /* id1 */) {
//     return (Orders + ...);
// }

// template <class... Ops>
// constexpr auto max_orders(std::tuple<Ops...> /* id1 */) {
//     constexpr std::array<std::size_t, sizeof...(Ops)> temp{order(Ops{})...};
//     auto result = std::max_element(temp.begin(), temp.end());
//     return *result;
// }

template <class... Ids, std::size_t... Orders, class Node>
constexpr auto
create_single_tuple(std::tuple<der::d<Orders, Ids>...> diff_operator, Node in) {

    std::size_t constexpr pos = get_first_type_idx(std::tuple<Ids...>{}, in);
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

template <class... Ders>
constexpr auto select_root_derivatives(std::tuple<Ders...> /* ids */) {
    return std::tuple_cat(
        std::conditional_t<is_derivative_input(Ders{}), std::tuple<Ders>,
                           std::tuple<>>{}...);
}

template <class... Ders>
constexpr auto select_non_root_derivatives(std::tuple<Ders...> /* ids */) {
    return std::tuple_cat(
        std::conditional_t<!is_derivative_input(Ders{}), std::tuple<Ders>,
                           std::tuple<>>{}...);
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
        constexpr auto idx =
            get_first_type_idx(std::tuple<InputsAndOutputsDers...>{}, in);
        static_assert(idx < sizeof...(InputsAndOutputsDers),
                      "derivative not found");
        return this->m_derivatives[idx];
    }

    template <class Input> auto inline set(Input in) -> double & {
        constexpr auto idx =
            get_first_type_idx(std::tuple<InputsAndOutputsDers...>{}, in);
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

        constexpr auto inputs =
            detail::select_root_derivatives(ordered_derivatives);
        constexpr auto outputs =
            detail::select_non_root_derivatives(ordered_derivatives);

        constexpr auto outputs_sorted =
            sort_differential_operators(outputs, primal_nodes_inverted);

        constexpr auto node_derivative_location = std::tuple_cat(
            std::array<detail::on_interface_t, size(outputs_sorted)>{});

        return detail::backpropagate_buffer_size<CalcTree>(
            node_derivative_location, outputs_sorted, ordered_derivatives,
            inputs);
    }

  public:
    template <class CalcTree> inline void backpropagate(CalcTree const &ct) {

        using PrimalNodesInverse = CalcTree::ValuesTupleInverse;
        constexpr auto primal_nodes_inverted = PrimalNodesInverse{};

        constexpr auto ordered_derivatives =
            std::tuple_cat(std::make_tuple(detail::order_differential_operator(
                InputsAndOutputsDers{}, primal_nodes_inverted))...);

        constexpr auto inputs =
            detail::select_root_derivatives(ordered_derivatives);
        constexpr auto outputs =
            detail::select_non_root_derivatives(ordered_derivatives);

        constexpr auto outputs_sorted =
            sort_differential_operators(outputs, primal_nodes_inverted);

        constexpr auto node_derivative_location = std::tuple_cat(
            std::array<detail::on_interface_t, size(outputs_sorted)>{});
        // constexpr auto nodes_derivative =
        //     expand_tree(nodes_value, ordered_roots,
        //     output_derivatives_ordered);

        // constexpr auto buffer_size = 30;

        constexpr auto buffer_size =
            detail::backpropagate_buffer_size<CalcTree>(
                node_derivative_location, outputs_sorted, ordered_derivatives,
                inputs);

        // constexpr auto buffer_size = 40;
        // constexpr auto buffer_size = tape_buffer_size(
        //     output_derivatives_ordered, ordered_inputs, primal_nodes);

        // std::cout << buffer_size << std::endl;

        std::array<double, buffer_size> buffer{};
        constexpr auto buffer_types =
            std::tuple_cat(std::array<detail::available_t, buffer_size>{});

        detail::backpropagate_aux(node_derivative_location, outputs_sorted, ct,
                                  ordered_derivatives, this->m_derivatives,
                                  buffer_types, buffer, inputs);
    }
};

} // namespace adhoc4

#endif // ADHOC4_BACKPROPAGATOR_HPP
