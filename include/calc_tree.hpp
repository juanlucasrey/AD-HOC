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

#ifndef CALC_TREE_HPP
#define CALC_TREE_HPP

#include "adhoc.hpp"
#include "constants_type.hpp"
#include "dependency.hpp"
#include "utils/tuple.hpp"

#include <array>
#include <tuple>

namespace adhoc {

namespace detail {

template <constants::ArgType N>
constexpr auto get_dependents(constants::CD<N> /* in*/) {
    return std::tuple<>{};
}

template <StringLiteral literal>
constexpr auto get_dependents(double_t<literal> /* in*/) {
    return std::tuple<>{};
}

template <template <class...> class Xvariate, class... Node>
constexpr auto get_dependents(Xvariate<Node...> /* in*/) {
    return std::tuple<Node...>{};
}

template <class Result, class ResultInputs>
constexpr auto calc_tree_t_aux(Result result, ResultInputs result_inputs,
                               std::tuple<> /* nodes */) {
    // we make sure inputs are first in the calculation tree
    return std::tuple_cat(result_inputs, result);
}

template <class... NonInputNodes, class... InputNodes, class Node,
          class... Nodes>
constexpr auto calc_tree_t_aux(std::tuple<NonInputNodes...> res1,
                               std::tuple<InputNodes...> res2,
                               std::tuple<Node, Nodes...> /* nodes */) {

    constexpr bool other_nodes_depend_on_this_node =
        (depends(Nodes{}, Node{}) || ...);

    constexpr bool is_const = is_constant(Node{});

    if constexpr (other_nodes_depend_on_this_node || is_const) {
        // this node will come up again because it is included on the Nodes
        return calc_tree_t_aux(res1, res2, std::tuple<Nodes...>{});
    } else {
        constexpr auto deps = get_dependents(Node{});
        constexpr auto new_nodes = std::tuple_cat(deps, std::tuple<Nodes...>{});
        if constexpr (is_input(Node{})) {
            // we add it to the input nodes
            return calc_tree_t_aux(std::tuple<NonInputNodes...>{},
                                   std::tuple<Node, InputNodes...>{},
                                   new_nodes);
        } else {
            // we add it to the non-input nodes
            return calc_tree_t_aux(std::tuple<Node, NonInputNodes...>{},
                                   std::tuple<InputNodes...>{}, new_nodes);
        }
    }
}

template <class... Outputs> constexpr auto calc_tree_t(Outputs... /* o */) {
    return calc_tree_t_aux(std::tuple<>{}, std::tuple<>{},
                           std::tuple<Outputs...>{});
}

} // namespace detail

template <class... Outputs> class CalcTree {
  public:
    using ValuesTuple = decltype(detail::calc_tree_t(Outputs{}...));
    using ValuesTupleInverse = decltype(invert(ValuesTuple{}));

  private:
    std::array<double, std::tuple_size_v<ValuesTuple>> m_values{};

    inline void evaluate_aux(std::tuple<> /* nodes */) {}

    template <detail::StringLiteral literal, class... NodesToCalc>
    inline void
    evaluate_aux(std::tuple<double_t<literal>, NodesToCalc...> /* nodes */) {
        // this node should be an input and filled by user
        this->evaluate_aux(std::tuple<NodesToCalc...>{});
    }

    template <template <class...> class Xvariate, class... Inputs,
              class... NodesToCalc>
    inline void
    evaluate_aux(std::tuple<Xvariate<Inputs...>, NodesToCalc...> /* nodes */) {
        using CurrentNode = Xvariate<Inputs...>;
        constexpr auto idx = find<ValuesTuple, CurrentNode>();
        this->m_values[idx] = CurrentNode::v(this->get(Inputs{})...);
        this->evaluate_aux(std::tuple<NodesToCalc...>{});
    }

    // why this special evaluation? because x/y is different from x * (1/y)
    template <class Node1, class Node2, class... NodesToCalc>
    inline void evaluate_aux(
        std::tuple<mul_t<Node1, inv_t<Node2>>, NodesToCalc...> /* nodes */) {
        using CurrentNode = mul_t<Node1, inv_t<Node2>>;
        constexpr auto idx = find<ValuesTuple, CurrentNode>();
        this->m_values[idx] = this->get(Node1{}) / this->get(Node2{});
        this->evaluate_aux(std::tuple<NodesToCalc...>{});
    }

  public:
    explicit CalcTree() = delete;
    explicit CalcTree(Outputs... /* out */) {}

    template <constants::ArgType D>
    auto inline get(constants::CD<D> /* in */) const -> double {
        return constants::CD<D>::v();
    }

    template <class Input> auto inline get(Input /* in */) const -> double {
        constexpr auto idx = find<ValuesTuple, Input>();
        return this->m_values[idx];
    }

    template <class Input> auto inline set(Input in) -> double & {
        static_assert(is_input(in), "only inputs are allowed to be set");
        constexpr auto idx = find<ValuesTuple, Input>();
        static_assert(idx < std::tuple_size_v<ValuesTuple>, "input not found");
        return this->m_values[idx];
    }

    inline void evaluate() { this->evaluate_aux(ValuesTuple{}); }
};

} // namespace adhoc

#endif // CALC_TREE_HPP
