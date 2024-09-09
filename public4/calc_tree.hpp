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

#ifndef ADHOC4_CALC_TREE_HPP
#define ADHOC4_CALC_TREE_HPP

#include "adhoc.hpp"
#include "constants_type.hpp"
#include "utils/tuple.hpp"

#include <array>
#include <tuple>

namespace adhoc4 {

namespace detail {

template <constants::ArgType N>
constexpr auto get_dependents(constants::CD<N> /* in*/) {
    return std::tuple<>{};
}

template <detail::StringLiteral literal>
constexpr auto get_dependents(double_t<literal> /* in*/) {
    return std::tuple<>{};
}

template <template <class...> class Xvariate, class... Node>
constexpr auto get_dependents(Xvariate<Node...> /* in*/) {
    return std::tuple<Node...>{};
}

template <class Res>
constexpr auto calc_tree_t_aux(Res res, std::tuple<> /* nodes */) {
    return res;
}

template <class... Res, class Node, class... Nodes>
constexpr auto calc_tree_t_aux(std::tuple<Res...> res,
                               std::tuple<Node, Nodes...> /* nodes */) {

    constexpr bool other_types_depend_on_this =
        (depends(Nodes{}, Node{}) || ...);

    constexpr bool is_const = is_constant(Node{});

    if constexpr (other_types_depend_on_this || is_const) {
        // this_type will come up again because it is included on Nodes
        return calc_tree_t_aux(res, std::tuple<Nodes...>{});
    } else {
        constexpr auto deps = get_dependents(Node{});
        constexpr auto new_nodes = std::tuple_cat(deps, std::tuple<Nodes...>{});
        return calc_tree_t_aux(std::tuple<Node, Res...>{}, new_nodes);
    }
}

template <class... Outputs> constexpr auto calc_tree_t(Outputs... /* o */) {
    return calc_tree_t_aux(std::tuple<>{}, std::tuple<Outputs...>{});
}

} // namespace detail

template <class... Outputs> class CalcTree {
  private:
    using ValuesTuple = decltype(detail::calc_tree_t(Outputs{}...));
    std::array<double, std::tuple_size<ValuesTuple>{}> m_values{};

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
        constexpr auto idx = get_first_type_idx(ValuesTuple{}, CurrentNode{});
        this->m_values[idx] = CurrentNode::v(this->val(Inputs{})...);
        this->evaluate_aux(std::tuple<NodesToCalc...>{});
    }

    // why this special evaluation? because x/y is different from x * (1/y)
    template <class Node1, class Node2, class... NodesToCalc>
    inline void evaluate_aux(
        std::tuple<mul_t<Node1, inv_t<Node2>>, NodesToCalc...> /* nodes */) {
        using CurrentNode = mul_t<Node1, inv_t<Node2>>;
        constexpr auto idx = get_first_type_idx(ValuesTuple{}, CurrentNode{});
        this->m_values[idx] = this->val(Node1{}) / this->val(Node2{});
        this->evaluate_aux(std::tuple<NodesToCalc...>{});
    }

  public:
    explicit CalcTree() = delete;
    explicit CalcTree(Outputs... /* out */) {}

    template <constants::ArgType D>
    auto inline val(constants::CD<D> /* in */) const -> double {
        return constants::CD<D>::v();
    }

    template <class Input> auto inline val(Input in) const -> double {
        constexpr auto idx = get_first_type_idx(ValuesTuple{}, in);
        return this->m_values[idx];
    }

    template <class Input> auto inline set(Input in) -> double & {
        static_assert(is_input(in), "only inputs are allowed to be set");
        constexpr auto idx = get_first_type_idx(ValuesTuple{}, in);
        static_assert(idx < std::tuple_size<ValuesTuple>{}, "input not found");
        return this->m_values[idx];
    }

    inline void evaluate() { this->evaluate_aux(ValuesTuple{}); }
};

} // namespace adhoc4

#endif // ADHOC4_CALC_TREE_HPP
