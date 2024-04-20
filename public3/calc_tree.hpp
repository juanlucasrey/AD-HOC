/*
 * AD-HOC, Automatic Differentiation for High Order Calculations
 *
 * This file is part of the AD-HOC distribution
 * (https://github.com/juanlucasrey/adhoc2).
 * Copyright (c) 2024 Juan Lucas Rey
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ADHOC3_CALC_TREE_HPP
#define ADHOC3_CALC_TREE_HPP

#include "calc_order.hpp"
#include "leafs.hpp"
#include "tuple_utils.hpp"

#include <tuple>

namespace adhoc3 {

template <class... Roots> class CalcTree {
  public:
    using leafs = decltype(detail::leafs_t<Roots...>::template call());
    using ValuesTuple = decltype(std::tuple_cat(
        leafs{},
        detail::calc_order_aux_t<false, true, Roots...>::template call()));

    using ValuesTupleInverse = decltype(std::tuple_cat(

        detail::calc_order_aux_t<false, false, Roots...>::template call(),
        leafs{}));

    std::array<double, std::tuple_size<ValuesTuple>{}> m_values{};

    inline void evaluate_fwd_aux(std::tuple<> /* nodes */) {}

    template <template <class> class Univariate, class Node,
              class... NodesToCalc>
    inline void
    evaluate_fwd_aux(std::tuple<Univariate<Node>, NodesToCalc...> /* nodes
    */) {
        using CurrentNode = Univariate<Node>;
        constexpr auto idx = get_idx_first2<CurrentNode>(ValuesTuple{});
        this->m_values[idx] = CurrentNode::v(this->val(Node{}));
        this->evaluate_fwd_aux(std::tuple<NodesToCalc...>{});
    }

    template <template <class, class> class Bivariate, class Node1, class Node2,
              class... NodesToCalc>
    inline void evaluate_fwd_aux(
        std::tuple<Bivariate<Node1, Node2>, NodesToCalc...> /* nodes */) {
        using CurrentNode = Bivariate<Node1, Node2>;
        constexpr auto idx = get_idx_first2<CurrentNode>(ValuesTuple{});
        this->m_values[idx] =
            CurrentNode::v(this->val(Node1{}), this->val(Node2{}));
        this->evaluate_fwd_aux(std::tuple<NodesToCalc...>{});
    }

    // why this special evaluation? because x/y is different from x * (1/y)
    template <class Node1, class Node2, class... NodesToCalc>
    inline void evaluate_fwd_aux(
        std::tuple<mul_t<Node1, inv_t<Node2>>, NodesToCalc...> /* nodes */) {
        using CurrentNode = mul_t<Node1, inv_t<Node2>>;
        constexpr auto idx = get_idx_first2<CurrentNode>(ValuesTuple{});
        this->m_values[idx] = this->val(Node1{}) / this->val(Node2{});
        this->evaluate_fwd_aux(std::tuple<NodesToCalc...>{});
    }

  public:
    explicit CalcTree() = default;
    explicit CalcTree(Roots... /* out */) {}

    template <constants::ArgType D>
    auto inline val(constants::CD<D> /* in */) const -> double {
        return constants::CD<D>::v();
    }

    template <class Input> auto inline val(Input /* in */) const -> double {
        constexpr auto idx = get_idx_first2<Input>(ValuesTuple{});
        return this->m_values[idx];
    }

    template <class Input> auto inline set(Input /* in */) -> double & {
        static_assert(has_type_tuple<Input, leafs>::value,
                      "only leafs are allowed to be set");
        constexpr auto idx = get_idx_first2<Input>(ValuesTuple{});
        return this->m_values[idx];
    }

    inline auto evaluate() {
        auto constexpr calcs =
            detail::calc_order_aux_t<false, true, Roots...>::template call();
        this->evaluate_fwd_aux(calcs);
    }
};

} // namespace adhoc3

#endif // ADHOC3_CALC_TREE_HPP
