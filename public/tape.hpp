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

#ifndef ADHOC_TAPE_HPP
#define ADHOC_TAPE_HPP

#include "calc_order.hpp"
#include "leafs.hpp"
#include "tape_size_bwd.hpp"
// #include "tape_size_fwd.hpp"
#include "tuple_utils.hpp"
#include "utils.hpp"

#include <array>
#include <iostream>
#include <tuple>

namespace adhoc {

template <class... ActiveLeafsAndRootsDerivatives> class Tape {

    using leafs = decltype(detail::leafs_t<
                           ActiveLeafsAndRootsDerivatives...>::template call());
    using ValuesTuple = decltype(std::tuple_cat(
        leafs{},
        detail::calc_order_aux_t<
            true, ActiveLeafsAndRootsDerivatives...>::template call()));

    std::array<double, std::tuple_size<ValuesTuple>{}> m_values{};
    std::array<double, sizeof...(ActiveLeafsAndRootsDerivatives)>
        m_derivatives{};

    using buffersize = std::integral_constant<
        std::size_t, tape_size_bwd<ActiveLeafsAndRootsDerivatives...>()>;
    std::array<double, buffersize::value> m_buffer{};

  public:
    explicit Tape() = default;
    explicit Tape(ActiveLeafsAndRootsDerivatives... /* out */) {}
    template <class Derived> auto inline val(Base<Derived> var) const -> double;
    template <class Derived> auto inline set(Base<Derived> var) -> double &;
    template <class Derived> auto inline der(Base<Derived> var) const -> double;
    template <class Derived> auto inline der(Base<Derived> var) -> double &;
    void inline reset_der() {
        std::fill(std::begin(this->m_derivatives),
                  std::end(this->m_derivatives), 0.);
    }

  private:
    template <class Derived>
    auto inline der_int(Base<Derived> var) const -> double;
    template <class Derived> auto inline der_int(Base<Derived> var) -> double &;

    template <constants::ArgType D>
    inline auto val_aux(constants::CD<D> /* node */) -> double {
        return constants::CD<D>::v();
    }

    template <class Node> inline auto val_aux(Node /* node */) -> double {
        constexpr auto idx = get_idx_first2<Node>(ValuesTuple{});
        return this->m_values[idx];
    }

    inline void evaluate_fwd_aux(std::tuple<> /* nodes */) {}

    template <template <class> class Univariate, class Node,
              class... NodesToCalc>
    inline void
    evaluate_fwd_aux(std::tuple<Univariate<Node>, NodesToCalc...> /* nodes */) {
        using CurrentNode = Univariate<Node>;
        constexpr auto idx = get_idx_first2<CurrentNode>(ValuesTuple{});
        this->m_values[idx] = CurrentNode::v(this->val_aux(Node{}));
        this->evaluate_fwd_aux(std::tuple<NodesToCalc...>{});
    }

    template <template <class, class> class Bivariate, class Node1, class Node2,
              class... NodesToCalc>
    inline void evaluate_fwd_aux(
        std::tuple<Bivariate<Node1, Node2>, NodesToCalc...> /* nodes */) {
        using CurrentNode = Bivariate<Node1, Node2>;
        constexpr auto idx = get_idx_first2<CurrentNode>(ValuesTuple{});
        this->m_values[idx] =
            CurrentNode::v(this->val_aux(Node1{}), this->val_aux(Node2{}));
        this->evaluate_fwd_aux(std::tuple<NodesToCalc...>{});
    }

  public:
    inline auto evaluate_fwd() {
        auto constexpr calcs = detail::calc_order_aux_t<
            true, ActiveLeafsAndRootsDerivatives...>::template call();
        this->evaluate_fwd_aux(calcs);
    }

  private:
    class available_t {};

    template <class CurrentNode, class NextNode, class... NodesAlive,
              class... NodesToCalc>
    inline void univariate_bwd(const double &current_node_d,
                               std::tuple<NodesAlive...> /* deriv_ids */,
                               std::tuple<NodesToCalc...> /* nodes */) {
        constexpr bool is_next_node_leaf =
            has_type<NextNode, ActiveLeafsAndRootsDerivatives...>();
        constexpr bool is_current_node_root =
            has_type<CurrentNode, ActiveLeafsAndRootsDerivatives...>();

        constexpr auto position_root =
            get_idx_first2<CurrentNode>(std::tuple<NodesAlive...>{});

        using NodesAlive2 = decltype(replace_first2<position_root, available_t>(
            std::tuple<NodesAlive...>{}));

        constexpr bool next_node_needs_new_storing =
            !is_next_node_leaf && !has_type<NextNode, NodesAlive...>();

        constexpr auto position_next_available =
            get_idx_first2<available_t>(NodesAlive2{});

        using NodesAlive3 = std::conditional_t<
            next_node_needs_new_storing,
            decltype(replace_first2<position_next_available, NextNode>(
                NodesAlive2{})),
            NodesAlive2>;

        constexpr auto position_next = get_idx_first2<NextNode>(NodesAlive3{});

        double &next_node_total_d = is_next_node_leaf
                                        ? this->der_int(NextNode{})
                                        : this->m_buffer[position_next];

        double const next_node_d =
            current_node_d * (is_current_node_root
                                  ? this->der_int(CurrentNode{})
                                  : this->m_buffer[position_root]);
        if constexpr (next_node_needs_new_storing) {
            next_node_total_d = next_node_d;
        } else {
            next_node_total_d += next_node_d;
        }

        this->evaluate_bwd_aux(NodesAlive3{}, std::tuple<NodesToCalc...>{});
    }

    template <class... NodesAlive>
    inline void
    evaluate_bwd_aux(std::tuple<NodesAlive...> /* nodes_intermediate */,
                     std::tuple<> /* nodes */) {}

    template <class... NodesAlive, template <class, class> class Bivariate,
              class NextNode1, class NextNode2, class... NodesToCalc>
    inline void evaluate_bwd_aux(std::tuple<NodesAlive...> /* deriv_ids */,
                                 std::tuple<Bivariate<NextNode1, NextNode2>,
                                            NodesToCalc...> /* nodes */) {
        using CurrentNode = Bivariate<NextNode1, NextNode2>;

        constexpr bool node1_has_0_deriv =
            !equal_or_depends_many<NextNode1,
                                   ActiveLeafsAndRootsDerivatives...>();
        constexpr bool node2_has_0_deriv =
            !equal_or_depends_many<NextNode2,
                                   ActiveLeafsAndRootsDerivatives...>();
        static_assert(!node1_has_0_deriv || !node2_has_0_deriv);

        // ideally std::is_same_v<Node1, Node2> should be differentiated as type
        if constexpr (std::is_same_v<NextNode1, NextNode2> ||
                      node1_has_0_deriv || node2_has_0_deriv) {
            double const current_node_d =
                std::is_same_v<NextNode1, NextNode2>
                    ? 2 * CurrentNode::d1 /* could be d1 or d2 */ (
                              this->val_aux(CurrentNode{}),
                              this->val_aux(NextNode1{}),
                              this->val_aux(NextNode2{}))
                : node1_has_0_deriv
                    ? CurrentNode::d2(this->val_aux(CurrentNode{}),
                                      this->val_aux(NextNode1{}),
                                      this->val_aux(NextNode2{}))
                    : CurrentNode::d1(this->val_aux(CurrentNode{}),
                                      this->val_aux(NextNode1{}),
                                      this->val_aux(NextNode2{}));

            using NextNodeAlive =
                std::conditional_t<node1_has_0_deriv, NextNode2, NextNode1>;

            this->univariate_bwd<CurrentNode, NextNodeAlive>(
                current_node_d, std::tuple<NodesAlive...>{},
                std::tuple<NodesToCalc...>{});
        } else {
            constexpr bool is_next_node1_leaf =
                has_type<NextNode1, ActiveLeafsAndRootsDerivatives...>();
            constexpr bool is_next_node2_leaf =
                has_type<NextNode2, ActiveLeafsAndRootsDerivatives...>();
            constexpr bool is_current_node_root =
                has_type<CurrentNode, ActiveLeafsAndRootsDerivatives...>();

            constexpr auto position_root =
                get_idx_first2<CurrentNode>(std::tuple<NodesAlive...>{});

            using NodesAlive2 =
                decltype(replace_first2<position_root, available_t>(
                    std::tuple<NodesAlive...>{}));

            constexpr bool next_node1_needs_new_storing =
                !is_next_node1_leaf && !has_type<NextNode1, NodesAlive...>();

            constexpr auto position_next_node1_available =
                get_idx_first2<available_t>(NodesAlive2{});

            using NodesAlive3 = std::conditional_t<
                next_node1_needs_new_storing,
                decltype(replace_first2<position_next_node1_available,
                                        NextNode1>(NodesAlive2{})),
                NodesAlive2>;

            constexpr bool next_node2_needs_new_storing =
                !is_next_node2_leaf && !has_type<NextNode2, NodesAlive...>();

            constexpr auto position_next_node2_available =
                get_idx_first2<available_t>(NodesAlive3{});

            using NodesAlive4 = std::conditional_t<
                next_node2_needs_new_storing,
                decltype(replace_first2<position_next_node2_available,
                                        NextNode2>(NodesAlive3{})),
                NodesAlive3>;

            constexpr auto position_next1 =
                get_idx_first2<NextNode1>(NodesAlive4{});
            constexpr auto position_next2 =
                get_idx_first2<NextNode2>(NodesAlive4{});

            double &next_node1_val = is_next_node1_leaf
                                         ? this->der_int(NextNode1{})
                                         : this->m_buffer[position_next1];

            double &next_node2_val = is_next_node2_leaf
                                         ? this->der_int(NextNode2{})
                                         : this->m_buffer[position_next2];

            double const current_node_d = is_current_node_root
                                              ? this->der_int(CurrentNode{})
                                              : this->m_buffer[position_root];

            auto d = CurrentNode::d(this->val_aux(CurrentNode{}),
                                    this->val_aux(NextNode1{}),
                                    this->val_aux(NextNode2{}));
            d[0] *= current_node_d;
            d[1] *= current_node_d;

            if constexpr (next_node1_needs_new_storing) {
                next_node1_val = d[0];
            } else {
                next_node1_val += d[0];
            }

            if constexpr (next_node2_needs_new_storing) {
                next_node2_val = d[1];
            } else {
                next_node2_val += d[1];
            }

            this->evaluate_bwd_aux(NodesAlive4{}, std::tuple<NodesToCalc...>{});
        }
    }

    template <class... NodesAlive, template <class> class Univariate,
              class Node, class... NodesToCalc>
    inline void
    evaluate_bwd_aux(std::tuple<NodesAlive...> /* deriv_ids */,
                     std::tuple<Univariate<Node>, NodesToCalc...> /* nodes */) {
        using CurrentNode = Univariate<Node>;

        double const d =
            CurrentNode::d(this->val_aux(CurrentNode{}), this->val_aux(Node{}));

        this->univariate_bwd<CurrentNode, Node>(d, std::tuple<NodesAlive...>{},
                                                std::tuple<NodesToCalc...>{});
    }

  public:
    inline void evaluate_bwd() {
        constexpr auto nodes_ordered =
            detail::calc_order_t<false, ActiveLeafsAndRootsDerivatives...>();
        this->evaluate_bwd_aux(
            typename generate_tuple_type<available_t,
                                         buffersize::value>::type{},
            nodes_ordered);
    }
};

template <class... ActiveLeafsAndRootsDerivatives>
template <class Derived>
auto Tape<ActiveLeafsAndRootsDerivatives...>::val(Base<Derived> /* in */) const
    -> double {
    constexpr auto idx = get_idx_first2<Derived>(ValuesTuple{});
    return this->m_values[idx];
}

template <class... ActiveLeafsAndRootsDerivatives>
template <class Derived>
auto Tape<ActiveLeafsAndRootsDerivatives...>::set(Base<Derived> /* in */)
    -> double & {
    static_assert(get_idx_first2<Derived>(leafs{}) != std::tuple_size<leafs>{},
                  "only leafs are allowed to be set");
    constexpr auto idx = get_idx_first2<Derived>(ValuesTuple{});
    return this->m_values[idx];
}

template <class... ActiveLeafsAndRootsDerivatives>
template <class Derived>
auto Tape<ActiveLeafsAndRootsDerivatives...>::der(Base<Derived> /* in */) const
    -> double {
    static_assert(has_type<Derived, ActiveLeafsAndRootsDerivatives...>(),
                  "Only nodes or leafs can have derivatives set or read");
    constexpr auto idx = get_idx_first2<Derived>(
        std::tuple<ActiveLeafsAndRootsDerivatives...>{});
    return this->m_derivatives[idx];
}

template <class... ActiveLeafsAndRootsDerivatives>
template <class Derived>
auto Tape<ActiveLeafsAndRootsDerivatives...>::der(Base<Derived> /* in */)
    -> double & {
    static_assert(has_type<Derived, ActiveLeafsAndRootsDerivatives...>(),
                  "Only nodes or leafs can have derivatives set or read");
    constexpr auto idx = get_idx_first2<Derived>(
        std::tuple<ActiveLeafsAndRootsDerivatives...>{});
    return this->m_derivatives[idx];
}

template <class... ActiveLeafsAndRootsDerivatives>
template <class Derived>
auto Tape<ActiveLeafsAndRootsDerivatives...>::der_int(
    Base<Derived> /* in */) const -> double {
    constexpr auto idx = get_idx_first2<Derived>(
        std::tuple<ActiveLeafsAndRootsDerivatives...>{});
    return this->m_derivatives[idx];
}

template <class... ActiveLeafsAndRootsDerivatives>
template <class Derived>
auto Tape<ActiveLeafsAndRootsDerivatives...>::der_int(Base<Derived> /* in */)
    -> double & {
    constexpr auto idx = get_idx_first2<Derived>(
        std::tuple<ActiveLeafsAndRootsDerivatives...>{});
    return this->m_derivatives[idx];
}

} // namespace adhoc

#endif // ADHOC_TAPE_HPP
