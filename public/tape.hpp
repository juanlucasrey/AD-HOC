#ifndef ADHOC_TAPE_HPP
#define ADHOC_TAPE_HPP

#include "calc_order.hpp"
#include "leafs.hpp"
#include "tape_size_bwd.hpp"
#include "tape_size_fwd.hpp"
#include "tuple_utils.hpp"
#include "utils.hpp"

#include <array>
#include <iostream>
#include <tuple>

namespace adhoc {

namespace detail {

template <class ValuesTuple, class ActiveLeafsAndRootsDerivativesTuple>
class Tape2;

template <class... Values, class... ActiveLeafsAndRootsDerivatives>
class Tape2<std::tuple<Values...>,
            std::tuple<ActiveLeafsAndRootsDerivatives...>> {
    std::array<double, sizeof...(Values)> m_values{};
    std::array<double, sizeof...(ActiveLeafsAndRootsDerivatives)>
        m_derivatives{};
    using buffersize = std::integral_constant<
        std::size_t, tape_size_bwd<ActiveLeafsAndRootsDerivatives...>()>;
    std::array<double, buffersize::value> m_buffer{};

  public:
    template <class Derived> auto inline val(Base<Derived> var) const -> double;
    template <class Derived> auto inline val(Base<Derived> var) -> double &;
    template <class Derived> auto inline der(Base<Derived> var) const -> double;
    template <class Derived> auto inline der(Base<Derived> var) -> double &;

    template <class Derived>
    auto inline val_check(Base<Derived> var) -> double &;
    template <class Derived>
    auto inline der_check(Base<Derived> var) -> double &;

  private:
    template <constants::ArgType D>
    inline auto val_aux(constants::CD<D> /* node */) -> double {
        return constants::CD<D>::v();
    }

    template <class Node> inline auto val_aux(Node /* node */) -> double {
        return this->val(Node{});
    }

    inline void evaluate_fwd_aux(std::tuple<> /* nodes */) {}

    template <template <class> class Univariate, class Node,
              class... NodesToCalc>
    inline void
    evaluate_fwd_aux(std::tuple<Univariate<Node>, NodesToCalc...> /* nodes */) {
        using CurrentNode = Univariate<Node>;
        this->val(CurrentNode{}) = CurrentNode::v(this->val_aux(Node{}));
        this->evaluate_fwd_aux(std::tuple<NodesToCalc...>{});
    }

    template <template <class, class> class Bivariate, class Node1, class Node2,
              class... NodesToCalc>
    inline void evaluate_fwd_aux(
        std::tuple<Bivariate<Node1, Node2>, NodesToCalc...> /* nodes */) {
        using CurrentNode = Bivariate<Node1, Node2>;
        this->val(CurrentNode{}) =
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
                                        ? this->der(NextNode{})
                                        : this->m_buffer[position_next];

        double const next_node_d =
            current_node_d * (is_current_node_root
                                  ? this->der(CurrentNode{})
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
                                         ? this->der(NextNode1{})
                                         : this->m_buffer[position_next1];

            double &next_node2_val = is_next_node2_leaf
                                         ? this->der(NextNode2{})
                                         : this->m_buffer[position_next2];

            double const current_node_d = is_current_node_root
                                              ? this->der(CurrentNode{})
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

template <class... Values, class... ActiveLeafsAndRootsDerivatives>
template <class Derived>
auto Tape2<std::tuple<Values...>,
           std::tuple<ActiveLeafsAndRootsDerivatives...>>::
    val(Base<Derived> /* in */) const -> double {
    if constexpr (has_type<Derived, Values...>()) {
        constexpr auto idx = idx_type<Derived, Values...>();
        return this->m_values[idx];
    } else {
        return 0;
    }
}

template <class... Values, class... ActiveLeafsAndRootsDerivatives>
template <class Derived>
auto Tape2<
    std::tuple<Values...>,
    std::tuple<ActiveLeafsAndRootsDerivatives...>>::val(Base<Derived> /* in */)
    -> double & {
    if constexpr (has_type<Derived, Values...>()) {
        constexpr auto idx = idx_type<Derived, Values...>();
        return this->m_values[idx];
    } else {
        return this->m_buffer[0];
    }
}

template <class... Values, class... ActiveLeafsAndRootsDerivatives>
template <class Derived>
auto Tape2<std::tuple<Values...>,
           std::tuple<ActiveLeafsAndRootsDerivatives...>>::
    der(Base<Derived> /* in */) const -> double {
    if constexpr (has_type<Derived, ActiveLeafsAndRootsDerivatives...>()) {
        constexpr auto idx =
            idx_type<Derived, ActiveLeafsAndRootsDerivatives...>();
        return this->m_derivatives[idx];
    } else {
        return 0;
    }
}

template <class... Values, class... ActiveLeafsAndRootsDerivatives>
template <class Derived>
auto Tape2<
    std::tuple<Values...>,
    std::tuple<ActiveLeafsAndRootsDerivatives...>>::der(Base<Derived> /* in */)
    -> double & {
    if constexpr (has_type<Derived, ActiveLeafsAndRootsDerivatives...>()) {
        constexpr auto idx =
            idx_type<Derived, ActiveLeafsAndRootsDerivatives...>();
        return this->m_derivatives[idx];
    } else {
        return this->m_buffer[0];
    }
}

template <class... Values, class... ActiveLeafsAndRootsDerivatives>
template <class Derived>
auto Tape2<std::tuple<Values...>,
           std::tuple<ActiveLeafsAndRootsDerivatives...>>::
    val_check(Base<Derived> /* in */) -> double & {
    static_assert(has_type<Derived, Values...>());
    constexpr auto idx = idx_type<Derived, Values...>();
    return this->m_values[idx];
}

template <class... Values, class... ActiveLeafsAndRootsDerivatives>
template <class Derived>
auto Tape2<std::tuple<Values...>,
           std::tuple<ActiveLeafsAndRootsDerivatives...>>::
    der_check(Base<Derived> /* in */) -> double & {
    static_assert(has_type<Derived, ActiveLeafsAndRootsDerivatives...>(),
                  "variable not active");
    constexpr auto idx = idx_type<Derived, ActiveLeafsAndRootsDerivatives...>();
    return this->m_derivatives[idx];
}

} // namespace detail

template <class... ActiveLeafsAndRoots>
auto inline Tape(ActiveLeafsAndRoots... /* out */) {
    constexpr auto leafs =
        detail::leafs_t<ActiveLeafsAndRoots...>::template call();
    constexpr auto nodes =
        detail::calc_order_aux_t<true, ActiveLeafsAndRoots...>::template call();
    return detail::Tape2<decltype(std::tuple_cat(leafs, nodes)),
                         std::tuple<ActiveLeafsAndRoots...>>();
}

} // namespace adhoc

#endif // ADHOC_TAPE_HPP