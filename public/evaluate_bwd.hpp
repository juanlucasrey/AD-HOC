#ifndef ADHOC_EVALUATE_BWD_HPP
#define ADHOC_EVALUATE_BWD_HPP

#include "calc_order.hpp"
#include "constants_type.hpp"
#include "tape.hpp"
#include "tape_size_bwd.hpp"
#include "tuple_utils.hpp"

#include <array>

namespace adhoc {

namespace detail {

class available_t {};

template <class... Values, class... ActiveLeafsAndRootsDerivatives,
          class... NodesAlive>
inline void evaluate_bwd(
    detail::Tape2<std::tuple<Values...>,
                  std::tuple<ActiveLeafsAndRootsDerivatives...>> & /* in */,
    std::tuple<NodesAlive...> /* nodes_intermediate */,
    std::tuple<> /* nodes */) {}

template <class... Values, class... ActiveLeafsAndRootsDerivatives,
          class... NodesAlive, template <class> class Univariate, class Node,
          class... NodesToCalc>
inline void
evaluate_bwd(detail::Tape2<std::tuple<Values...>,
                           std::tuple<ActiveLeafsAndRootsDerivatives...>> &in,
             std::tuple<NodesAlive...> /* deriv_ids */,
             std::tuple<Univariate<Node>, NodesToCalc...> /* nodes */);

template <class... Values, class... ActiveLeafsAndRootsDerivatives,
          class... NodesAlive, template <class, class> class Bivariate,
          class NextNode1, class NextNode2, class... NodesToCalc>
inline void evaluate_bwd(
    detail::Tape2<std::tuple<Values...>,
                  std::tuple<ActiveLeafsAndRootsDerivatives...>> &in,
    std::tuple<NodesAlive...> /* deriv_ids */,
    std::tuple<Bivariate<NextNode1, NextNode2>, NodesToCalc...> /* nodes*/);

template <class CurrentNode, class NextNode, class... Values,
          class... ActiveLeafsAndRootsDerivatives, class... NodesAlive,
          class... NodesToCalc>
inline void univariate_bwd(
    Tape2<std::tuple<Values...>, std::tuple<ActiveLeafsAndRootsDerivatives...>>
        &in,
    const double &current_node_d, std::tuple<NodesAlive...> /* deriv_ids */,
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

    using NodesAlive4 =
        std::conditional_t<next_node_needs_new_storing,
                           decltype(replace_first2<position_next_available,
                                                   NextNode>(NodesAlive2{})),
                           NodesAlive2>;

    constexpr auto position_next = get_idx_first2<NextNode>(NodesAlive4{});

    double &next_node_total_d =
        is_next_node_leaf ? in.der(NextNode{}) : in.m_buffer[position_next];

    double const next_node_d =
        current_node_d * (is_current_node_root ? in.der(CurrentNode{})
                                               : in.m_buffer[position_root]);
    if constexpr (next_node_needs_new_storing) {
        next_node_total_d = next_node_d;
    } else {
        next_node_total_d += next_node_d;
    }

    evaluate_bwd(in, NodesAlive4{}, std::tuple<NodesToCalc...>{});
}

template <class... Values, class... ActiveLeafsAndRootsDerivatives,
          class... NodesAlive, template <class, class> class Bivariate,
          class NextNode1, class NextNode2, class... NodesToCalc>
inline void evaluate_bwd(
    detail::Tape2<std::tuple<Values...>,
                  std::tuple<ActiveLeafsAndRootsDerivatives...>> &in,
    std::tuple<NodesAlive...> /* deriv_ids */,
    std::tuple<Bivariate<NextNode1, NextNode2>, NodesToCalc...> /* nodes */) {
    using CurrentNode = Bivariate<NextNode1, NextNode2>;

    constexpr bool node1_has_0_deriv =
        !equal_or_depends_many<NextNode1, ActiveLeafsAndRootsDerivatives...>();
    constexpr bool node2_has_0_deriv =
        !equal_or_depends_many<NextNode2, ActiveLeafsAndRootsDerivatives...>();
    static_assert(!node1_has_0_deriv || !node2_has_0_deriv);

    // ideally std::is_same_v<Node1, Node2> should be differentiated as type
    if constexpr (std::is_same_v<NextNode1, NextNode2> || node1_has_0_deriv ||
                  node2_has_0_deriv) {
        double const current_node_d =
            std::is_same_v<NextNode1, NextNode2>
                ? 2 * CurrentNode::d1 /* could be d1 or d2 */ (
                          val(in, CurrentNode{}), val(in, NextNode1{}),
                          val(in, NextNode2{}))
            : node1_has_0_deriv
                ? CurrentNode::d2(val(in, CurrentNode{}), val(in, NextNode1{}),
                                  val(in, NextNode2{}))
                : CurrentNode::d1(val(in, CurrentNode{}), val(in, NextNode1{}),
                                  val(in, NextNode2{}));

        using NextNodeAlive =
            std::conditional_t<node1_has_0_deriv, NextNode2, NextNode1>;

        univariate_bwd<CurrentNode, NextNodeAlive>(
            in, current_node_d, std::tuple<NodesAlive...>{},
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

        using NodesAlive2 = decltype(replace_first2<position_root, available_t>(
            std::tuple<NodesAlive...>{}));

        constexpr bool next_node1_needs_new_storing =
            !is_next_node1_leaf && !has_type<NextNode1, NodesAlive...>();

        constexpr auto position_next_node1_available =
            get_idx_first2<available_t>(NodesAlive2{});

        using NodesAlive3 = std::conditional_t<
            next_node1_needs_new_storing,
            decltype(replace_first2<position_next_node1_available, NextNode1>(
                NodesAlive2{})),
            NodesAlive2>;

        constexpr bool next_node2_needs_new_storing =
            !is_next_node2_leaf && !has_type<NextNode2, NodesAlive...>();

        constexpr auto position_next_node2_available =
            get_idx_first2<available_t>(NodesAlive3{});

        using NodesAlive4 = std::conditional_t<
            next_node2_needs_new_storing,
            decltype(replace_first2<position_next_node2_available, NextNode2>(
                NodesAlive3{})),
            NodesAlive3>;

        constexpr auto position_next1 =
            get_idx_first2<NextNode1>(NodesAlive4{});
        constexpr auto position_next2 =
            get_idx_first2<NextNode2>(NodesAlive4{});

        double &next_node1_val = is_next_node1_leaf
                                     ? in.der(NextNode1{})
                                     : in.m_buffer[position_next1];

        double &next_node2_val = is_next_node2_leaf
                                     ? in.der(NextNode2{})
                                     : in.m_buffer[position_next2];

        double const current_node_d = is_current_node_root
                                          ? in.der(CurrentNode{})
                                          : in.m_buffer[position_root];

        auto d = CurrentNode::d(val(in, CurrentNode{}), val(in, NextNode1{}),
                                val(in, NextNode2{}));
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

        evaluate_bwd(in, NodesAlive4{}, std::tuple<NodesToCalc...>{});
    }
}

template <class... Values, class... ActiveLeafsAndRootsDerivatives,
          class... NodesAlive, template <class> class Univariate, class Node,
          class... NodesToCalc>
inline void
evaluate_bwd(detail::Tape2<std::tuple<Values...>,
                           std::tuple<ActiveLeafsAndRootsDerivatives...>> &in,
             std::tuple<NodesAlive...> /* deriv_ids */,
             std::tuple<Univariate<Node>, NodesToCalc...> /* nodes */) {
    using CurrentNode = Univariate<Node>;

    double const d = CurrentNode::d(val(in, CurrentNode{}), val(in, Node{}));

    univariate_bwd<CurrentNode, Node>(in, d, std::tuple<NodesAlive...>{},
                                      std::tuple<NodesToCalc...>{});
}

} // namespace detail

template <class... Values, class... ActiveLeafsAndRootsDerivatives>
inline void
evaluate_bwd(detail::Tape2<std::tuple<Values...>,
                           std::tuple<ActiveLeafsAndRootsDerivatives...>> &in) {
    constexpr auto nodes_ordered =
        detail::calc_order_t<false, ActiveLeafsAndRootsDerivatives...>();
    detail::evaluate_bwd(
        in,
        typename generate_tuple_type<
            detail::available_t,
            tape_size_bwd<ActiveLeafsAndRootsDerivatives...>()>::type{},
        nodes_ordered);
}

} // namespace adhoc

#endif // ADHOC_EVALUATE_BWD_HPP