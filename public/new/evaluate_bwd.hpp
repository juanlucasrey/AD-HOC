#ifndef ADHOC_EVALUATE_BWD_HPP
#define ADHOC_EVALUATE_BWD_HPP

// #include "adhoc.hpp"
#include "calc_order.hpp"
#include "constants_type.hpp"
#include "tape_nodes.hpp"
#include "tape_size_bwd.hpp"
#include "tuple_utils.hpp"

#include <array>

namespace adhoc2 {

namespace detail {

// consolidate later on separate file
template <constants::ArgType D, class... InputTypes, class... IntermediateTypes>
inline auto get2(Tape<InputTypes...> const &,
                 Tape<IntermediateTypes...> const &, constants::CD<D> /* in */)
    -> double {
    return constants::CD<D>::v();
}

template <class Input, class... InputTypes, class... IntermediateTypes>
inline auto get2(Tape<InputTypes...> const &in,
                Tape<IntermediateTypes...> const &intermediate, Input /* in
                */)
    -> double {
    constexpr bool is_input = has_type<Input, InputTypes...>();
    if constexpr (is_input) {
        return in.get(Input{});
    } else {
        return intermediate.get(Input{});
    }
}

class available_t {};

template <class... RootsAndLeafs, class... IntermediateNodes,
          class... ActiveRootsAndLeafs, std::size_t N, class... NodesAlive>
inline void evaluate_bwd(Tape<RootsAndLeafs...> const & /* in */,
                         Tape<IntermediateNodes...> const & /* intermediate */,
                         Tape<ActiveRootsAndLeafs...> & /* derivs */,
                         std::array<double, N> & /* derivs_intermediate */,
                         std::tuple<NodesAlive...> /* nodes_intermediate */,
                         std::tuple<> /* nodes */) {}

template <class... RootsAndLeafs, class... IntermediateNodes,
          class... ActiveRootsAndLeafs, std::size_t N, class... NodesAlive,
          template <class> class Univariate, class NextNode,
          class... NodesToCalc>
inline void
evaluate_bwd(Tape<RootsAndLeafs...> const &in,
             Tape<IntermediateNodes...> const &intermediate,
             Tape<ActiveRootsAndLeafs...> &derivs,
             std::array<double, N> &deriv_vals,
             std::tuple<NodesAlive...> /* deriv_ids */,
             std::tuple<Univariate<NextNode>, NodesToCalc...> /* nodes */);

template <class... RootsAndLeafs, class... IntermediateNodes,
          class... ActiveRootsAndLeafs, std::size_t N, class... NodesAlive,
          template <class, class> class Bivariate, class NextNode1,
          class NextNode2, class... NodesToCalc>
inline void evaluate_bwd(
    Tape<RootsAndLeafs...> const &in,
    Tape<IntermediateNodes...> const &intermediate,
    Tape<ActiveRootsAndLeafs...> &derivs, std::array<double, N> &deriv_vals,
    std::tuple<NodesAlive...> /* deriv_ids */,
    std::tuple<Bivariate<NextNode1, NextNode2>, NodesToCalc...> /* nodes */);

template <class CurrentNode, class NextNode, class... RootsAndLeafs,
          class... IntermediateNodes, class... ActiveRootsAndLeafs,
          std::size_t N, class... NodesAlive, class... NodesToCalc>
inline void univariate_bwd(Tape<RootsAndLeafs...> const &in,
                           Tape<IntermediateNodes...> const &intermediate,
                           Tape<ActiveRootsAndLeafs...> &derivs,
                           std::array<double, N> &deriv_vals,
                           const double &current_node_d,
                           std::tuple<NodesAlive...> /* deriv_ids */,
                           std::tuple<NodesToCalc...> /* nodes */) {
    constexpr bool is_next_node_leaf =
        has_type<NextNode, ActiveRootsAndLeafs...>();
    constexpr bool is_current_node_root =
        has_type<CurrentNode, ActiveRootsAndLeafs...>();

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
        is_next_node_leaf ? derivs.get2(NextNode{}) : deriv_vals[position_next];

    double const next_node_d =
        current_node_d * (is_current_node_root ? derivs.get2(CurrentNode{})
                                               : deriv_vals[position_root]);
    if constexpr (next_node_needs_new_storing) {
        next_node_total_d = next_node_d;
    } else {
        next_node_total_d += next_node_d;
    }

    evaluate_bwd(in, intermediate, derivs, deriv_vals, NodesAlive4{},
                 std::tuple<NodesToCalc...>{});
}

template <class... RootsAndLeafs, class... IntermediateNodes,
          class... ActiveRootsAndLeafs, std::size_t N, class... NodesAlive,
          template <class, class> class Bivariate, class NextNode1,
          class NextNode2, class... NodesToCalc>
inline void evaluate_bwd(
    Tape<RootsAndLeafs...> const &in,
    Tape<IntermediateNodes...> const &intermediate,
    Tape<ActiveRootsAndLeafs...> &derivs, std::array<double, N> &deriv_vals,
    std::tuple<NodesAlive...> /* deriv_ids */,
    std::tuple<Bivariate<NextNode1, NextNode2>, NodesToCalc...> /* nodes */) {
    using CurrentNode = Bivariate<NextNode1, NextNode2>;

    constexpr bool node1_has_0_deriv =
        !equal_or_depends_many<NextNode1, ActiveRootsAndLeafs...>();
    constexpr bool node2_has_0_deriv =
        !equal_or_depends_many<NextNode2, ActiveRootsAndLeafs...>();
    static_assert(!node1_has_0_deriv || !node2_has_0_deriv);

    // ideally std::is_same_v<Node1, Node2> should be differentiated as type
    if constexpr (std::is_same_v<NextNode1, NextNode2> || node1_has_0_deriv ||
                  node2_has_0_deriv) {
        double const current_node_d =
            std::is_same_v<NextNode1, NextNode2>
                ? 2 * CurrentNode::d1 /* could be d1 or d2 */ (
                          get2(in, intermediate, CurrentNode{}),
                          get2(in, intermediate, NextNode1{}),
                          get2(in, intermediate, NextNode2{}))
            : node1_has_0_deriv
                ? CurrentNode::d2(get2(in, intermediate, CurrentNode{}),
                                  get2(in, intermediate, NextNode1{}),
                                  get2(in, intermediate, NextNode2{}))
                : CurrentNode::d1(get2(in, intermediate, CurrentNode{}),
                                  get2(in, intermediate, NextNode1{}),
                                  get2(in, intermediate, NextNode2{}));

        using NextNodeAlive =
            std::conditional_t<node1_has_0_deriv, NextNode2, NextNode1>;

        univariate_bwd<CurrentNode, NextNodeAlive>(
            in, intermediate, derivs, deriv_vals, current_node_d,
            std::tuple<NodesAlive...>{}, std::tuple<NodesToCalc...>{});
    } else {
        constexpr bool is_next_node1_leaf =
            has_type<NextNode1, ActiveRootsAndLeafs...>();
        constexpr bool is_next_node2_leaf =
            has_type<NextNode2, ActiveRootsAndLeafs...>();
        constexpr bool is_current_node_root =
            has_type<CurrentNode, ActiveRootsAndLeafs...>();

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
                                     ? derivs.get2(NextNode1{})
                                     : deriv_vals[position_next1];

        double &next_node2_val = is_next_node2_leaf
                                     ? derivs.get2(NextNode2{})
                                     : deriv_vals[position_next2];

        double const current_node_d = is_current_node_root
                                          ? derivs.get2(CurrentNode{})
                                          : deriv_vals[position_root];

        auto d = CurrentNode::d(get2(in, intermediate, CurrentNode{}),
                                get2(in, intermediate, NextNode1{}),
                                get2(in, intermediate, NextNode2{}));
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

        evaluate_bwd(in, intermediate, derivs, deriv_vals, NodesAlive4{},
                     std::tuple<NodesToCalc...>{});
    }
}

template <class... RootsAndLeafs, class... IntermediateNodes,
          class... ActiveRootsAndLeafs, std::size_t N, class... NodesAlive,
          template <class> class Univariate, class Node, class... NodesToCalc>
inline void
evaluate_bwd(Tape<RootsAndLeafs...> const &in,
             Tape<IntermediateNodes...> const &intermediate,
             Tape<ActiveRootsAndLeafs...> &derivs,
             std::array<double, N> &deriv_vals,
             std::tuple<NodesAlive...> /* deriv_ids */,
             std::tuple<Univariate<Node>, NodesToCalc...> /* nodes */) {
    using CurrentNode = Univariate<Node>;

    double const d = CurrentNode::d(get2(in, intermediate, CurrentNode{}),
                                    get2(in, intermediate, Node{}));

    univariate_bwd<CurrentNode, Node>(in, intermediate, derivs, deriv_vals, d,
                                      std::tuple<NodesAlive...>{},
                                      std::tuple<NodesToCalc...>{});
}

} // namespace detail

template <class... RootsAndLeafs, class... IntermediateNodes,
          class... ActiveRootsAndLeafs>
inline void evaluate_bwd(Tape<RootsAndLeafs...> const &in,
                         Tape<IntermediateNodes...> const &intermediate,
                         Tape<ActiveRootsAndLeafs...> &derivs) {
    constexpr auto intermediate_deriv_tape_size =
        tape_size_bwd<ActiveRootsAndLeafs...>();
    std::array<double, intermediate_deriv_tape_size> derivs_intermediate{};
    constexpr auto nodes_ordered =
        detail::calc_order_t<false, ActiveRootsAndLeafs...>();
    detail::evaluate_bwd(
        in, intermediate, derivs, derivs_intermediate,
        typename generate_tuple_type<detail::available_t,
                                     intermediate_deriv_tape_size>::type{},
        nodes_ordered);
}

} // namespace adhoc2

#endif // ADHOC_EVALUATE_BWD_HPP