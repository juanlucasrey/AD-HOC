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
          template <class> class Univariate, class Node, class... NodesToCalc>
inline void
evaluate_bwd(Tape<RootsAndLeafs...> const &in,
             Tape<IntermediateNodes...> const &intermediate,
             Tape<ActiveRootsAndLeafs...> &derivs,
             std::array<double, N> &deriv_vals,
             std::tuple<NodesAlive...> /* deriv_ids */,
             std::tuple<Univariate<Node>, NodesToCalc...> /* nodes */);

template <class... RootsAndLeafs, class... IntermediateNodes,
          class... ActiveRootsAndLeafs, std::size_t N, class... NodesAlive,
          template <class, class> class Bivariate, class Node1, class Node2,
          class... NodesToCalc>
inline void
evaluate_bwd(Tape<RootsAndLeafs...> const &in,
             Tape<IntermediateNodes...> const &intermediate,
             Tape<ActiveRootsAndLeafs...> &derivs,
             std::array<double, N> &deriv_vals,
             std::tuple<NodesAlive...> /* deriv_ids */,
             std::tuple<Bivariate<Node1, Node2>, NodesToCalc...> /* nodes */);

template <class CurrentNode, class NextNode, class... RootsAndLeafs,
          class... IntermediateNodes, class... ActiveRootsAndLeafs,
          std::size_t N, class... NodesAlive, class... NodesToCalc>
inline void univariate_bwd(Tape<RootsAndLeafs...> const &in,
                           Tape<IntermediateNodes...> const &intermediate,
                           Tape<ActiveRootsAndLeafs...> &derivs,
                           std::array<double, N> &deriv_vals, const double &d,
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

    double &next_node_val =
        is_next_node_leaf ? derivs.get2(NextNode{}) : deriv_vals[position_next];

    double const current_node_d = is_current_node_root
                                      ? derivs.get2(CurrentNode{})
                                      : deriv_vals[position_root];
    if constexpr (next_node_needs_new_storing) {
        next_node_val = current_node_d * d;
    } else {
        next_node_val += current_node_d * d;
    }

    evaluate_bwd(in, intermediate, derivs, deriv_vals, NodesAlive4{},
                 std::tuple<NodesToCalc...>{});
}

template <class... RootsAndLeafs, class... IntermediateNodes,
          class... ActiveRootsAndLeafs, std::size_t N, class... NodesAlive,
          template <class, class> class Bivariate, class Node1, class Node2,
          class... NodesToCalc>
inline void
evaluate_bwd(Tape<RootsAndLeafs...> const &in,
             Tape<IntermediateNodes...> const &intermediate,
             Tape<ActiveRootsAndLeafs...> &derivs,
             std::array<double, N> &deriv_vals,
             std::tuple<NodesAlive...> /* deriv_ids */,
             std::tuple<Bivariate<Node1, Node2>, NodesToCalc...> /* nodes */) {
    using this_type = Bivariate<Node1, Node2>;

    constexpr bool node1_has_0_deriv =
        !equal_or_depends_many<Node1, ActiveRootsAndLeafs...>();
    constexpr bool node2_has_0_deriv =
        !equal_or_depends_many<Node2, ActiveRootsAndLeafs...>();
    static_assert(!node1_has_0_deriv || !node2_has_0_deriv);

    // ideally this should be differentiated at the type level
    if constexpr (std::is_same_v<Node1, Node2>) {
        static_assert(!node1_has_0_deriv && !node2_has_0_deriv);
        // could be d1 or d2
        double const d = 2 * this_type::d1(get2(in, intermediate, this_type{}),
                                           get2(in, intermediate, Node1{}),
                                           get2(in, intermediate, Node2{}));

        // could be Node1 or Node2
        univariate_bwd<this_type, Node1>(in, intermediate, derivs, deriv_vals,
                                         d, std::tuple<NodesAlive...>{},
                                         std::tuple<NodesToCalc...>{});
    } else if constexpr (!node1_has_0_deriv && !node2_has_0_deriv) {
        constexpr bool is_next_node1_leaf =
            has_type<Node1, ActiveRootsAndLeafs...>();
        constexpr bool is_next_node2_leaf =
            has_type<Node2, ActiveRootsAndLeafs...>();
        constexpr bool is_current_node_root =
            has_type<this_type, ActiveRootsAndLeafs...>();

        double const d1 = this_type::d1(get2(in, intermediate, this_type{}),
                                        get2(in, intermediate, Node1{}),
                                        get2(in, intermediate, Node2{}));

        double const d2 = this_type::d2(get2(in, intermediate, this_type{}),
                                        get2(in, intermediate, Node1{}),
                                        get2(in, intermediate, Node2{}));

        auto constexpr position_root2 =
            get_idx_first2<this_type>(std::tuple<NodesAlive...>{});
        double const current_node_d = is_current_node_root
                                          ? derivs.get2(this_type{})
                                          : deriv_vals[position_root2];

        if constexpr (is_current_node_root) {
            if constexpr (is_next_node1_leaf) {
                derivs.get(Node1{}) += current_node_d * d1;
                if constexpr (is_next_node2_leaf) {
                    derivs.get(Node2{}) += current_node_d * d2;
                    evaluate_bwd(in, intermediate, derivs, deriv_vals,
                                 std::tuple<NodesAlive...>{},
                                 std::tuple<NodesToCalc...>{});
                } else {
                    constexpr bool is_next_node2_stored =
                        has_type<Node2, NodesAlive...>();
                    if constexpr (is_next_node2_stored) {
                        auto constexpr position_node2 =
                            get_idx_first<Node2>(std::tuple<NodesAlive...>{});
                        deriv_vals[position_node2] +=
                            derivs.get(this_type{}) * d2;
                        evaluate_bwd(in, intermediate, derivs, deriv_vals,
                                     std::tuple<NodesAlive...>{},
                                     std::tuple<NodesToCalc...>{});
                    } else {
                        auto constexpr position_node2 =
                            get_idx_first<available_t>(
                                std::tuple<NodesAlive...>{});
                        auto constexpr deriv_ids_new =
                            replace_first<available_t, Node2>(
                                std::tuple<NodesAlive...>{});
                        deriv_vals[position_node2] =
                            derivs.get(this_type{}) * d2;
                        evaluate_bwd(in, intermediate, derivs, deriv_vals,
                                     deriv_ids_new,
                                     std::tuple<NodesToCalc...>{});
                    }
                }
            } else {
                if constexpr (is_next_node2_leaf) {
                    derivs.get(Node2{}) += derivs.get(this_type{}) * d2;

                    constexpr bool is_next_node1_stored =
                        has_type<Node1, NodesAlive...>();
                    if constexpr (is_next_node1_stored) {
                        auto constexpr position_node1 =
                            get_idx_first<Node1>(std::tuple<NodesAlive...>{});
                        deriv_vals[position_node1] +=
                            derivs.get(this_type{}) * d1;
                        evaluate_bwd(in, intermediate, derivs, deriv_vals,
                                     std::tuple<NodesAlive...>{},
                                     std::tuple<NodesToCalc...>{});
                    } else {
                        auto constexpr position_node1 =
                            get_idx_first<available_t>(
                                std::tuple<NodesAlive...>{});
                        auto constexpr deriv_ids_new =
                            replace_first<available_t, Node1>(
                                std::tuple<NodesAlive...>{});
                        deriv_vals[position_node1] =
                            derivs.get(this_type{}) * d1;
                        evaluate_bwd(in, intermediate, derivs, deriv_vals,
                                     deriv_ids_new,
                                     std::tuple<NodesToCalc...>{});
                    }
                } else {
                    constexpr bool is_next_node1_stored =
                        has_type<Node1, NodesAlive...>();
                    constexpr bool is_next_node2_stored =
                        has_type<Node2, NodesAlive...>();
                    if constexpr (is_next_node1_stored) {
                        auto constexpr position_node1 =
                            get_idx_first<Node1>(std::tuple<NodesAlive...>{});
                        deriv_vals[position_node1] +=
                            derivs.get(this_type{}) * d1;
                        if constexpr (is_next_node2_stored) {
                            auto constexpr position_node2 =
                                get_idx_first<Node2>(
                                    std::tuple<NodesAlive...>{});
                            deriv_vals[position_node2] +=
                                derivs.get(this_type{}) * d2;
                            evaluate_bwd(in, intermediate, derivs, deriv_vals,
                                         std::tuple<NodesAlive...>{},
                                         std::tuple<NodesToCalc...>{});
                        } else {
                            auto constexpr position_node2 =
                                get_idx_first<available_t>(
                                    std::tuple<NodesAlive...>{});
                            auto constexpr deriv_ids_new =
                                replace_first<available_t, Node2>(
                                    std::tuple<NodesAlive...>{});
                            deriv_vals[position_node2] =
                                derivs.get(this_type{}) * d2;
                            evaluate_bwd(in, intermediate, derivs, deriv_vals,
                                         deriv_ids_new,
                                         std::tuple<NodesToCalc...>{});
                        }
                    } else {
                        auto constexpr position_node1 =
                            get_idx_first<available_t>(
                                std::tuple<NodesAlive...>{});
                        auto constexpr deriv_ids_new =
                            replace_first<available_t, Node1>(
                                std::tuple<NodesAlive...>{});
                        deriv_vals[position_node1] =
                            derivs.get(this_type{}) * d1;
                        if constexpr (is_next_node2_stored) {
                            auto constexpr position_node2 =
                                get_idx_first<Node2>(
                                    std::tuple<NodesAlive...>{});
                            deriv_vals[position_node2] +=
                                derivs.get(this_type{}) * d2;
                            evaluate_bwd(in, intermediate, derivs, deriv_vals,
                                         deriv_ids_new,
                                         std::tuple<NodesToCalc...>{});
                        } else {
                            auto constexpr position_node2 =
                                get_idx_first<available_t>(deriv_ids_new);
                            auto constexpr deriv_ids_new2 =
                                replace_first<available_t, Node2>(
                                    deriv_ids_new);
                            deriv_vals[position_node2] =
                                derivs.get(this_type{}) * d2;
                            evaluate_bwd(in, intermediate, derivs, deriv_vals,
                                         deriv_ids_new2,
                                         std::tuple<NodesToCalc...>{});
                        }
                    }
                }
            }
        } else {
            auto constexpr position_root =
                get_idx_first<this_type>(std::tuple<NodesAlive...>{});
            auto constexpr deriv_ids_new =
                replace_first<this_type, available_t>(
                    std::tuple<NodesAlive...>{});
            if constexpr (is_next_node1_leaf) {
                derivs.get(Node1{}) += deriv_vals[position_root] * d1;
                if constexpr (is_next_node2_leaf) {
                    derivs.get(Node2{}) += deriv_vals[position_root] * d2;
                    evaluate_bwd(in, intermediate, derivs, deriv_vals,
                                 deriv_ids_new, std::tuple<NodesToCalc...>{});
                } else {
                    constexpr bool is_next_node2_stored =
                        has_type<Node2, NodesAlive...>();
                    if constexpr (is_next_node2_stored) {
                        auto constexpr position_node2 =
                            get_idx_first<Node2>(deriv_ids_new);
                        deriv_vals[position_node2] +=
                            deriv_vals[position_root] * d2;
                        evaluate_bwd(in, intermediate, derivs, deriv_vals,
                                     deriv_ids_new,
                                     std::tuple<NodesToCalc...>{});
                    } else {
                        auto constexpr position_node2 =
                            get_idx_first<available_t>(deriv_ids_new);
                        auto constexpr deriv_ids_new2 =
                            replace_first<available_t, Node2>(deriv_ids_new);
                        deriv_vals[position_node2] =
                            deriv_vals[position_root] * d2;
                        evaluate_bwd(in, intermediate, derivs, deriv_vals,
                                     deriv_ids_new2,
                                     std::tuple<NodesToCalc...>{});
                    }
                }
            } else {
                if constexpr (is_next_node2_leaf) {
                    derivs.get(Node2{}) += deriv_vals[position_root] * d2;
                    constexpr bool is_next_node1_stored =
                        has_type<Node1, NodesAlive...>();
                    if constexpr (is_next_node1_stored) {
                        auto constexpr position_node1 =
                            get_idx_first<Node1>(deriv_ids_new);
                        deriv_vals[position_node1] +=
                            deriv_vals[position_root] * d1;
                        evaluate_bwd(in, intermediate, derivs, deriv_vals,
                                     deriv_ids_new,
                                     std::tuple<NodesToCalc...>{});
                    } else {
                        auto constexpr position_node1 =
                            get_idx_first<available_t>(deriv_ids_new);
                        auto constexpr deriv_ids_new2 =
                            replace_first<available_t, Node1>(deriv_ids_new);
                        deriv_vals[position_node1] =
                            deriv_vals[position_root] * d1;
                        evaluate_bwd(in, intermediate, derivs, deriv_vals,
                                     deriv_ids_new2,
                                     std::tuple<NodesToCalc...>{});
                    }
                } else {
                    constexpr bool is_next_node1_stored =
                        has_type<Node1, NodesAlive...>();
                    constexpr bool is_next_node2_stored =
                        has_type<Node2, NodesAlive...>();
                    if constexpr (is_next_node1_stored) {
                        auto constexpr position_node1 =
                            get_idx_first<Node1>(deriv_ids_new);
                        deriv_vals[position_node1] +=
                            deriv_vals[position_root] * d1;
                        if constexpr (is_next_node2_stored) {
                            auto constexpr position_node2 =
                                get_idx_first<Node2>(deriv_ids_new);
                            deriv_vals[position_node2] +=
                                deriv_vals[position_root] * d2;
                            evaluate_bwd(in, intermediate, derivs, deriv_vals,
                                         deriv_ids_new,
                                         std::tuple<NodesToCalc...>{});
                        } else {
                            auto constexpr position_node2 =
                                get_idx_first<available_t>(deriv_ids_new);
                            auto constexpr deriv_ids_new2 =
                                replace_first<available_t, Node2>(
                                    deriv_ids_new);
                            deriv_vals[position_node2] =
                                deriv_vals[position_root] * d2;
                            evaluate_bwd(in, intermediate, derivs, deriv_vals,
                                         deriv_ids_new2,
                                         std::tuple<NodesToCalc...>{});
                        }
                    } else {
                        auto constexpr position_node1 =
                            get_idx_first<available_t>(deriv_ids_new);
                        auto constexpr deriv_ids_new2 =
                            replace_first<available_t, Node1>(deriv_ids_new);
                        if constexpr (is_next_node2_stored) {
                            auto constexpr position_node2 =
                                get_idx_first<Node2>(deriv_ids_new2);
                            deriv_vals[position_node2] +=
                                deriv_vals[position_root] * d2;

                            if constexpr (position_node1 == position_root) {
                                deriv_vals[position_node1] *= d1;
                            } else {
                                deriv_vals[position_node1] =
                                    deriv_vals[position_root] * d1;
                            }

                            evaluate_bwd(in, intermediate, derivs, deriv_vals,
                                         deriv_ids_new2,
                                         std::tuple<NodesToCalc...>{});
                        } else {
                            auto constexpr position_node2 =
                                get_idx_first<available_t>(deriv_ids_new2);
                            auto constexpr deriv_ids_new3 =
                                replace_first<available_t, Node2>(
                                    deriv_ids_new2);
                            deriv_vals[position_node2] =
                                deriv_vals[position_root] * d2;

                            if constexpr (position_node1 == position_root) {
                                deriv_vals[position_node1] *= d1;
                            } else {
                                deriv_vals[position_node1] =
                                    deriv_vals[position_root] * d1;
                            }
                            evaluate_bwd(in, intermediate, derivs, deriv_vals,
                                         deriv_ids_new3,
                                         std::tuple<NodesToCalc...>{});
                        }
                    }
                }
            }
        }
    } else if constexpr (node1_has_0_deriv) {
        double const d = this_type::d2(get2(in, intermediate, this_type{}),
                                       get2(in, intermediate, Node1{}),
                                       get2(in, intermediate, Node2{}));
        univariate_bwd<this_type, Node2>(in, intermediate, derivs, deriv_vals,
                                         d, std::tuple<NodesAlive...>{},
                                         std::tuple<NodesToCalc...>{});
    } else if constexpr (node2_has_0_deriv) {
        double const d = this_type::d1(get2(in, intermediate, this_type{}),
                                       get2(in, intermediate, Node1{}),
                                       get2(in, intermediate, Node2{}));
        univariate_bwd<this_type, Node1>(in, intermediate, derivs, deriv_vals,
                                         d, std::tuple<NodesAlive...>{},
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
    using this_type = Univariate<Node>;

    double const d = this_type::d(get2(in, intermediate, this_type{}),
                                  get2(in, intermediate, Node{}));

    univariate_bwd<this_type, Node>(in, intermediate, derivs, deriv_vals, d,
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