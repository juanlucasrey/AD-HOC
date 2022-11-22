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
    constexpr bool is_input = has_type2<Input, InputTypes...>();
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

template <class this_type, class Node, class... RootsAndLeafs,
          class... IntermediateNodes, class... ActiveRootsAndLeafs,
          std::size_t N, class... NodesAlive, class... NodesToCalc>
inline void univariate_bwd(Tape<RootsAndLeafs...> const &in,
                           Tape<IntermediateNodes...> const &intermediate,
                           Tape<ActiveRootsAndLeafs...> &derivs,
                           std::array<double, N> &deriv_vals, const double &d,
                           std::tuple<NodesAlive...> /* deriv_ids */,
                           std::tuple<NodesToCalc...> /* nodes */) {
    constexpr bool is_next_node_leaf =
        has_type2<Node, ActiveRootsAndLeafs...>();
    constexpr bool is_current_node_root =
        has_type2<this_type, ActiveRootsAndLeafs...>();

    if constexpr (is_next_node_leaf) {
        if constexpr (is_current_node_root) {
            derivs.get(Node{}) += derivs.get(this_type{}) * d;
            evaluate_bwd(in, intermediate, derivs, deriv_vals,
                         std::tuple<NodesAlive...>{},
                         std::tuple<NodesToCalc...>{});
        } else {
            auto constexpr pos =
                get_idx_first<this_type>(std::tuple<NodesAlive...>{});
            auto constexpr deriv_ids_new =
                replace_first<this_type, available_t>(
                    std::tuple<NodesAlive...>{});
            derivs.get(Node{}) += deriv_vals[pos] * d;
            evaluate_bwd(in, intermediate, derivs, deriv_vals, deriv_ids_new,
                         std::tuple<NodesToCalc...>{});
        }
    } else {
        if constexpr (is_current_node_root) {
            auto constexpr pos =
                get_idx_first<available_t>(std::tuple<NodesAlive...>{});
            auto constexpr deriv_ids_new =
                replace_first<available_t, Node>(std::tuple<NodesAlive...>{});
            deriv_vals[pos] = derivs.get(this_type{}) * d;
            evaluate_bwd(in, intermediate, derivs, deriv_vals, deriv_ids_new,
                         std::tuple<NodesToCalc...>{});
        } else {
            auto constexpr pos =
                get_idx_first<this_type>(std::tuple<NodesAlive...>{});
            auto constexpr deriv_ids_new =
                replace_first<this_type, Node>(std::tuple<NodesAlive...>{});
            deriv_vals[pos] *= d;
            evaluate_bwd(in, intermediate, derivs, deriv_vals, deriv_ids_new,
                         std::tuple<NodesToCalc...>{});
        }
    }
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

    constexpr bool is_current_node_root =
        has_type2<this_type, ActiveRootsAndLeafs...>();

    if constexpr (!node1_has_0_deriv && !node2_has_0_deriv) {
        double const d1 = this_type::d1(get2(in, intermediate, this_type{}),
                                        get2(in, intermediate, Node1{}),
                                        get2(in, intermediate, Node2{}));

        double const d2 = this_type::d2(get2(in, intermediate, this_type{}),
                                        get2(in, intermediate, Node1{}),
                                        get2(in, intermediate, Node2{}));
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
    // constexpr bool is_next_node_leaf =
    //     has_type2<Node, ActiveRootsAndLeafs...>();
    // constexpr bool is_current_node_root =
    //     has_type2<this_type, ActiveRootsAndLeafs...>();

    // if constexpr (is_next_node_leaf) {
    //     if constexpr (is_current_node_root) {
    //         derivs.get(Node{}) += derivs.get(this_type{}) * d;
    //         evaluate_bwd(in, intermediate, derivs, deriv_vals,
    //                      std::tuple<NodesAlive...>{},
    //                      std::tuple<NodesToCalc...>{});
    //     } else {
    //         auto constexpr pos =
    //             get_idx_first<this_type>(std::tuple<NodesAlive...>{});
    //         auto constexpr deriv_ids_new =
    //             replace_first<this_type, available_t>(
    //                 std::tuple<NodesAlive...>{});
    //         derivs.get(Node{}) += deriv_vals[pos] * d;
    //         evaluate_bwd(in, intermediate, derivs, deriv_vals, deriv_ids_new,
    //                      std::tuple<NodesToCalc...>{});
    //     }
    // } else {
    //     if constexpr (is_current_node_root) {
    //         auto constexpr pos =
    //             get_idx_first<available_t>(std::tuple<NodesAlive...>{});
    //         auto constexpr deriv_ids_new =
    //             replace_first<available_t,
    //             Node>(std::tuple<NodesAlive...>{});
    //         deriv_vals[pos] = derivs.get(this_type{}) * d;
    //         evaluate_bwd(in, intermediate, derivs, deriv_vals, deriv_ids_new,
    //                      std::tuple<NodesToCalc...>{});
    //     } else {
    //         auto constexpr pos =
    //             get_idx_first<this_type>(std::tuple<NodesAlive...>{});
    //         auto constexpr deriv_ids_new =
    //             replace_first<this_type, Node>(std::tuple<NodesAlive...>{});
    //         deriv_vals[pos] *= d;
    //         evaluate_bwd(in, intermediate, derivs, deriv_vals, deriv_ids_new,
    //                      std::tuple<NodesToCalc...>{});
    //     }
    // }
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