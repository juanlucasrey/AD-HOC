#ifndef ADHOC_EVALUATE_HPP
#define ADHOC_EVALUATE_HPP

#include "adhoc.hpp"
#include "calc_order.hpp"
#include "constants_type.hpp"
#include "tape_nodes.hpp"
#include "tape_size_fwd.hpp"

namespace adhoc2 {

namespace detail {

template <constants::ArgType D, class... RootsAndLeafs,
          class... IntermediateNodes>
inline auto get(Tape<RootsAndLeafs...> const & /* in */,
                Tape<IntermediateNodes...> const & /* intermediate */,
                constants::CD<D> /* node */) -> double {
    return constants::CD<D>::v();
}

template <class Node, class... RootsAndLeafs, class... IntermediateNodes>
inline auto get(Tape<RootsAndLeafs...> const &in,
                Tape<IntermediateNodes...> const &intermediate, Node /* node */)
    -> double {
    constexpr bool is_root_or_leaf = has_type2<Node, RootsAndLeafs...>();
    if constexpr (is_root_or_leaf) {
        return in.get(Node{});
    } else {
        return intermediate.get(Node{});
    }
}

template <class... RootsAndLeafs, class... IntermediateNodes>
inline void evaluate_fwd(Tape<RootsAndLeafs...> & /* in */,
                         Tape<IntermediateNodes...> & /* intermediate */,
                         std::tuple<> /* nodes */) {}

template <class... RootsAndLeafs, class... IntermediateNodes,
          template <class, class> class Bivariate, class Node1, class Node2,
          class... NodesToCalc>
inline void
evaluate_fwd(Tape<RootsAndLeafs...> &in,
             Tape<IntermediateNodes...> &intermediate,
             std::tuple<Bivariate<Node1, Node2>, NodesToCalc...> nodes);

template <class... RootsAndLeafs, class... IntermediateNodes,
          template <class> class Univariate, class Node,
          class... IntermediateTypesToCalc>
inline void evaluate_fwd(
    Tape<RootsAndLeafs...> &in, Tape<IntermediateNodes...> &intermediate,
    std::tuple<Univariate<Node>, IntermediateTypesToCalc...> /* nodes */) {
    using this_type = Univariate<Node>;
    constexpr bool is_root_or_leaf = has_type2<this_type, RootsAndLeafs...>();
    if constexpr (is_root_or_leaf) {
        in.set(Univariate<Node>{},
               Univariate<Node>::v(get(in, intermediate, Node{})));
    } else {
        intermediate.set(Univariate<Node>{},
                         Univariate<Node>::v(get(in, intermediate, Node{})));
    }

    evaluate_fwd(in, intermediate, std::tuple<IntermediateTypesToCalc...>{});
}

template <class... RootsAndLeafs, class... IntermediateNodes,
          template <class, class> class Bivariate, class Node1, class Node2,
          class... NodesToCalc>
inline void
evaluate_fwd(Tape<RootsAndLeafs...> &in,
             Tape<IntermediateNodes...> &intermediate,
             std::tuple<Bivariate<Node1, Node2>, NodesToCalc...> /* nodes */) {
    using this_type = Bivariate<Node1, Node2>;
    constexpr bool is_root_or_leaf = has_type2<this_type, RootsAndLeafs...>();
    if constexpr (is_root_or_leaf) {
        in.set(Bivariate<Node1, Node2>{},
               Bivariate<Node1, Node2>::v(get(in, intermediate, Node1{}),
                                          get(in, intermediate, Node2{})));
    } else {
        intermediate.set(
            Bivariate<Node1, Node2>{},
            Bivariate<Node1, Node2>::v(get(in, intermediate, Node1{}),
                                       get(in, intermediate, Node2{})));
    }

    evaluate_fwd(in, intermediate, std::tuple<NodesToCalc...>{});
}

template <class... RootsAndLeafs, class... IntermediateNodes>
inline void evaluate_aux(Tape<RootsAndLeafs...> &in,
                         std::tuple<IntermediateNodes...> /* types */) {
    auto constexpr calcs =
        detail::calc_order_aux_t<true, RootsAndLeafs...>::template call();
    auto intermediate = Tape(IntermediateNodes{}...);
    evaluate_fwd(in, intermediate, calcs);
}

template <class... RootsAndLeafs, class... IntermediateNodes>
inline auto
evaluate_fwd_return_vals_aux(Tape<RootsAndLeafs...> &in,
                             std::tuple<IntermediateNodes...> /* types */) {
    auto constexpr calcs =
        detail::calc_order_aux_t<true, RootsAndLeafs...>::template call();
    auto intermediate = Tape(IntermediateNodes{}...);
    evaluate_fwd(in, intermediate, calcs);
    return intermediate;
}

} // namespace detail

template <class... RootsAndLeafs>
inline void evaluate(Tape<RootsAndLeafs...> &in) {
    auto constexpr nodes_inside = tape_size_fwd_t(RootsAndLeafs{}...);
    detail::evaluate_aux(in, nodes_inside);
}

template <class... RootsAndLeafs>
inline auto evaluate_fwd_return_vals(Tape<RootsAndLeafs...> &in) {
    auto constexpr nodes_inside = tape_size_fwd_t(RootsAndLeafs{}...);
    return detail::evaluate_fwd_return_vals_aux(in, nodes_inside);
}

} // namespace adhoc2

#endif // ADHOC_EVALUATE_HPP