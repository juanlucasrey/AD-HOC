#ifndef ADHOC_EVALUATE_FWD_HPP
#define ADHOC_EVALUATE_FWD_HPP

#include "adhoc.hpp"
#include "calc_order.hpp"
#include "constants_type.hpp"
#include "tape.hpp"
#include "tape_size_fwd.hpp"

namespace adhoc {

namespace detail {

template <class Tape>
inline void evaluate_fwd2(Tape & /* in */, std::tuple<> /* nodes */) {}

template <class Tape, template <class, class> class Bivariate, class Node1,
          class Node2, class... NodesToCalc>
inline void
evaluate_fwd2(Tape &in,
              std::tuple<Bivariate<Node1, Node2>, NodesToCalc...> nodes);

template <class Tape, template <class> class Univariate, class Node,
          class... NodesToCalc>
inline void
evaluate_fwd2(Tape &in,
              std::tuple<Univariate<Node>, NodesToCalc...> /* nodes */) {
    using CurrentNode = Univariate<Node>;
    in.val(CurrentNode{}) = CurrentNode::v(val(in, Node{}));
    evaluate_fwd2(in, std::tuple<NodesToCalc...>{});
}

template <class Tape, template <class, class> class Bivariate, class Node1,
          class Node2, class... NodesToCalc>
inline void
evaluate_fwd2(Tape &in,
              std::tuple<Bivariate<Node1, Node2>, NodesToCalc...> /* nodes */) {
    using CurrentNode = Bivariate<Node1, Node2>;
    in.val(CurrentNode{}) = CurrentNode::v(val(in, Node1{}), val(in, Node2{}));
    evaluate_fwd2(in, std::tuple<NodesToCalc...>{});
}

} // namespace detail

template <class... Values, class... ActiveLeafsAndRootsDerivatives,
          std::size_t MaxWidth>
inline auto
evaluate_fwd(detail::Tape2<std::tuple<Values...>,
                           std::tuple<ActiveLeafsAndRootsDerivatives...>,
                           MaxWidth> &in) {
    auto constexpr calcs = detail::calc_order_aux_t<
        true, ActiveLeafsAndRootsDerivatives...>::template call();
    evaluate_fwd2(in, calcs);
}

} // namespace adhoc

#endif // ADHOC_EVALUATE_FWD_HPP