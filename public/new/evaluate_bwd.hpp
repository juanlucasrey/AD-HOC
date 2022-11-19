#ifndef ADHOC_EVALUATE_BWD_HPP
#define ADHOC_EVALUATE_BWD_HPP

// #include "adhoc.hpp"
// #include "calc_order.hpp"
// #include "constants_type.hpp"
#include "tape_nodes.hpp"
#include "tape_size_bwd.hpp"

#include <array>

namespace adhoc2 {

namespace detail {} // namespace detail

template <class... RootsAndLeafs, class... IntermediateNodes,
          class... ActiveRootsAndLeafs>
inline void evaluate_bwd(Tape<RootsAndLeafs...> const & /* in */,
                         Tape<IntermediateNodes...> const & /* intermediate */,
                         Tape<ActiveRootsAndLeafs...> & /* derivs */) {
    constexpr auto intermediate_deriv_tape_size =
        tape_size_bwd<ActiveRootsAndLeafs...>();
    std::array<double, intermediate_deriv_tape_size> intermediateDerivs{};
}

} // namespace adhoc2

#endif // ADHOC_EVALUATE_BWD_HPP