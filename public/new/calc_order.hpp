#ifndef ADHOC_TAPE_SIZE_FWD_HPP
#define ADHOC_TAPE_SIZE_FWD_HPP

#include "adhoc.hpp"
#include "constants_type.hpp"
#include "dependency.hpp"
#include "utils.hpp"

#include <tuple>

namespace adhoc2 {

namespace detail {

template <typename... NodesAlive> struct calc_order_aux_t {
    template <typename... Operations> constexpr static auto call() noexcept;
};

template <> struct calc_order_aux_t<> {
    template <typename... Operations> constexpr static auto call() noexcept {
        return std::tuple<Operations...>{};
    }
};

#if __cplusplus >= 202002L

template <constants::detail::AsTemplateArg<double> D, typename... NodesAlive>
struct calc_order_aux_t<constants::CD<D>, NodesAlive...> {
    template <typename... Operations> constexpr static auto call() noexcept {
        // we don't add anything to operations because it's a constant
        return calc_order_aux_t<NodesAlive...>::template call<Operations...>();
    }
};

#else

template <std::uint64_t D, typename... NodesAlive>
struct calc_order_aux_t<constants::CD<D>, NodesAlive...> {
    template <typename... Operations> constexpr static auto call() noexcept {
        // we don't add anything to operations because it's a constant
        return calc_order_aux_t<NodesAlive...>::template call<Operations...>();
    }
};

#endif

template <std::size_t N, typename... NodesAlive>
struct calc_order_aux_t<double_t<N>, NodesAlive...> {
    template <typename... Operations> constexpr static auto call() noexcept {
        // we don't add anything to operations because it's an input
        return calc_order_aux_t<NodesAlive...>::template call<Operations...>();
    }
};

// template <std::size_t N, typename... NodesAlive>
// struct calc_order_aux_t<double_t<N> const, NodesAlive...> {
//     template <typename... Operations> constexpr static auto call() noexcept {
//         // we don't add anything to operations because it's an input
//         return calc_order_aux_t<NodesAlive...>::template
//         call<Operations...>();
//     }
// };

template <template <class...> class Xvariate, class... Node,
          typename... NodesAlive>
struct calc_order_aux_t<Xvariate<Node...>, NodesAlive...> {
    template <typename... Operations> constexpr static auto call() noexcept {
        using this_type = Xvariate<Node...>;
        // it is possible that some intermediate nodes are sent as input
        // static_assert(!has_type2<this_type, NodesAlive...>());

        constexpr bool other_types_depend_on_this =
            (depends<NodesAlive, this_type>::call() || ...);

        if constexpr (other_types_depend_on_this) {
            // this_type will come up again because it is included on NodesAlive
            return calc_order_aux_t<NodesAlive...>::template call<
                Operations...>();
        } else {
            return calc_order_aux_t<Node..., NodesAlive...>::template call<
                this_type, Operations...>();
        }
    }
};

// template <template <class...> class Xvariate, class... Node,
//           typename... NodesAlive>
// struct calc_order_aux_t<Xvariate<Node...> const, NodesAlive...> {
//     template <typename... Operations> constexpr static auto call() noexcept {
//         using this_type = Xvariate<Node...>;
//         static_assert(!has_type2<this_type, NodesAlive...>());

//         constexpr bool other_types_depend_on_this =
//             (depends<NodesAlive, this_type>::call() || ...);

//         if constexpr (other_types_depend_on_this) {
//             // this_type will come up again because it is included on
//             NodesAlive return calc_order_aux_t<NodesAlive...>::template call<
//                 Operations...>();
//         } else {
//             return calc_order_aux_t<Node..., NodesAlive...>::template call<
//                 this_type, Operations...>();
//         }
//     }
// };

} // namespace detail

template <class... Roots> constexpr auto calc_order_t(Roots const &.../* o */) {
    return detail::calc_order_aux_t<Roots...>::template call();
}

} // namespace adhoc2

#endif // ADHOC_TAPE_SIZE_FWD_HPP