#ifndef ADHOC_LEAFS_HPP
#define ADHOC_LEAFS_HPP

#include "adhoc.hpp"
#include "constants_type.hpp"
// #include "dependency.hpp"
#include "utils.hpp"

#include <tuple>

namespace adhoc {

namespace detail {

template <typename... NodesAlive> struct leafs_t {
    template <typename... Leafs> constexpr static auto call() noexcept;
};

template <> struct leafs_t<> {
    template <typename... Leafs> constexpr static auto call() noexcept {
        return std::tuple<Leafs...>{};
    }
};

template <constants::ArgType D, typename... NodesAlive>
struct leafs_t<constants::CD<D>, NodesAlive...> {
    template <typename... Leafs> constexpr static auto call() noexcept {
        // we don't add anything because a constant is not a leaf
        return leafs_t<NodesAlive...>::template call<Leafs...>();
    }
};

template <std::size_t N, typename... NodesAlive>
struct leafs_t<double_t<N>, NodesAlive...> {
    template <typename... Leafs> constexpr static auto call() noexcept {
        using this_type = double_t<N>;

        constexpr bool already_included = has_type<this_type, Leafs...>();

        if constexpr (already_included) {
            // It's already been added, so we don't add the leaf
            return leafs_t<NodesAlive...>::template call<Leafs...>();
        } else {
            // we add the leaf
            return leafs_t<NodesAlive...>::template call<Leafs..., this_type>();
        }
    }
};

template <template <class...> class Xvariate, class... Node,
          typename... NodesAlive>
struct leafs_t<Xvariate<Node...>, NodesAlive...> {
    template <typename... Leafs> constexpr static auto call() noexcept {
        return leafs_t<Node..., NodesAlive...>::template call<Leafs...>();
    }
};

} // namespace detail

template <class... Roots> constexpr auto Leafs(Roots const &.../* o */) {
    return detail::leafs_t<Roots...>::template call();
}

} // namespace adhoc

#endif // ADHOC_LEAFS_HPP