#ifndef ADHOC_TAPE_HPP
#define ADHOC_TAPE_HPP

#include "calc_order.hpp"
#include "leafs.hpp"
#include "tape_size_bwd.hpp"
#include "tape_size_fwd.hpp"
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

  public:
    std::array<double, tape_size_bwd<ActiveLeafsAndRootsDerivatives...>()>
        m_buffer{};

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

template <constants::ArgType D, class Tape>
inline auto val(Tape const & /* in */, constants::CD<D> /* node */) -> double {
    return constants::CD<D>::v();
}

template <class Node, class Tape>
inline auto val(Tape const &in, Node /* node */) -> double {
    return in.val(Node{});
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