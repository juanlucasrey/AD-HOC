#ifndef ADHOC3_CALC_TREE_HPP
#define ADHOC3_CALC_TREE_HPP

#include "calc_order.hpp"
#include "leafs.hpp"
#include "tuple_utils.hpp"

#include <tuple>

namespace adhoc3 {

template <class... Roots> class CalcTree {
    using leafs = decltype(detail::leafs_t<Roots...>::template call());
    using ValuesTuple = decltype(std::tuple_cat(
        leafs{},
        detail::calc_order_aux_t<false, true, Roots...>::template call()));

    std::array<double, std::tuple_size<ValuesTuple>{}> m_values{};

  public:
    explicit CalcTree() = default;
    explicit CalcTree(Roots... /* out */) {}
    template <class Derived> auto inline val(Base<Derived> var) const -> double;
    template <class Input> auto inline set(Input var) -> double &;

  private:
    template <constants::ArgType D>
    inline auto val_aux(constants::CD<D> /* node */) -> double {
        return constants::CD<D>::v();
    }

    template <class Node> inline auto val_aux(Node /* node */) -> double {
        constexpr auto idx = get_idx_first2<Node>(ValuesTuple{});
        return this->m_values[idx];
    }

    inline void evaluate_fwd_aux(std::tuple<> /* nodes */) {}

    template <template <class> class Univariate, class Node,
              class... NodesToCalc>
    inline void
    evaluate_fwd_aux(std::tuple<Univariate<Node>, NodesToCalc...> /* nodes
    */) {
        using CurrentNode = Univariate<Node>;
        constexpr auto idx = get_idx_first2<CurrentNode>(ValuesTuple{});
        this->m_values[idx] = CurrentNode::v(this->val_aux(Node{}));
        this->evaluate_fwd_aux(std::tuple<NodesToCalc...>{});
    }

    template <template <class, class> class Bivariate, class Node1, class Node2,
              class... NodesToCalc>
    inline void evaluate_fwd_aux(
        std::tuple<Bivariate<Node1, Node2>, NodesToCalc...> /* nodes */) {
        using CurrentNode = Bivariate<Node1, Node2>;
        constexpr auto idx = get_idx_first2<CurrentNode>(ValuesTuple{});
        this->m_values[idx] =
            CurrentNode::v(this->val_aux(Node1{}), this->val_aux(Node2{}));
        this->evaluate_fwd_aux(std::tuple<NodesToCalc...>{});
    }

  public:
    inline auto evaluate() {
        auto constexpr calcs =
            detail::calc_order_aux_t<false, true, Roots...>::template call();
        this->evaluate_fwd_aux(calcs);
    }
};

template <class... Roots>
template <class Derived>
auto CalcTree<Roots...>::val(Base<Derived> /* in */) const -> double {
    constexpr auto idx = get_idx_first2<Derived>(ValuesTuple{});
    return this->m_values[idx];
}

template <class... Roots>
template <class Input>
auto CalcTree<Roots...>::set(Input /* in */) -> double & {
    static_assert(has_type_tuple<Input, leafs>::value,
                  "only leafs are allowed to be set");
    constexpr auto idx = get_idx_first2<Input>(ValuesTuple{});
    return this->m_values[idx];
}

} // namespace adhoc3

#endif // ADHOC3_CALC_TREE_HPP