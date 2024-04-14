#ifndef ADHOC3_TAPE2_HPP
#define ADHOC3_TAPE2_HPP

#include <calc_tree.hpp>
#include <tape_backpropagator.hpp>

#include <differential_operator/expand_tree.hpp>
#include <differential_operator/order_differential_operator.hpp>
#include <differential_operator/select_root_derivatives.hpp>

#include "../tests3/type_name.hpp"

#include <array>
#include <tuple>

namespace adhoc3 {

template <class... OutputsAndDerivatives> class Tape2 {
    //   private:
  public:
    using roots =
        decltype(select_root_derivatives2<OutputsAndDerivatives...>());

    using leaves =
        decltype(select_non_root_derivatives2<OutputsAndDerivatives...>());

    std::array<double, sizeof...(OutputsAndDerivatives)> m_derivatives{};

  public:
    explicit Tape2(OutputsAndDerivatives... /* out */) {
        static_assert(std::tuple_size_v<roots> != 0, "roots derivaties empty!");
        static_assert(std::tuple_size_v<leaves> != 0,
                      "leaves derivaties empty!");
        // more checks to be included here later
    }

    template <class D> auto inline get(D var) const -> double;
    template <class D> auto inline set(D var) -> double &;

    void inline reset_der() {
        std::fill(std::begin(this->m_derivatives),
                  std::end(this->m_derivatives), 0.);
    }

    template <class... Roots> void backpropagate(CalcTree<Roots...> ct);
};

template <class... OutputsAndDerivatives>
template <class D>
auto Tape2<OutputsAndDerivatives...>::get(D /* in */) const -> double {
    static_assert(has_type<D, OutputsAndDerivatives...>(),
                  "derivative not on tape");

    return this
        ->m_derivatives[detail::get_index<D, OutputsAndDerivatives...>::value];
}

template <class... OutputsAndDerivatives>
template <class D>
auto Tape2<OutputsAndDerivatives...>::set(D /* in */) -> double & {
    static_assert(has_type<D, OutputsAndDerivatives...>(),
                  "derivative not on tape");

    return this
        ->m_derivatives[detail::get_index<D, OutputsAndDerivatives...>::value];
}

template <class... OutputsAndDerivatives>
template <class... Roots>
void Tape2<OutputsAndDerivatives...>::backpropagate(CalcTree<Roots...> ct) {
    typename decltype(ct)::ValuesTupleInverse co;
    constexpr auto ordered_derivatives = std::tuple_cat(std::make_tuple(
        order_differential_operator(OutputsAndDerivatives{}, co))...);
    constexpr auto ordered_roots = select_root_derivatives(ordered_derivatives);
    constexpr auto ordered_leaves = order_differential_operators(
        select_non_root_derivatives(ordered_derivatives), co);

    constexpr auto deriv_calc_tree =
        expand_tree(co, ordered_roots, ordered_leaves);

    auto tb = TapeBackpropagator(co, deriv_calc_tree, ordered_leaves);
    tb.backpropagate(ct, ordered_derivatives, this->m_derivatives);
}

} // namespace adhoc3

#endif // ADHOC3_TAPE2_HPP
