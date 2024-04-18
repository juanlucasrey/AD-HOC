#ifndef ADHOC3_TAPE2_HPP
#define ADHOC3_TAPE2_HPP

#include <calc_tree.hpp>
#include <tape_backpropagator.hpp>

#include <differential_operator/expand_tree.hpp>
#include <differential_operator/order.hpp>
#include <differential_operator/select_root_derivatives.hpp>
#include <differential_operator/sort.hpp>

#include "../tests3/type_name.hpp"

#include <array>
#include <tuple>
#include <utility>

namespace adhoc3 {

template <class... Ids, std::size_t... Orders, std::size_t... Powers>
constexpr auto
multiplicity(std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...> /* in */)
    -> std::size_t {
    return MultinomialCoeff2(std::index_sequence<Powers...>{});
}

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
    template <class D> auto inline get_d(D var) const -> double;
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
auto Tape2<OutputsAndDerivatives...>::get_d(D in) const -> double {
    static_assert(has_type<D, OutputsAndDerivatives...>(),
                  "derivative not on tape");

    constexpr auto mult = multiplicity(in);
    return this->m_derivatives
               [detail::get_index<D, OutputsAndDerivatives...>::value] /
           mult;
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
    using NodesValue = decltype(ct)::ValuesTupleInverse;
    constexpr auto nodes_value = NodesValue{};
    constexpr auto ordered_derivatives = std::tuple_cat(std::make_tuple(
        order_differential_operator(OutputsAndDerivatives{}, nodes_value))...);
    constexpr auto ordered_roots = select_root_derivatives(ordered_derivatives);
    constexpr auto output_derivatives_ordered = order_differential_operators(
        select_non_root_derivatives(ordered_derivatives), nodes_value);
    constexpr auto nodes_derivative =
        expand_tree(nodes_value, ordered_roots, output_derivatives_ordered);
    constexpr auto buffer_size = tape_buffer_size(nodes_value, nodes_derivative,
                                                  output_derivatives_ordered);

    std::array<double, buffer_size> buffer{};

    constexpr auto order = max_orders(ordered_derivatives);
    std::array<double, order> univariate_buffer{};

    detail::backpropagate_aux(
        nodes_derivative, ct, ordered_derivatives, this->m_derivatives,
        std::tuple_cat(std::array<detail::available_t, buffer_size>{}), buffer,
        detail::available_t{}, univariate_buffer);
}

} // namespace adhoc3

#endif // ADHOC3_TAPE2_HPP
