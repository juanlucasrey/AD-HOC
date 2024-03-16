#ifndef ADHOC3_TAPE2_HPP
#define ADHOC3_TAPE2_HPP

#include "calc_tree.hpp"

#include <array>

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
void Tape2<OutputsAndDerivatives...>::backpropagate(
    CalcTree<Roots...> /* ct */) {

    // exact buffer size to be determined at compile time
    std::array<double, 10> buffer{};
}

} // namespace adhoc3

#endif // ADHOC3_TAPE2_HPP
