/*
 * AD-HOC, Automatic Differentiation for High Order Calculations
 *
 * This file is part of the AD-HOC distribution
 * (https://github.com/juanlucasrey/AD-HOC).
 * Copyright (c) 2024 Juan Lucas Rey
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ADHOC3_TAPE_HPP
#define ADHOC3_TAPE_HPP

#include "tape_args_tuple.hpp"
#include "tuple_utils.hpp"

namespace adhoc3 {

template <class... OutputsAndDerivatives> class Tape {
  private:
    using outputs =
        decltype(detail::outputs_t<OutputsAndDerivatives...>::template call());

    using inputs = decltype(Inputs(outputs{}));

    using derivatives_inputs =
        decltype(detail::derivatives_t<
                 OutputsAndDerivatives...>::template call());

    using orders = decltype(Orders(derivatives_inputs{}));

    using derivatives_outputs = decltype(ExpandOutputs(outputs{}, orders{}));

    using derivatives_outputs_and_inputs =
        decltype(std::tuple_cat(derivatives_inputs{}, derivatives_outputs{}));

    std::array<double, std::tuple_size_v<derivatives_outputs_and_inputs>>
        m_derivatives{};

  public:
    explicit Tape(OutputsAndDerivatives... /* out */) {}

    template <class D> auto inline der(D var) const -> double;
    template <class D> auto inline der(D var) -> double &;
    void inline reset_der() {
        std::fill(std::begin(this->m_derivatives),
                  std::end(this->m_derivatives), 0.);
    }
};

template <class... OutputsAndDerivatives>
template <class D>
auto Tape<OutputsAndDerivatives...>::der(D /* in */) const -> double {
    constexpr auto idx = get_idx_first2<D>(derivatives_outputs_and_inputs{});
    static_assert(idx != std::tuple_size_v<derivatives_outputs_and_inputs>,
                  "derivative not on tape");
    return this->m_derivatives[idx];
}

template <class... OutputsAndDerivatives>
template <class D>
auto Tape<OutputsAndDerivatives...>::der(D /* in */) -> double & {
    constexpr auto idx = get_idx_first2<D>(derivatives_outputs_and_inputs{});
    static_assert(idx != std::tuple_size_v<derivatives_outputs_and_inputs>,
                  "derivative not on tape");
    return this->m_derivatives[idx];
}

} // namespace adhoc3

#endif // ADHOC3_TAPE_HPP
