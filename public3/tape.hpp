/*
 * AD-HOC, Automatic Differentiation for High Order Calculations
 *
 * This file is part of the AD-HOC distribution
 * (https://github.com/juanlucasrey/AD-HOC).
 * Copyright (c) 2024 Juan Lucas Rey
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
