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

#ifndef ADHOC4_TAPE_BUFFER_SIZE_HPP
#define ADHOC4_TAPE_BUFFER_SIZE_HPP

#include "utils/tuple.hpp"

#include <cstddef>
#include <dependency.hpp>

namespace adhoc4 {

namespace detail {

// template <std::size_t MaxWidth, class Output, class Input>
// constexpr auto tape_buffer_size_aux(Output /* output */, Input /* input */,
//                                     std::tuple<> /* nodes */) -> std::size_t
//                                     {
//     return MaxWidth;
// }

// template <std::size_t MaxWidth, class Output, class Input>
// constexpr auto tape_buffer_size_aux(Output /* output */, Input /* input */,
//                                     std::tuple<> /* nodes */) -> std::size_t
//                                     {
//     return MaxWidth;
// }

template <std::size_t MaxWidth, class Output, class Input, class Nodes>
constexpr auto tape_buffer_size_aux(Output output, Input input, Nodes nodes)
    -> std::size_t {
    constexpr auto current_node = head(nodes);
    if constexpr (is_input(current_node)) {
        return MaxWidth;
    } else {
        constexpr auto next_nodes = tail(nodes);
        return tape_buffer_size_aux<MaxWidth + 1>(output, input, next_nodes);
    }
}

} // namespace detail
template <class Output, class Input, class Nodes>
constexpr auto tape_buffer_size(Output output, Input input, Nodes nodes)
    -> std::size_t {
    return detail::tape_buffer_size_aux<0U>(output, input, nodes);
}

} // namespace adhoc4

#endif // ADHOC4_TAPE_BUFFER_SIZE_HPP
