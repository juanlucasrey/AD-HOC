/*
 * AD-HOC, Automatic Differentiation for High Order Calculations
 *
 * This file is part of the AD-HOC distribution
 * (https://github.com/juanlucasrey/AD-HOC).
 * Copyright (c) 2026 Juan Lucas Rey
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

#ifndef ADHOC_HASH_HPP
#define ADHOC_HASH_HPP

#include "position_impl.hpp"

#include <functional>

namespace adhoc {

template<class T>
auto
hash_span(T const& v, std::size_t start, std::size_t end) -> std::size_t
{
    std::hash<typename T::value_type> hasher;
    std::size_t seed = 0;
    for (std::size_t i = start; i < end; ++i) {
        seed ^= hasher(v[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
}

template<class TapeDataType>
inline auto
hash(PositionImpl const& pos, TapeDataType const& data) -> std::size_t
{
    auto const h1 = hash_span(data.ids, pos.id_position, data.ids.size());
    auto const h2 = hash_span(data.ops, pos.op_position, data.ops.size());
    auto const h = h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
    return h;
}

} // namespace adhoc

#endif // ADHOC_HASH_HPP
