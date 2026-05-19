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
#include "tape_data.hpp"

#include <functional>
#include <span>

namespace adhoc {

template<typename T>
auto
hash_span(std::span<const T> v) -> std::size_t
{
    std::hash<T> hasher;
    std::size_t seed = 0;
    for (const auto& x : v) {
        seed ^= hasher(x) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
}

inline auto
hash(PositionImpl const& pos, TapeData const& data) -> std::size_t
{
    std::span<std::size_t const> ids_span(data.ids.data() + pos.id_position, data.ids.size() - pos.id_position);
    std::span<OpCode const> ops_span(data.ops.data() + pos.op_position, data.ops.size() - pos.op_position);
    auto const h1 = hash_span(ids_span);
    auto const h2 = hash_span(ops_span);
    auto const h = h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
    return h;
}

} // namespace adhoc

#endif // ADHOC_HASH_HPP
