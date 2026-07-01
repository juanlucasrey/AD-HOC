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

#ifndef ADHOC_VECTOR_SIZE_OF_HPP
#define ADHOC_VECTOR_SIZE_OF_HPP

#include "vector_enum.hpp"
#include "vector_idx.hpp"

#include <cstddef>
#include <vector>

namespace adhoc {

template<class T>
auto
vector_size_of(std::vector<T> const& vec, bool capacity = false) -> std::size_t
{
    std::size_t size = 0;
    size += sizeof(T) * (capacity ? vec.capacity() : vec.size());
    return size;
}

template<class T>
auto
vector_size_of(vector_idx<T> const& vec, bool capacity) -> std::size_t
{
    return vec.size_of(capacity);
}

template<class T, std::size_t Bits, class StoringType>
auto
vector_size_of(vector_enum<T, Bits, StoringType> const& vec, bool capacity) -> std::size_t
{
    return vec.size_of(capacity);
}

template<class T, std::size_t NumValues, class StoringType>
auto
vector_size_of(vector_enum2<T, NumValues, StoringType> const& vec, bool capacity) -> std::size_t
{
    return vec.size_of(capacity);
}

} // namespace adhoc

#endif // ADHOC_VECTOR_SIZE_OF_HPP
