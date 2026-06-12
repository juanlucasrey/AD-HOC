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

#ifndef ADHOC_MAPS_HPP
#define ADHOC_MAPS_HPP

#include "map/ankerl_unordered_dense.h"
#include "map/boost_unordered.hpp"
#include <map>
#include <unordered_map>

namespace adhoc {

enum class MapType {
    STD_MAP,
    STD_UNORDERED_MAP,
    ANKERL_UNORDERED_DENSE,
    BOOST_UNORDERED_MAP,
};

template<MapType maptype, class Key, class T>
using map_t = std::conditional_t<maptype == MapType::STD_MAP,
                                 std::map<Key, T>,
                                 std::conditional_t<maptype == MapType::STD_UNORDERED_MAP,
                                                    std::unordered_map<Key, T>,
                                                    std::conditional_t<maptype == MapType::ANKERL_UNORDERED_DENSE,
                                                                       ankerl::unordered_dense::map<Key, T>,
                                                                       boost::unordered_flat_map<Key, T> > > >;

template<MapType maptype>
consteval auto
has_stable_references() -> bool
{
    if constexpr (maptype == MapType::STD_MAP) {
        return true;
    }
    else if constexpr (maptype == MapType::STD_UNORDERED_MAP) {
        return true;
    }
    else if constexpr (maptype == MapType::ANKERL_UNORDERED_DENSE) {
        return false;
    }
    else if constexpr (maptype == MapType::BOOST_UNORDERED_MAP) {
        return false;
    }
}

} // namespace adhoc

#endif // ADHOC_MAPS_HPP
