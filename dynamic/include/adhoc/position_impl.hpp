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

#ifndef ADHOC_POSITION_IMPL_HPP
#define ADHOC_POSITION_IMPL_HPP

#include <cstddef>

namespace adhoc {

struct PositionImpl {
    std::size_t op_position{ 0 };
    std::size_t id_position{ 0 };
    std::size_t val_position{ 0 };
};

} // namespace adhoc

#endif // ADHOC_POSITION_IMPL_HPP
