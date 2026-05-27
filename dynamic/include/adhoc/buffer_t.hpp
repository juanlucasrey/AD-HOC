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

#ifndef ADHOC_BUFFER_T_HPP
#define ADHOC_BUFFER_T_HPP

#include <type_traits>
#include <vector>

namespace adhoc {

template<class T = std::nullopt_t, bool CheckBounds = false>
class buffer_t {
  private:
    std::size_t m_size{ 0 };
    std::size_t m_num_lanes{ 1 };
    std::vector<std::size_t> free_positions;

    struct empty_t {};
    std::conditional_t<std::is_same_v<T, std::nullopt_t>, empty_t, std::vector<T> > m_data;

  public:
    buffer_t() = default;
    explicit buffer_t(std::size_t lanes)
      : m_num_lanes(lanes)
    {
    }

    template<bool Allocate = true>
    auto get_new_loc() -> std::size_t
    {
        if (this->free_positions.empty()) {
            std::size_t pos = this->m_size;
            ++this->m_size;
            if constexpr (!std::is_same_v<T, std::nullopt_t> && Allocate) {
                this->m_data.resize(this->m_size * this->m_num_lanes);
            }
            return pos;
        }

        std::size_t pos = this->free_positions.back();
        this->free_positions.pop_back();
        return pos;
    }

    auto free_loc(std::size_t pos) -> void
    {
        if constexpr (CheckBounds) {
            if (pos >= this->m_size) {
                // should not happen
                throw;
            }
        }
        this->free_positions.push_back(pos);
    }

    auto operator[](std::size_t pos) -> T&
    {
        static_assert(!std::is_same_v<T, empty_t>, "buffer_t does not hold values");
        if constexpr (CheckBounds) {
            if (pos >= this->m_size) {
                // should not happen
                throw;
            }
        }
        return this->m_data[pos * this->m_num_lanes];
    }

    auto data() -> std::vector<T>&
    {
        static_assert(!std::is_same_v<T, empty_t>, "buffer_t does not hold values");
        return this->m_data;
    }

    auto data() const -> const std::vector<T>&
    {
        static_assert(!std::is_same_v<T, empty_t>, "buffer_t does not hold values");
        return this->m_data;
    }

    auto size_of(bool capacity = false) const -> std::size_t
    {
        std::size_t size = 0;
        size += 2 * sizeof(std::size_t); // m_num_lanes, size
        size += sizeof(std::size_t) * (capacity ? this->free_positions.capacity() : this->free_positions.size());
        if constexpr (!std::is_same_v<T, empty_t>) {
            size += sizeof(T) * (capacity ? this->m_data.capacity() : this->m_data.size());
        }
        return size;
    }

    auto resize(std::size_t new_size) -> void
    {
        if constexpr (!std::is_same_v<T, empty_t>) {
            this->m_data.resize(new_size * this->m_num_lanes);
        }
        this->m_size = new_size;
    }

    auto size() const -> std::size_t { return this->m_size; }

    auto allocate() { this->m_data.resize(this->m_size * this->m_num_lanes); }
};

} // namespace adhoc

#endif // ADHOC_BUFFER_T_HPP
