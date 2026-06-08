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

#ifndef ADHOC_BUFFER_T2_HPP
#define ADHOC_BUFFER_T2_HPP

#include <optional>
#include <span>
#include <type_traits>
#include <vector>

namespace adhoc {

template<class T = std::nullopt_t, bool Vectorised = false>
class buffer_t2 {
  private:
    std::size_t m_size{ 0 };
    std::size_t m_allocated_size{ 0 };
    std::size_t m_num_lanes{ 1 };
    std::vector<std::size_t> free_positions;
    std::vector<T> m_data;

  public:
    buffer_t2() = default;
    explicit buffer_t2(std::size_t lanes)
      : m_num_lanes(lanes)
    {

        if constexpr (!Vectorised) {
            if (this->m_num_lanes != 1) {
                // This backpropagator is not designed for multiple lanes
                throw;
            }
        }
    }

    template<bool Reset = false, bool Allocate = true>
    auto get_new_loc() -> std::size_t
    {
        if (this->free_positions.empty()) {
            std::size_t pos = this->m_size;
            ++this->m_size;
            if constexpr (Allocate) {
                this->m_data.resize(this->m_size * this->m_num_lanes);
                this->m_allocated_size = this->m_size;
            }
            return pos;
        }

        std::size_t pos = this->free_positions.back();
        this->free_positions.pop_back();

        if constexpr (Reset) {
            if constexpr (Vectorised) {
                double* dest = &this->m_data[pos * this->m_num_lanes];
#pragma omp simd
                for (std::size_t i = 0; i < this->m_num_lanes; ++i) {
                    dest[i] = 0.;
                }
            }
            else {
                this->m_data[pos] = 0.;
            }
        }
        return pos;
    }

    auto free_loc(std::size_t pos) -> void { this->free_positions.push_back(pos); }

    using buffer_result_t = std::conditional_t<Vectorised, std::span<T>, T&>;

    auto operator[](std::size_t pos) -> buffer_result_t
    {
        if constexpr (Vectorised) {
            return std::span<double>(&this->m_data[pos * this->m_num_lanes], this->m_num_lanes);
        }
        else {
            return this->m_data[pos * this->m_num_lanes];
        }
    }

    using buffer_result_const_t = std::conditional_t<Vectorised, std::span<T const>, T>;

    auto operator[](std::size_t pos) const -> const buffer_result_const_t
    {
        if constexpr (Vectorised) {
            return std::span<double const>(&this->m_data[pos * this->m_num_lanes], this->m_num_lanes);
        }
        else {
            return this->m_data[pos * this->m_num_lanes];
        }
    }

    auto size_of(bool capacity = false) const -> std::size_t
    {
        std::size_t size = 0;
        size += 2 * sizeof(std::size_t); // m_num_lanes, size
        size += sizeof(std::size_t) * (capacity ? this->free_positions.capacity() : this->free_positions.size());
        size += sizeof(double) * (capacity ? this->m_data.capacity() : this->m_data.size());

        return size;
    }

    auto size() const -> std::size_t { return this->m_size; }

    void reserve(std::size_t reserve_size) { this->m_data.reserve(this->m_num_lanes * reserve_size); }

    auto zero() { std::fill(this->m_data.begin(), this->m_data.end(), 0.0); }

    auto lanes() const -> std::size_t { return this->m_num_lanes; }
};

} // namespace adhoc

#endif // ADHOC_BUFFER_T_HPP
