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

#ifndef ADHOC_VECTOR_IDX_HPP
#define ADHOC_VECTOR_IDX_HPP

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace adhoc {

template<std::unsigned_integral T = std::size_t>
class vector_idx {
  private:
    std::vector<std::uint8_t> m_vec8;
    std::vector<std::uint16_t> m_vec16;
    std::vector<std::uint32_t> m_vec32;
    std::vector<std::uint64_t> m_vec64;

    std::size_t m_size{ 0 };

  public:
    using value_type = T;

    struct iterator {
        using iterator_category = std::bidirectional_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = T;

        iterator() = default;

        iterator(std::vector<std::uint8_t>::iterator vec8_it,
                 std::vector<std::uint16_t>::iterator vec16_it,
                 std::vector<std::uint32_t>::iterator vec32_it,
                 std::vector<std::uint64_t>::iterator vec64_it,
                 std::vector<std::uint8_t>::iterator vec8_begin,
                 std::vector<std::uint16_t>::iterator vec16_begin,
                 std::vector<std::uint32_t>::iterator vec32_begin,
                 std::vector<std::uint64_t>::iterator vec64_begin,
                 std::vector<std::uint8_t>::iterator vec8_end,
                 std::vector<std::uint16_t>::iterator vec16_end,
                 std::vector<std::uint32_t>::iterator vec32_end,
                 std::vector<std::uint64_t>::iterator vec64_end)
          : m_vec8_it(vec8_it)
          , m_vec16_it(vec16_it)
          , m_vec32_it(vec32_it)
          , m_vec64_it(vec64_it)
          , m_vec8_begin(vec8_begin)
          , m_vec16_begin(vec16_begin)
          , m_vec32_begin(vec32_begin)
          , m_vec64_begin(vec64_begin)
          , m_vec8_end(vec8_end)
          , m_vec16_end(vec16_end)
          , m_vec32_end(vec32_end)
          , m_vec64_end(vec64_end)
        {
        }
        auto operator*() const -> value_type
        {
            if (m_vec8_it != m_vec8_end) {
                return static_cast<std::size_t>(*m_vec8_it);
            }

            if (m_vec16_it != m_vec16_end) {
                return static_cast<std::size_t>(*m_vec16_it);
            }

            if (m_vec32_it != m_vec32_end) {
                return static_cast<std::size_t>(*m_vec32_it);
            }

            if (m_vec64_it == m_vec64_end) {
                throw;
            }
            return static_cast<std::size_t>(*m_vec64_it);
        }

        auto operator++() -> iterator&
        {
            if (m_vec8_it != m_vec8_end) {
                ++m_vec8_it;
                return *this;
            }

            if (m_vec16_it != m_vec16_end) {
                ++m_vec16_it;
                return *this;
            }

            if (m_vec32_it != m_vec32_end) {
                ++m_vec32_it;
                return *this;
            }

            if (m_vec64_it == m_vec64_end) {
                throw;
            }

            ++m_vec64_it;
            return *this;
        }

        auto operator--() -> iterator&
        {
            if (m_vec64_it != m_vec64_begin) {
                --m_vec64_it;
                return *this;
            }

            if (m_vec32_it != m_vec32_begin) {
                --m_vec32_it;
                return *this;
            }

            if (m_vec16_it != m_vec16_begin) {
                --m_vec16_it;
                return *this;
            }

            if (m_vec8_it == m_vec8_begin) {
                throw;
            }

            --m_vec8_it;
            return *this;
        }

        // creates a copy. should be avoided!
        auto operator++(int) -> iterator
        {
            iterator tmp = *this;
            ++(*this);
            return tmp;
        }
        // creates a copy. should be avoided!
        auto operator--(int) -> iterator
        {
            iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        auto operator==(const iterator& rhs) const -> bool
        {
            return this->m_vec8_it == rhs.m_vec8_it && this->m_vec16_it == rhs.m_vec16_it &&
                   this->m_vec32_it == rhs.m_vec32_it && this->m_vec64_it == rhs.m_vec64_it;
        };
        auto operator!=(const iterator& rhs) const -> bool { return !(this->operator==(rhs)); };

      private:
        std::vector<std::uint8_t>::iterator m_vec8_it;
        std::vector<std::uint16_t>::iterator m_vec16_it;
        std::vector<std::uint32_t>::iterator m_vec32_it;
        std::vector<std::uint64_t>::iterator m_vec64_it;

        std::vector<std::uint8_t>::iterator m_vec8_begin;
        std::vector<std::uint16_t>::iterator m_vec16_begin;
        std::vector<std::uint32_t>::iterator m_vec32_begin;
        std::vector<std::uint64_t>::iterator m_vec64_begin;

        std::vector<std::uint8_t>::iterator m_vec8_end;
        std::vector<std::uint16_t>::iterator m_vec16_end;
        std::vector<std::uint32_t>::iterator m_vec32_end;
        std::vector<std::uint64_t>::iterator m_vec64_end;
    };

    static_assert(std::bidirectional_iterator<iterator>);

    auto operator[](std::size_t pos) const -> value_type
    {
        if (pos < this->m_vec8.size()) {
            return static_cast<std::size_t>(this->m_vec8[pos]);
        }
        pos -= this->m_vec8.size();

        if (pos < this->m_vec16.size()) {
            return static_cast<std::size_t>(this->m_vec16[pos]);
        }
        pos -= this->m_vec16.size();

        if (pos < this->m_vec32.size()) {
            return static_cast<std::size_t>(this->m_vec32[pos]);
        }
        pos -= this->m_vec32.size();

        return static_cast<std::size_t>(this->m_vec64[pos]);
    }

    auto size() const -> std::size_t { return this->m_size; }

    auto begin() -> iterator
    {
        return { m_vec8.begin(),  m_vec16.begin(), m_vec32.begin(), m_vec64.begin(), m_vec8.begin(), m_vec16.begin(),
                 m_vec32.begin(), m_vec64.begin(), m_vec8.end(),    m_vec16.end(),   m_vec32.end(),  m_vec64.end() };
    }
    auto end() -> iterator
    {
        return { m_vec8.end(),    m_vec16.end(),   m_vec32.end(), m_vec64.end(), m_vec8.begin(), m_vec16.begin(),
                 m_vec32.begin(), m_vec64.begin(), m_vec8.end(),  m_vec16.end(), m_vec32.end(),  m_vec64.end() };
    }

    void push_back(value_type val)
    {
        constexpr value_type two_power_8 = 1ULL << 8U;
        constexpr value_type two_power_16 = 1ULL << 16U;
        constexpr value_type two_power_32 = 1ULL << 32U;

        if (this->m_vec64.size() > 0 || val >= two_power_32) {
            this->m_vec64.push_back(static_cast<std::uint64_t>(val));
        }
        else if (this->m_vec32.size() > 0 || val >= two_power_16) {
            this->m_vec32.push_back(static_cast<std::uint32_t>(val));
        }
        else if (this->m_vec16.size() > 0 || val >= two_power_8) {
            this->m_vec16.push_back(static_cast<std::uint16_t>(val));
        }
        else {
            this->m_vec8.push_back(static_cast<std::uint8_t>(val));
        }

        ++this->m_size;
    }

    void resize(std::size_t count)
    {
        if (count > this->m_size) {
            // resize only allowed for reduction of size
            throw;
        }

        std::size_t reduction = this->m_size - count;
        std::size_t reduction_64 = std::min(reduction, this->m_vec64.size());
        this->m_vec64.resize(this->m_vec64.size() - reduction_64);
        reduction -= reduction_64;

        std::size_t reduction_32 = std::min(reduction, this->m_vec32.size());
        this->m_vec32.resize(this->m_vec32.size() - reduction_32);
        reduction -= reduction_32;

        std::size_t reduction_16 = std::min(reduction, this->m_vec16.size());
        this->m_vec16.resize(this->m_vec16.size() - reduction_16);
        reduction -= reduction_16;

        std::size_t reduction_8 = std::min(reduction, this->m_vec8.size());
        this->m_vec8.resize(this->m_vec8.size() - reduction_8);
        reduction -= reduction_8;

        this->m_size = count;
    }

    void clear()
    {
        this->m_vec8.clear();
        this->m_vec16.clear();
        this->m_vec32.clear();
        this->m_vec64.clear();
        this->m_size = 0;
    }

    auto size_of(bool capacity = false) const -> std::size_t
    {
        std::size_t size = 0;
        size += sizeof(std::size_t); // m_size
        size += sizeof(std::uint8_t) * (capacity ? this->m_vec8.capacity() : this->m_vec8.size());
        size += sizeof(std::uint16_t) * (capacity ? this->m_vec16.capacity() : this->m_vec16.size());
        size += sizeof(std::uint32_t) * (capacity ? this->m_vec32.capacity() : this->m_vec32.size());
        size += sizeof(std::uint64_t) * (capacity ? this->m_vec64.capacity() : this->m_vec64.size());
        return size;
    }
};

// template<std::unsigned_integral T = std::size_t>
// class vector_idx2 {
//   private:
//     static constexpr bool use_branching = false;
//     std::vector<std::uint64_t> m_vec;
//     std::size_t m_size{ 0 };

//   public:
//     using value_type = T;

//     struct iterator {
//         using iterator_category = std::bidirectional_iterator_tag;
//         using difference_type = std::ptrdiff_t;
//         using value_type = T;

//         iterator() = default;

//         explicit iterator(std::vector<std::uint64_t>::const_iterator vec_it, std::size_t pos)
//           : m_pos(pos)
//           , m_real_pos(0)
//           , m_offset(0)
//           , m_offset_increase(8)
//           , m_mask(0xFFU)
//           , m_vec_it(vec_it)
//         {
//             constexpr std::size_t two_power_8 = 1ULL << 8U;
//             constexpr std::size_t two_power_16 = 1ULL << 16U;
//             constexpr std::size_t two_power_32 = 1ULL << 32U;

//             if (pos < two_power_8) {
//                 this->m_real_pos = pos / 8;
//                 this->m_offset = (pos % 8) * 8;
//                 this->m_offset_increase = 8;
//                 this->m_mask = 0xFFU;
//             }
//             else if (pos < two_power_16) {
//                 constexpr std::size_t size_previous = two_power_8 / 8;
//                 std::size_t const size_remaining = pos - two_power_8;

//                 this->m_real_pos = size_previous + (size_remaining / 4);
//                 this->m_offset = (size_remaining % 4) * 16;
//                 this->m_offset_increase = 16;
//                 this->m_mask = 0xFFFFU;
//             }
//             else if (pos < two_power_32) {
//                 constexpr std::size_t size_previous = (two_power_8 / 8) + ((two_power_16 - two_power_8) / 4);
//                 std::size_t const size_remaining = pos - two_power_16;

//                 this->m_real_pos = size_previous + (size_remaining / 2);
//                 this->m_offset = (size_remaining % 2) * 32;
//                 this->m_offset_increase = 32;
//                 this->m_mask = 0xFFFFFFFFU;
//             }
//             else {
//                 constexpr std::size_t size_previous =
//                   (two_power_8 / 8) + ((two_power_16 - two_power_8) / 4) + ((two_power_32 - two_power_16) / 2);
//                 std::size_t const size_remaining = pos - two_power_32;
//                 this->m_real_pos = size_previous + size_remaining;
//                 this->m_offset = 0;
//                 this->m_offset_increase = 64;
//                 this->m_mask = 0xFFFFFFFFFFFFFFFFU;
//             }
//         }

//         auto operator*() const -> value_type { return ((*m_vec_it) >> m_offset) & m_mask; }

//         auto operator++() -> iterator&
//         {
//             if constexpr (use_branching) {
//                 constexpr std::size_t two_power_8 = 1ULL << 8U;
//                 constexpr std::size_t two_power_16 = 1ULL << 16U;
//                 constexpr std::size_t two_power_32 = 1ULL << 32U;

//                 this->m_offset += this->m_offset_increase;
//                 this->m_offset %= 64;
//                 if (!static_cast<bool>(this->m_offset)) {
//                     ++this->m_real_pos;
//                     ++this->m_vec_it;
//                 }

//                 ++this->m_pos;

//                 if (this->m_pos == two_power_8) {
//                     this->m_offset_increase = 16;
//                     this->m_mask = 0xFFFFU;
//                 }
//                 else if (this->m_pos == two_power_16) {
//                     this->m_offset_increase = 32;
//                     this->m_mask = 0xFFFFFFFFU;
//                 }
//                 else if (this->m_pos == two_power_32) {
//                     this->m_offset_increase = 64;
//                     this->m_mask = 0xFFFFFFFFFFFFFFFFU;
//                 }

//                 return *this;
//             }
//             else {
//                 this->m_offset += this->m_offset_increase;
//                 this->m_offset %= 64;

//                 // Branchless conditional increment
//                 std::size_t const offset_is_zero = !static_cast<bool>(this->m_offset);
//                 this->m_real_pos += offset_is_zero;
//                 this->m_vec_it += offset_is_zero;

//                 ++this->m_pos;

//                 // Branchless lookup table for offset_increase and mask updates
//                 std::size_t const index =
//                   ((this->m_pos >> 8) != 0) + ((this->m_pos >> 16) != 0) + ((this->m_pos >> 32) != 0);

//                 static constexpr std::array<std::size_t, 4> offsets = { 8U, 16U, 32U, 64U };
//                 static constexpr std::array<std::size_t, 4> masks = {
//                     0xFFU, 0xFFFFU, 0xFFFFFFFFU, 0xFFFFFFFFFFFFFFFFU
//                 };

//                 this->m_offset_increase = offsets[index];
//                 this->m_mask = masks[index];

//                 return *this;
//             }
//         }

//         auto operator--() -> iterator&
//         {
//             if constexpr (use_branching) {
//                 constexpr std::size_t two_power_8 = 1ULL << 8U;
//                 constexpr std::size_t two_power_16 = 1ULL << 16U;
//                 constexpr std::size_t two_power_32 = 1ULL << 32U;

//                 if (this->m_pos == two_power_8) {
//                     this->m_offset_increase = 8;
//                     this->m_mask = 0xFFU;
//                 }
//                 else if (this->m_pos == two_power_16) {
//                     this->m_offset_increase = 16;
//                     this->m_mask = 0xFFFFU;
//                 }
//                 else if (this->m_pos == two_power_32) {
//                     this->m_offset_increase = 32;
//                     this->m_mask = 0xFFFFFFFFU;
//                 }

//                 --this->m_pos;

//                 if (!static_cast<bool>(this->m_offset)) {
//                     --this->m_real_pos;
//                     --this->m_vec_it;
//                     this->m_offset = 64;
//                 }

//                 this->m_offset -= this->m_offset_increase;

//                 return *this;
//             }
//             else {
//                 --this->m_pos;

//                 // Branchless lookup table for offset_increase and mask updates
//                 std::size_t const index =
//                   ((this->m_pos >> 8) != 0) + ((this->m_pos >> 16) != 0) + ((this->m_pos >> 32) != 0);

//                 static constexpr std::array<std::size_t, 4> offsets = { 8U, 16U, 32U, 64U };
//                 static constexpr std::array<std::size_t, 4> masks = {
//                     0xFFU, 0xFFFFU, 0xFFFFFFFFU, 0xFFFFFFFFFFFFFFFFU
//                 };

//                 this->m_offset_increase = offsets[index];
//                 this->m_mask = masks[index];

//                 // Branchless conditional increment
//                 std::size_t const offset_is_zero = !static_cast<bool>(this->m_offset);
//                 this->m_real_pos -= offset_is_zero;
//                 this->m_vec_it -= offset_is_zero;

//                 this->m_offset = offset_is_zero ? 64 : this->m_offset;
//                 this->m_offset -= this->m_offset_increase;

//                 return *this;
//             }
//         }

//         // creates a copy. should be avoided!
//         auto operator++(int) -> iterator
//         {
//             iterator tmp = *this;
//             ++(*this);
//             return tmp;
//         }
//         // creates a copy. should be avoided!
//         auto operator--(int) -> iterator
//         {
//             iterator tmp = *this;
//             ++(*this);
//             return tmp;
//         }

//         auto operator==(const iterator& rhs) const -> bool { return this->m_pos == rhs.m_pos; };
//         auto operator!=(const iterator& rhs) const -> bool { return !(this->operator==(rhs)); };

//       private:
//         std::size_t m_pos;
//         std::size_t m_real_pos;
//         std::size_t m_offset;
//         std::size_t m_offset_increase;
//         std::size_t m_mask;
//         std::vector<std::uint64_t>::const_iterator m_vec_it;
//     };

//     static_assert(std::bidirectional_iterator<iterator>);

//     auto operator[](std::size_t pos) const -> value_type
//     {
//         constexpr std::size_t two_power_8 = 1ULL << 8U;
//         constexpr std::size_t two_power_16 = 1ULL << 16U;
//         constexpr std::size_t two_power_32 = 1ULL << 32U;

//         if (pos < two_power_8) {
//             std::size_t const pos_vec = pos / 8;
//             std::size_t const offset = (pos % 8) * 8;
//             return static_cast<value_type>((this->m_vec[pos_vec] >> offset) & 0xFFU);
//         }

//         if (pos < two_power_16) {
//             constexpr std::size_t size_previous = two_power_8 / 8;
//             std::size_t const size_remaining = pos - two_power_8;

//             std::size_t const pos_vec = size_previous + (size_remaining / 4);
//             std::size_t const offset = (size_remaining % 4) * 16;
//             return static_cast<value_type>((this->m_vec[pos_vec] >> offset) & 0xFFFFU);
//         }

//         if (pos < two_power_32) {
//             constexpr std::size_t size_previous = (two_power_8 / 8) + ((two_power_16 - two_power_8) / 4);
//             std::size_t const size_remaining = pos - two_power_16;

//             std::size_t const pos_vec = size_previous + (size_remaining / 2);
//             std::size_t const offset = (size_remaining % 2) * 32;
//             return static_cast<value_type>((this->m_vec[pos_vec] >> offset) & 0xFFFFFFFFU);
//         }

//         constexpr std::size_t size_previous =
//           (two_power_8 / 8) + ((two_power_16 - two_power_8) / 4) + ((two_power_32 - two_power_16) / 2);
//         std::size_t const size_remaining = pos - two_power_32;
//         std::size_t const pos_vec = size_previous + size_remaining;
//         return static_cast<value_type>(this->m_vec[pos_vec]);
//     }

//     auto size() const -> std::size_t { return this->m_size; }

//     auto begin() -> iterator { return iterator{ this->m_vec.cbegin(), 0 }; }
//     auto end() -> iterator { return iterator{ this->m_vec.cend(), this->m_size }; }

//     void push_back(value_type val)
//     {
//         if (val > this->m_size) {
//             // should not happen
//             throw;
//         }

//         constexpr std::size_t two_power_8 = 1ULL << 8U;
//         constexpr std::size_t two_power_16 = 1ULL << 16U;
//         constexpr std::size_t two_power_32 = 1ULL << 32U;

//         if (this->m_size < two_power_8) {
//             std::size_t const pos = this->m_size / 8;
//             std::size_t const offset = (this->m_size % 8) * 8;
//             if (pos >= this->m_vec.size()) {
//                 this->m_vec.push_back(0);
//             }
//             this->m_vec[pos] |= (static_cast<std::uint64_t>(val) << offset);
//         }
//         else if (this->m_size < two_power_16) {
//             constexpr std::size_t size_previous = two_power_8 / 8;
//             std::size_t const size_remaining = this->m_size - two_power_8;

//             std::size_t const pos = size_previous + (size_remaining / 4);
//             std::size_t const offset = (size_remaining % 4) * 16;
//             if (pos >= this->m_vec.size()) {
//                 this->m_vec.push_back(0);
//             }
//             this->m_vec[pos] |= (static_cast<std::uint64_t>(val) << offset);
//         }
//         else if (this->m_size < two_power_32) {
//             constexpr std::size_t size_previous = (two_power_8 / 8) + ((two_power_16 - two_power_8) / 4);
//             std::size_t const size_remaining = this->m_size - two_power_16;

//             std::size_t const pos = size_previous + (size_remaining / 2);
//             std::size_t const offset = (size_remaining % 2) * 32;
//             if (pos >= this->m_vec.size()) {
//                 this->m_vec.push_back(0);
//             }
//             this->m_vec[pos] |= (static_cast<std::uint64_t>(val) << offset);
//         }
//         else {
//             this->m_vec.push_back(val);
//         }

//         ++this->m_size;
//     }

//     void resize(std::size_t count)
//     {
//         if (count > this->m_size) {
//             // resize only allowed for reduction of size
//             throw;
//         }

//         constexpr std::size_t two_power_8 = 1ULL << 8U;
//         constexpr std::size_t two_power_16 = 1ULL << 16U;
//         constexpr std::size_t two_power_32 = 1ULL << 32U;

//         if (count < two_power_8) {
//             std::size_t const size_remaining = count;
//             std::size_t const size_remaining_div = size_remaining / 8;
//             std::size_t size_remaining_rem = size_remaining % 8;
//             bool const has_remaining = size_remaining_rem > 0;
//             this->m_vec.resize(size_remaining_div + static_cast<std::size_t>(has_remaining));
//             if (has_remaining) {
//                 auto const local_mask = mask<std::uint64_t>(size_remaining_rem * 8);
//                 this->m_vec.back() &= local_mask;
//             }
//         }
//         else if (count < two_power_16) {
//             constexpr std::size_t size_previous = two_power_8 / 8;
//             std::size_t const size_remaining = count - two_power_8;
//             std::size_t const size_remaining_div = size_remaining / 4;
//             std::size_t size_remaining_rem = size_remaining % 4;
//             bool const has_remaining = size_remaining_rem > 0;
//             this->m_vec.resize(size_previous + size_remaining_div + static_cast<std::size_t>(has_remaining));
//             if (has_remaining) {
//                 auto const local_mask = mask<std::uint64_t>(size_remaining_rem * 16);
//                 this->m_vec.back() &= local_mask;
//             }
//         }
//         else if (count < two_power_32) {
//             constexpr std::size_t size_previous = (two_power_8 / 8) + ((two_power_16 - two_power_8) / 4);
//             std::size_t const size_remaining = count - two_power_16;

//             std::size_t const size_remaining_div = size_remaining / 2;
//             std::size_t size_remaining_rem = size_remaining % 2;
//             bool const has_remaining = size_remaining_rem > 0;

//             this->m_vec.resize(size_previous + size_remaining_div + static_cast<std::size_t>(has_remaining));

//             if (has_remaining) {
//                 auto const local_mask = mask<std::uint64_t>(size_remaining_rem * 32);
//                 this->m_vec.back() &= local_mask;
//             }
//         }
//         else {
//             constexpr std::size_t size_previous =
//               (two_power_8 / 8) + ((two_power_16 - two_power_8) / 4) + ((two_power_32 - two_power_16) / 2);
//             std::size_t const real_size = (count - two_power_32) + size_previous;
//             this->m_vec.resize(real_size);
//         }

//         this->m_size = count;
//     }

//     void clear()
//     {
//         this->m_vec.clear();
//         this->m_size = 0;
//     }

//     auto size_of(bool capacity = false) const -> std::size_t
//     {
//         std::size_t size = 0;
//         size += sizeof(std::size_t); // m_size
//         size += sizeof(std::uint64_t) * (capacity ? this->m_vec.capacity() : this->m_vec.size());
//         return size;
//     }
// };

// template<std::unsigned_integral T = std::size_t>
// class vector_idx3 {
//   private:
//     std::vector<std::uint8_t> m_vec8;
//     std::vector<std::uint16_t> m_vec16;
//     std::vector<std::uint32_t> m_vec32;
//     std::vector<std::uint64_t> m_vec64;

//     std::size_t m_size{ 0 };

//   public:
//     using value_type = T;

//     struct iterator {
//         using iterator_category = std::bidirectional_iterator_tag;
//         using difference_type = std::ptrdiff_t;
//         using value_type = T;

//         iterator() = default;

//         iterator(std::vector<std::uint8_t>::iterator vec8_it,
//                  std::vector<std::uint16_t>::iterator vec16_it,
//                  std::vector<std::uint32_t>::iterator vec32_it,
//                  std::vector<std::uint64_t>::iterator vec64_it,
//                  std::vector<std::uint8_t>::iterator vec8_begin,
//                  std::vector<std::uint16_t>::iterator vec16_begin,
//                  std::vector<std::uint32_t>::iterator vec32_begin,
//                  std::vector<std::uint64_t>::iterator vec64_begin,
//                  std::vector<std::uint8_t>::iterator vec8_end,
//                  std::vector<std::uint16_t>::iterator vec16_end,
//                  std::vector<std::uint32_t>::iterator vec32_end,
//                  std::vector<std::uint64_t>::iterator vec64_end)
//           : m_vec8_it(vec8_it)
//           , m_vec16_it(vec16_it)
//           , m_vec32_it(vec32_it)
//           , m_vec64_it(vec64_it)
//           , m_vec8_begin(vec8_begin)
//           , m_vec16_begin(vec16_begin)
//           , m_vec32_begin(vec32_begin)
//           , m_vec64_begin(vec64_begin)
//           , m_vec8_end(vec8_end)
//           , m_vec16_end(vec16_end)
//           , m_vec32_end(vec32_end)
//           , m_vec64_end(vec64_end)
//         {
//         }
//         auto operator*() const -> value_type
//         {
//             if (m_vec8_it != m_vec8_end) {
//                 return static_cast<std::size_t>(*m_vec8_it);
//             }

//             if (m_vec16_it != m_vec16_end) {
//                 return static_cast<std::size_t>(*m_vec16_it);
//             }

//             if (m_vec32_it != m_vec32_end) {
//                 return static_cast<std::size_t>(*m_vec32_it);
//             }

//             if (m_vec64_it == m_vec64_end) {
//                 throw;
//             }
//             return static_cast<std::size_t>(*m_vec64_it);
//         }

//         auto operator++() -> iterator&
//         {
//             if (m_vec8_it != m_vec8_end) {
//                 ++m_vec8_it;
//                 return *this;
//             }

//             if (m_vec16_it != m_vec16_end) {
//                 ++m_vec16_it;
//                 return *this;
//             }

//             if (m_vec32_it != m_vec32_end) {
//                 ++m_vec32_it;
//                 return *this;
//             }

//             if (m_vec64_it == m_vec64_end) {
//                 throw;
//             }

//             ++m_vec64_it;
//             return *this;
//         }

//         auto operator--() -> iterator&
//         {
//             if (m_vec64_it != m_vec64_begin) {
//                 --m_vec64_it;
//                 return *this;
//             }

//             if (m_vec32_it != m_vec32_begin) {
//                 --m_vec32_it;
//                 return *this;
//             }

//             if (m_vec16_it != m_vec16_begin) {
//                 --m_vec16_it;
//                 return *this;
//             }

//             if (m_vec8_it == m_vec8_begin) {
//                 throw;
//             }

//             --m_vec8_it;
//             return *this;
//         }

//         // creates a copy. should be avoided!
//         auto operator++(int) -> iterator
//         {
//             iterator tmp = *this;
//             ++(*this);
//             return tmp;
//         }
//         // creates a copy. should be avoided!
//         auto operator--(int) -> iterator
//         {
//             iterator tmp = *this;
//             ++(*this);
//             return tmp;
//         }

//         auto operator==(const iterator& rhs) const -> bool
//         {
//             return this->m_vec8_it == rhs.m_vec8_it && this->m_vec16_it == rhs.m_vec16_it &&
//                    this->m_vec32_it == rhs.m_vec32_it && this->m_vec64_it == rhs.m_vec64_it;
//         };
//         auto operator!=(const iterator& rhs) const -> bool { return !(this->operator==(rhs)); };

//       private:
//         std::vector<std::uint8_t>::iterator m_vec8_it;
//         std::vector<std::uint16_t>::iterator m_vec16_it;
//         std::vector<std::uint32_t>::iterator m_vec32_it;
//         std::vector<std::uint64_t>::iterator m_vec64_it;

//         std::vector<std::uint8_t>::iterator m_vec8_begin;
//         std::vector<std::uint16_t>::iterator m_vec16_begin;
//         std::vector<std::uint32_t>::iterator m_vec32_begin;
//         std::vector<std::uint64_t>::iterator m_vec64_begin;

//         std::vector<std::uint8_t>::iterator m_vec8_end;
//         std::vector<std::uint16_t>::iterator m_vec16_end;
//         std::vector<std::uint32_t>::iterator m_vec32_end;
//         std::vector<std::uint64_t>::iterator m_vec64_end;
//     };

//     static_assert(std::bidirectional_iterator<iterator>);

//     auto operator[](std::size_t pos) const -> value_type
//     {
//         constexpr std::size_t two_power_8 = 1ULL << 8U;
//         constexpr std::size_t two_power_16 = 1ULL << 16U;
//         constexpr std::size_t two_power_32 = 1ULL << 32U;

//         if (pos < two_power_8) {
//             return static_cast<std::size_t>(this->m_vec8[pos]);
//         }

//         if (pos < two_power_16) {
//             return static_cast<std::size_t>(this->m_vec16[pos - two_power_8]);
//         }

//         if (pos < two_power_32) {
//             return static_cast<std::size_t>(this->m_vec32[pos - two_power_16]);
//         }

//         return static_cast<std::size_t>(this->m_vec64[pos - two_power_32]);
//     }

//     auto size() const -> std::size_t { return this->m_size; }

//     auto begin() -> iterator
//     {
//         return { m_vec8.begin(),  m_vec16.begin(), m_vec32.begin(), m_vec64.begin(), m_vec8.begin(), m_vec16.begin(),
//                  m_vec32.begin(), m_vec64.begin(), m_vec8.end(),    m_vec16.end(),   m_vec32.end(),  m_vec64.end() };
//     }
//     auto end() -> iterator
//     {
//         return { m_vec8.end(),    m_vec16.end(),   m_vec32.end(), m_vec64.end(), m_vec8.begin(), m_vec16.begin(),
//                  m_vec32.begin(), m_vec64.begin(), m_vec8.end(),  m_vec16.end(), m_vec32.end(),  m_vec64.end() };
//     }

//     void push_back(value_type val)
//     {
//         if (val > this->m_size) {
//             // should not happen
//             throw;
//         }

//         constexpr std::size_t two_power_8 = 1ULL << 8U;
//         constexpr std::size_t two_power_16 = 1ULL << 16U;
//         constexpr std::size_t two_power_32 = 1ULL << 32U;

//         if (this->m_size < two_power_8) {
//             m_vec8.push_back(static_cast<std::uint8_t>(val));
//         }
//         else if (this->m_size < two_power_16) {
//             m_vec16.push_back(static_cast<std::uint16_t>(val));
//         }
//         else if (this->m_size < two_power_32) {
//             m_vec32.push_back(static_cast<std::uint32_t>(val));
//         }
//         else {
//             m_vec64.push_back(static_cast<std::uint64_t>(val));
//         }

//         ++this->m_size;
//     }

//     void resize(std::size_t count)
//     {
//         if (count > this->m_size) {
//             // resize only allowed for reduction of size
//             throw;
//         }

//         std::size_t reduction = this->m_size - count;
//         std::size_t reduction_64 = std::min(reduction, this->m_vec64.size());
//         this->m_vec64.resize(this->m_vec64.size() - reduction_64);
//         reduction -= reduction_64;

//         std::size_t reduction_32 = std::min(reduction, this->m_vec32.size());
//         this->m_vec32.resize(this->m_vec32.size() - reduction_32);
//         reduction -= reduction_32;

//         std::size_t reduction_16 = std::min(reduction, this->m_vec16.size());
//         this->m_vec16.resize(this->m_vec16.size() - reduction_16);
//         reduction -= reduction_16;

//         std::size_t reduction_8 = std::min(reduction, this->m_vec8.size());
//         this->m_vec8.resize(this->m_vec8.size() - reduction_8);
//         reduction -= reduction_8;

//         this->m_size = count;
//     }

//     void clear()
//     {
//         this->m_vec8.clear();
//         this->m_vec16.clear();
//         this->m_vec32.clear();
//         this->m_vec64.clear();
//         this->m_size = 0;
//     }

//     auto size_of(bool capacity = false) const -> std::size_t
//     {
//         std::size_t size = 0;
//         size += sizeof(std::size_t); // m_size
//         size += sizeof(std::uint8_t) * (capacity ? this->m_vec8.capacity() : this->m_vec8.size());
//         size += sizeof(std::uint16_t) * (capacity ? this->m_vec16.capacity() : this->m_vec16.size());
//         size += sizeof(std::uint32_t) * (capacity ? this->m_vec32.capacity() : this->m_vec32.size());
//         size += sizeof(std::uint64_t) * (capacity ? this->m_vec64.capacity() : this->m_vec64.size());
//         return size;
//     }
// };

} // namespace adhoc

#endif // ADHOC_VECTOR_IDX_HPP
