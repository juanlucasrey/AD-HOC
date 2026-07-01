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

#ifndef ADHOC_VECTOR_ENUM_HPP
#define ADHOC_VECTOR_ENUM_HPP

#include "mask.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace adhoc {

template<class T, std::size_t NumValues, class StoringType = std::uint64_t>
class vector_enum {
  private:
    static constexpr bool use_branching = false;

    template<class Type>
    static constexpr auto get_max() -> std::size_t
    {
        if constexpr (std::is_enum_v<Type>) {
            return std::numeric_limits<std::underlying_type_t<Type> >::max();
        }
        else {
            return std::numeric_limits<Type>::max();
        }
    }

    static_assert(get_max<T>() >= (NumValues - 1), "return type not large enough");
    static_assert(get_max<StoringType>() >= (NumValues - 1), "storing type not large enough");

    static constexpr std::size_t storing_bits = sizeof(StoringType) * 8;
    static constexpr std::size_t enum_bits = []() {
        std::size_t bits = 0;
        while ((1ULL << bits) < NumValues) {
            ++bits;
        }
        return bits;
    }();
    static constexpr std::size_t items_per_value = storing_bits / enum_bits;

    std::vector<StoringType> m_data;
    std::size_t m_size{ 0 };

  public:
    using value_type = T;

    struct iterator {
        using iterator_category = std::bidirectional_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = T;

        iterator() = default;

        explicit iterator(std::vector<StoringType>::iterator data_it, std::size_t pos)
          : m_data_it(data_it)
          , m_offset((pos % items_per_value) * enum_bits)
        {
            std::size_t const idx = pos / items_per_value;
            this->m_data_it += idx;
        }

        auto operator*() const -> value_type
        {
            auto constexpr local_mask = mask<StoringType>(enum_bits);
            return static_cast<T>((*this->m_data_it >> this->m_offset) & local_mask);
        }

        auto operator++() -> iterator&
        {
            constexpr std::size_t limit_bits = items_per_value * enum_bits;
            this->m_offset += enum_bits;
            if (this->m_offset >= limit_bits) {
                this->m_offset = 0;
                ++this->m_data_it;
            }
            return *this;
        }

        auto operator--() -> iterator&
        {
            constexpr std::size_t limit_bits = items_per_value * enum_bits;
            if (this->m_offset == 0) {
                this->m_offset = limit_bits;
                --this->m_data_it;
            }

            this->m_offset -= enum_bits;

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
            return this->m_data_it == rhs.m_data_it && this->m_offset == rhs.m_offset;
        };
        auto operator!=(const iterator& rhs) const -> bool { return !(this->operator==(rhs)); };

      private:
        std::vector<StoringType>::iterator m_data_it;
        std::size_t m_offset{ 0 };
    };

    static_assert(std::bidirectional_iterator<iterator>);

    auto operator[](std::size_t pos) const -> T
    {
        std::size_t const idx = pos / items_per_value;
        std::size_t const offset = (pos % items_per_value) * enum_bits;
        auto constexpr local_mask = mask<StoringType>(enum_bits);
        return static_cast<T>((this->m_data[idx] >> offset) & local_mask);
    }

    auto size() const -> std::size_t { return this->m_size; }
    auto begin() -> iterator { return iterator{ this->m_data.begin(), 0 }; }
    auto end() -> iterator { return iterator{ this->m_data.begin(), this->m_size }; }

    void push_back(value_type val)
    {
        std::size_t const idx = this->m_size / items_per_value;
        std::size_t const offset = (this->m_size % items_per_value) * enum_bits;

        if constexpr (use_branching) {
            if (offset == 0) {
                this->m_data.push_back(static_cast<StoringType>(val));
            }
            else {
                this->m_data[idx] |= static_cast<StoringType>(val) << offset;
            }
        }
        else {
            this->m_data.resize(idx + 1);
            this->m_data[idx] |= static_cast<StoringType>(val) << offset;
        }

        ++this->m_size;
    }

    void resize(std::size_t count)
    {
        std::size_t const idx = count / items_per_value;
        std::size_t const rem = count % items_per_value;
        bool const has_remaining = rem > 0;
        this->m_data.resize(idx + static_cast<std::size_t>(has_remaining));
        if (has_remaining) {
            auto const local_mask = mask<StoringType>(rem * enum_bits);
            this->m_data.back() &= local_mask;
        }

        this->m_size = count;
    }

    void clear()
    {
        this->m_data.clear();
        this->m_size = 0;
    }

    auto size_of(bool capacity = false) const -> std::size_t
    {
        std::size_t size = 0;
        size += sizeof(std::size_t); // m_size
        size += sizeof(StoringType) * (capacity ? this->m_data.capacity() : this->m_data.size());
        return size;
    }
};

template<class T, std::size_t NumValues, class StoringType = std::uint64_t>
class vector_enum2 {
  private:
    static constexpr bool use_branching = false;

    static_assert(NumValues > 1, "NumValues must be greater than 1");
    static_assert(std::is_integral_v<StoringType>&& std::is_unsigned_v<StoringType>,
                  "StoringType must be an unsigned integer type");

    static constexpr std::size_t items_per_value = []() {
        std::size_t count = 0;
        std::size_t encoded_max = 0;

        while (encoded_max <= (std::numeric_limits<StoringType>::max() - (NumValues - 1)) / NumValues) {
            encoded_max = (encoded_max * NumValues) + (NumValues - 1);
            ++count;
        }

        return count;
    }();

    static constexpr auto uint_pow(std::size_t b, std::size_t e) -> std::size_t
    {
        return (e == 0) ? 1 : b * uint_pow(b, e - 1);
    }

    static constexpr std::array<StoringType, items_per_value> precomputed_powers = []() {
        std::array<StoringType, items_per_value> arr{};
        for (std::size_t i = 0; i < items_per_value; ++i) {
            arr[i] = uint_pow(static_cast<StoringType>(NumValues), i);
        }
        return arr;
    }();

    std::vector<StoringType> m_data;
    std::size_t m_size{ 0 };

  public:
    using value_type = T;

    struct iterator {
        using iterator_category = std::bidirectional_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = T;

        iterator() = default;

        explicit iterator(std::vector<StoringType>::iterator data_it, std::size_t pos)
          : m_data_it(data_it)
          , m_rem(pos % items_per_value)
        {
            std::size_t const idx = pos / items_per_value;
            this->m_data_it += idx;
        }

        auto operator*() const -> value_type
        {
            return static_cast<value_type>(((*m_data_it) / precomputed_powers[m_rem]) %
                                           static_cast<StoringType>(NumValues));
        }

        auto operator++() -> iterator&
        {
            ++this->m_rem;
            if (this->m_rem == items_per_value) {
                this->m_rem = 0;
                ++this->m_data_it;
            }
            return *this;
        }

        auto operator--() -> iterator&
        {
            if (this->m_rem == 0) {
                this->m_rem = items_per_value;
                --this->m_data_it;
            }
            --this->m_rem;

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
            return this->m_data_it == rhs.m_data_it && this->m_rem == rhs.m_rem;
        };
        auto operator!=(const iterator& rhs) const -> bool { return !(this->operator==(rhs)); };

      private:
        std::vector<StoringType>::iterator m_data_it;
        std::size_t m_rem{ 0 };
    };

    static_assert(std::bidirectional_iterator<iterator>);

    auto operator[](std::size_t pos) const -> T
    {
        std::size_t const idx = pos / items_per_value;
        std::size_t const rem = pos % items_per_value;
        return static_cast<T>((this->m_data[idx] / precomputed_powers[rem]) % static_cast<StoringType>(NumValues));
    }

    auto size() const -> std::size_t { return this->m_size; }
    auto begin() -> iterator { return iterator{ this->m_data.begin(), 0 }; }
    auto end() -> iterator { return iterator{ this->m_data.begin(), this->m_size }; }

    void push_back(value_type val)
    {
        std::size_t const idx = this->m_size / items_per_value;
        std::size_t const rem = this->m_size % items_per_value;
        if constexpr (use_branching) {
            if (rem == 0) {
                this->m_data.push_back(static_cast<StoringType>(val));
            }
            else {
                this->m_data[idx] += static_cast<StoringType>(val) * precomputed_powers[rem];
            }
        }
        else {
            this->m_data.resize(idx + 1);
            this->m_data[idx] += static_cast<StoringType>(val) * precomputed_powers[rem];
        }

        ++this->m_size;
    }

    void resize(std::size_t count)
    {
        std::size_t const idx = count / items_per_value;
        std::size_t const rem = count % items_per_value;
        bool const has_remaining = rem > 0;
        this->m_data.resize(idx + static_cast<std::size_t>(has_remaining));
        if (has_remaining) {
            auto const local_mask = precomputed_powers[rem];
            this->m_data.back() %= local_mask;
        }

        this->m_size = count;
    }

    void clear()
    {
        this->m_data.clear();
        this->m_size = 0;
    }

    auto size_of(bool capacity = false) const -> std::size_t
    {
        std::size_t size = 0;
        size += sizeof(std::size_t); // m_size
        size += sizeof(StoringType) * (capacity ? this->m_data.capacity() : this->m_data.size());
        return size;
    }
};

} // namespace adhoc

#endif // ADHOC_VECTOR_ENUM_HPP
