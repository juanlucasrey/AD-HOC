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

#ifndef BRNG_TOOLS_CIRCULAR_BUFFER
#define BRNG_TOOLS_CIRCULAR_BUFFER

#include "mask.hpp"

#include <array>
#include <bit>

namespace adhoc {

template <class T, std::size_t N> class circular_buffer final {
    static_assert(N > 0);

  public:
    using reference = T &;
    using const_reference = const T &;

    auto data() -> std::array<T, N> & { return this->data_; }
    auto index() -> std::size_t { return this->index_; }
    auto static constexpr size() -> std::size_t { return N; }

    auto operator++() {
        if constexpr (std::has_single_bit(N)) {
            constexpr auto m = mask<std::size_t>(std::bit_width(N) - 1);
            this->index_++;
            this->index_ &= m;
        } else {
            this->index_ = this->index_ == (N - 1) ? 0 : this->index_ + 1;
        }
    }

    auto operator--() {
        if constexpr (std::has_single_bit(N)) {
            constexpr auto m = mask<std::size_t>(std::bit_width(N) - 1);
            this->index_ += m;
            this->index_ &= m;
        } else {
            this->index_ = this->index_ == 0 ? (N - 1) : this->index_ - 1;
        }
    }

    auto operator+=(std::size_t incr) {
        this->index_ += incr;
        if constexpr (std::has_single_bit(N)) {
            constexpr auto m = mask<std::size_t>(std::bit_width(N) - 1);
            this->index_ &= m;
        } else {
            this->index_ %= N;
        }
    }

    auto operator-=(std::size_t incr) {
        this->index_ -= incr;
        if constexpr (std::has_single_bit(N)) {
            constexpr auto m = mask<std::size_t>(std::bit_width(N) - 1);
            this->index_ &= m;
        } else {
            this->index_ %= N;
        }
    }

    template <std::size_t Offset = 0> auto at() -> reference {
        constexpr std::size_t ReducedOffset = Offset % N;
        if constexpr (ReducedOffset == 0) {
            return this->data_[this->index_];
        } else if constexpr (std::has_single_bit(N)) {
            constexpr auto m = mask<std::size_t>(std::bit_width(N) - 1);
            return this->data_[(this->index_ + ReducedOffset) & m];
        } else {
            std::size_t shifted_index = this->index_ + ReducedOffset;
            shifted_index =
                shifted_index >= N ? shifted_index - N : shifted_index;
            return this->data_[shifted_index];
        }
    }

    template <std::size_t Offset = 0> auto at() const -> const_reference {
        constexpr std::size_t ReducedOffset = Offset % N;
        if constexpr (ReducedOffset == 0) {
            return this->data_[this->index_];
        } else if constexpr (std::has_single_bit(N)) {
            constexpr auto m = mask<std::size_t>(std::bit_width(N) - 1);
            return this->data_[(this->index_ + ReducedOffset) & m];
        } else {
            std::size_t shifted_index = this->index_ + ReducedOffset;
            shifted_index =
                shifted_index >= N ? shifted_index - N : shifted_index;
            return this->data_[shifted_index];
        }
    }

    auto operator==(const circular_buffer &rhs) const -> bool {
        if (this == &rhs) {
            return true;
        }

        if (this->index_ == rhs.index_ && this->data_ == rhs.data_) {
            return true;
        }

        auto get = [](const std::array<T, N> &data, std::size_t index,
                      std::size_t offset) -> const T & {
            return data[(index + offset) % N];
        };

        for (std::size_t i = 0; i < N; ++i) {
            if (get(this->data_, this->index_, i) !=
                get(rhs.data_, rhs.index_, i)) {
                return false;
            }
        }

        return true;
    }

    auto operator!=(const circular_buffer &rhs) -> bool {
        return !(this->operator==(rhs));
    }

  private:
    std::array<T, N> data_{};
    std::size_t index_{0};
};

} // namespace adhoc

#endif // BRNG_TOOLS_CIRCULAR_BUFFER
