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

#ifndef BMC_TOOLS_MDSPAN
#define BMC_TOOLS_MDSPAN

#include "array_operators.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <utility>

namespace adhoc {

namespace detail {

// Multicombination (stars and bars): C(n+k-1, k)
// This computes the binomial coefficient (n+k-1 choose k)
// Represents the number of ways to choose k items from n types with repetition
// allowed Also known as: multiset coefficient, combinations with repetition
// Formula: C(n+k-1, k) = (n+k-1)! / (k! * (n-1)!)
constexpr auto
multicombination(std::size_t n, std::size_t k) noexcept -> std::size_t
{
    if (k == 0 || n == 0) {
        return (k == 0) ? 1 : 0;
    }

    // Optimization: C(n+k-1, k) = C(n+k-1, n-1), choose the smaller value
    std::size_t numerator_count = k;
    if (n - 1 < k) {
        numerator_count = n - 1;
    }

    // Calculate C(n+k-1, numerator_count) iteratively to avoid overflow
    // We compute: (n+k-1) * (n+k-2) * ... * (n) / (1 * 2 * ... *
    // numerator_count)
    std::size_t result = 1;
    for (std::size_t i = 0; i < numerator_count; ++i) {
        result *= (n + k - 1 - i);
        result /= (i + 1);
    }

    return result;
}

// Variadic version: convert sorted multi-index to linear index for symmetric
// tensor Works for any number of dimensions Indices must be sorted: idx[0] <=
// idx[1] <= ... <= idx[N-1]
template<typename... Args>
auto
multiindex_to_linear(std::size_t dim, Args... indices) -> std::size_t
{
    constexpr std::size_t N = sizeof...(Args);
    std::array<std::size_t, N> idx_array{ static_cast<std::size_t>(indices)... };

    // Ensure indices are sorted
    std::sort(idx_array.begin(), idx_array.end());

    std::size_t linear_idx = 0;

    // For each position in the sorted multi-index
    for (std::size_t pos = 0; pos < N; ++pos) {
        std::size_t curr_val = idx_array[pos];
        std::size_t prev_val = (pos == 0) ? 0 : idx_array[pos - 1];

        // Count all combinations with smaller values at this position
        for (std::size_t val = prev_val; val < curr_val; ++val) {
            // Remaining dimensions after fixing current position
            std::size_t remaining_dims = N - pos - 1;
            // Available range for remaining indices [val, dim-1]
            std::size_t remaining_range = dim - val;

            // Add number of combinations: C(remaining_range + remaining_dims -
            // 1, remaining_dims)
            if (remaining_dims > 0) {
                linear_idx += multicombination(remaining_range, remaining_dims);
            }
            else {
                linear_idx += 1; // Last dimension, just count
            }
        }
    }

    return linear_idx;
}

// Convert linear index back to sorted multi-index for any dimension
template<std::size_t N>
auto
linear_to_multiindex(std::size_t dim, std::size_t linear_idx) -> std::array<std::size_t, N>
{
    std::array<std::size_t, N> result{};
    std::size_t remaining = linear_idx;
    std::size_t prev_val = 0;

    for (std::size_t pos = 0; pos < N; ++pos) {
        // Find the value at this position
        for (std::size_t val = prev_val; val < dim; ++val) {
            std::size_t remaining_dims = N - pos - 1;
            std::size_t remaining_range = dim - val;

            std::size_t block_size = (remaining_dims > 0) ? multicombination(remaining_range, remaining_dims) : 1;

            if (remaining < block_size) {
                result[pos] = val;
                prev_val = val;
                break;
            }
            remaining -= block_size;
        }
    }

    return result;
}

} // namespace detail

template<class T, std::size_t S = 1>
class mdspan {
  public:
    template<typename... Args>
    explicit mdspan(T* data = nullptr, Args... dims)
      : m_data(data)
      , m_dims({ static_cast<std::size_t>(dims)... })
    {
        static_assert(S > 0, "number of dimensions cannot be 0");
        static_assert(sizeof...(Args) <= S, "number of arguments has to be at most number of dimensions");

        this->sync_offset();
    }

    explicit mdspan(T* data, std::array<std::size_t, S> const& dims)
      : m_data(data)
      , m_dims(dims)
    {
        this->sync_offset();
    }

    template<class T2>
    explicit mdspan(mdspan<T2, S> const& other)
      : m_data(other.m_data)
      , m_dims(other.m_dims)
      , m_offsets(other.m_offsets)
    {
        static_assert(std::is_same_v<std::remove_cv_t<T>, std::remove_cv_t<T2> >, "Invalid span type.");
    }

    // Constructor from container (only for 1D span_tensor)
    // Requires: S == 1, and Cont has .data() and .size() members
    template<class Cont>
    explicit mdspan(Cont& cont) requires(S == 1 &&
                                         requires(Cont & container) {
                                             { container.data() }->std::convertible_to<T*>;
                                             { container.size() }->std::convertible_to<std::size_t>;
                                         })
      : mdspan(cont.data(), cont.size())
    {
    }

    template<typename... Args>
    inline auto operator()(Args... indices) const -> T& requires(sizeof...(Args) == S)
    {
        return this->m_data[this->index_calculator(static_cast<std::size_t>(indices)...)];
    }

    template<typename... Args>
    inline auto operator()(Args... indices) const requires(sizeof...(Args) < S)
    {
        constexpr std::size_t new_dim = S - sizeof...(Args);
        return mdspan<T, new_dim>(&this->m_data[this->index_calculator(static_cast<std::size_t>(indices)...)],
                                  tail<new_dim>(this->m_dims),
                                  tail<new_dim>(this->m_offsets));
    }

    // note: operator[] will only work with 1 argument
    auto operator[](std::size_t index) const -> T& requires(S == 1) { return this->operator()(index); }

    inline auto operator[](std::size_t index) const -> mdspan<T, S - 1> requires(1 < S)
    {
        return this->operator()(index);
    }

    // at() operator is same as operator() but with boundary checks.
    template<typename... Args>
    auto at(Args... indices) const
    {

        static_assert(sizeof...(Args) <= S, "number of arguments has to be at most number of dimensions");

        // Lambda for boundary checking - ensures all indices are within bounds
        // using fold expressions and index sequences for modern C++ style
        auto is_in_bounds = [this]<std::size_t... Idx>(std::index_sequence<Idx...>, auto... idx_vals) constexpr {
            std::array<std::size_t, sizeof...(Args)> idx_array{ static_cast<std::size_t>(idx_vals)... };
            return ((this->m_dims[Idx] > idx_array[Idx]) && ...);
        };

        if (!is_in_bounds(std::make_index_sequence<sizeof...(Args)>{}, indices...)) {
            throw std::out_of_range("index too large");
        }

        return this->operator()(indices...);
    }

    auto data() const { return this->m_data; }
    auto size() const { return this->m_offsets[0]; }
    auto cbegin() const { return this->m_data; }
    auto cend() const { return this->m_data + this->size(); }
    auto begin() const { return this->m_data; }
    auto end() const { return this->m_data + this->size(); }

    auto dims() const -> std::array<std::size_t, S> const& { return m_dims; }

    template<unsigned int D>
    auto dim() const -> std::size_t
    {
        static_assert(D < S, "Multi Dimension array does not contain that many dimensions");
        return this->m_dims[D];
    }

    template<typename... Args>
    static auto storage_size(Args... dims) -> std::size_t
    {
        static_assert(sizeof...(Args) == S, "number of arguments has to be number of dimensions");
        return (static_cast<std::size_t>(dims) * ...);
    }

  private:
    T* m_data;
    std::array<std::size_t, S> m_dims;
    std::array<std::size_t, S> m_offsets;

    void sync_offset()
    {
        // Lambda to recursively initialize offsets using compile-time recursion
        // Using Y-combinator pattern for recursive lambda
        auto init_offset = [this]<std::size_t C>(auto&& self, std::size_t previous = 1) -> void {
            std::get<C>(this->m_offsets) = std::get<C>(this->m_dims) * previous;
            std::size_t new_previous = std::get<C>(this->m_offsets);
            if constexpr (C != 0) {
                self.template operator()<C - 1>(self, new_previous);
            }
        };

        init_offset.template operator()<S - 1>(init_offset);
    }

    template<unsigned int C = 0>
    inline auto index_calculator(std::size_t first) const -> std::size_t
    {
        if constexpr ((C + 1) >= S) {
            return first;
        }
        else {
            return first * std::get<C + 1>(this->m_offsets);
        }
    }

    template<unsigned int C = 0, typename... Args>
    inline auto index_calculator(std::size_t first, Args... args) const -> std::size_t
    {
        return (first * std::get<C + 1>(this->m_offsets)) + index_calculator<C + 1>(args...);
    }

    // this friend declaration is to that the last costructor can be called by
    // the operators. partial template specializations are not allowed on friend
    // declarations, so this the only choice
    template<class T2, std::size_t S2>
    friend class mdspan;

    mdspan(T* data, std::array<std::size_t, S> const& dims, std::array<std::size_t, S> const& offsets)
      : m_data(data)
      , m_dims{ dims }
      , m_offsets{ offsets }
    {
    }
};

// Symmetric mdspan - assumes symmetry across all indices
// For 2D: m(i,j) = m(j,i), stores only upper triangular
// For 3D: m(i,j,k) = m(j,k,i) = m(k,i,j) = ... (all permutations equal)
template<class T, std::size_t S = 2>
class mdspan_sym {
  public:
    // Constructor with dimensions - all dimensions must be equal
    template<typename... Args>
    explicit mdspan_sym(T* data = nullptr, std::size_t dim = 0)
      : m_data(data)
      , m_dim(dim)
    {
        static_assert(S >= 1, "mdspan_sym requires at least 1 dimension");
    }

    // Access operator - automatically handles symmetry
    template<typename... Args>
    auto operator()(Args... indices) const -> T& requires(sizeof...(Args) == S)
    {
        std::array<std::size_t, S> idx_array{ static_cast<std::size_t>(indices)... };

        // Ensure indices are sorted
        std::sort(idx_array.begin(), idx_array.end());

        std::size_t linear_idx = 0;

        // For each position in the sorted multi-index
        for (std::size_t pos = 0; pos < S; ++pos) {
            std::size_t curr_val = idx_array[pos];
            std::size_t prev_val = (pos == 0) ? 0 : idx_array[pos - 1];

            // Count all combinations with smaller values at this position
            for (std::size_t val = prev_val; val < curr_val; ++val) {
                // Remaining dimensions after fixing current position
                std::size_t remaining_dims = S - pos - 1;
                // Available range for remaining indices [val, dim-1]
                std::size_t remaining_range = this->m_dim - val;

                // Add number of combinations: C(remaining_range +
                // remaining_dims - 1, remaining_dims)
                if (remaining_dims > 0) {
                    linear_idx += adhoc::detail::multicombination(remaining_range, remaining_dims);
                }
                else {
                    linear_idx += 1; // Last dimension, just count
                }
            }
        }
        return this->m_data[linear_idx];
    }

    // Get the underlying dimension
    auto dim() const -> std::size_t { return m_dim; }

    auto size() const { return this->m_dim; }

    // Size of storage needed for symmetric tensor
    // Uses multicombination (stars and bars): C(dim+S-1, S)
    // This represents choosing S indices from dim values with repetition
    // allowed
    static auto storage_size(std::size_t dim) -> std::size_t
    {
        // For symmetric tensor of order S with dimension dim,
        // we need C(dim+S-1, S) storage locations
        return detail::multicombination(dim, S);
    }

    auto data() const -> T* { return m_data; }

  private:
    T* m_data;
    // All dimensions are equal for symmetric tensor
    std::size_t m_dim;
};

} // namespace adhoc

#endif // BMC_TOOLS_MDSPAN
