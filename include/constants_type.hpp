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

#ifndef CONSTANTS_TYPE_HPP
#define CONSTANTS_TYPE_HPP

#include "base.hpp"

namespace adhoc4::constants {

template <class Derived> struct ConstBase {
    operator double() const { return Derived::v(); }
};

namespace detail {
template <typename T> struct AsTemplateArg {
    T buffer{};
    constexpr AsTemplateArg(T t) : buffer(t) {}
    constexpr operator T() const { return this->buffer; }
};
} // namespace detail

template <typename T, detail::AsTemplateArg<T> F>
struct C : public Base<C<T, F>> {
    constexpr operator T() const { return F; }
    constexpr auto static v() -> T { return F; }
};

template <detail::AsTemplateArg<float> F> using CF = C<float, F>;
template <detail::AsTemplateArg<double> F> using CD = C<double, F>;
template <detail::AsTemplateArg<long double> F> using CLD = C<long double, F>;

using ArgType = constants::detail::AsTemplateArg<double>;

} // namespace adhoc4::constants

#endif // CONSTANTS_TYPE_HPP
