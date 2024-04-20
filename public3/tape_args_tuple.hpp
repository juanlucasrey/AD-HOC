/*
 * AD-HOC, Automatic Differentiation for High Order Calculations
 *
 * This file is part of the AD-HOC distribution
 * (https://github.com/juanlucasrey/AD-HOC).
 * Copyright (c) 2024 Juan Lucas Rey
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ADHOC3_TAPE_ARGS_TUPLE_HPP
#define ADHOC3_TAPE_ARGS_TUPLE_HPP

#include "adhoc.hpp"
#include "constants_type.hpp"
#include "ders.hpp"
#include "differential_operator.hpp"
#include "utils.hpp"

#include <tuple>

namespace adhoc3 {

namespace detail {

template <class... Ders, std::size_t... Orders>
constexpr auto d_order(std::tuple<der::p<der::d<Ders, 1>, Orders>...> /* in */)
    -> std::size_t {
    return (Orders + ...);
}

template <typename... OutputsAndDerivativesAlive> struct outputs_t {
    template <typename... Outputs> constexpr static auto call() noexcept;
};

template <> struct outputs_t<> {
    template <typename... Outputs> constexpr static auto call() noexcept {
        return std::tuple<Outputs...>{};
    }
};

template <typename... D, typename... OutputsAndDerivativesAlive>
struct outputs_t<std::tuple<D...>, OutputsAndDerivativesAlive...> {
    template <typename... Outputs> constexpr static auto call() noexcept {
        // we don't add anything because a derivative is not an output
        return outputs_t<OutputsAndDerivativesAlive...>::template call<
            Outputs...>();
    }
};

template <class D, typename... OutputsAndDerivativesAlive>
struct outputs_t<D, OutputsAndDerivativesAlive...> {
    template <typename... Outputs> constexpr static auto call() noexcept {
        return outputs_t<OutputsAndDerivativesAlive...>::template call<
            Outputs..., D>();
    }
};

template <typename... OutputsAndDerivativesAlive> struct derivatives_t {
    template <typename... Derivatives> constexpr static auto call() noexcept;
};

template <> struct derivatives_t<> {
    template <typename... Derivatives> constexpr static auto call() noexcept {
        return std::tuple<Derivatives...>{};
    }
};

template <typename... D, typename... OutputsAndDerivativesAlive>
struct derivatives_t<std::tuple<D...>, OutputsAndDerivativesAlive...> {
    template <typename... Derivatives> constexpr static auto call() noexcept {
        return derivatives_t<OutputsAndDerivativesAlive...>::template call<
            Derivatives..., std::tuple<D...>>();
    }
};

template <class D, typename... OutputsAndDerivativesAlive>
struct derivatives_t<D, OutputsAndDerivativesAlive...> {
    template <typename... Derivatives> constexpr static auto call() noexcept {
        // we don't add anything because an output is not a derivative
        return derivatives_t<OutputsAndDerivativesAlive...>::template call<
            Derivatives...>();
    }
};

template <typename... NodesAlive> struct inputs_t {
    template <typename... Inputs> constexpr static auto call() noexcept;
};

template <> struct inputs_t<> {
    template <typename... Inputs> constexpr static auto call() noexcept {
        return std::tuple<Inputs...>{};
    }
};

template <constants::ArgType D, typename... NodesAlive>
struct inputs_t<constants::CD<D>, NodesAlive...> {
    template <typename... Inputs> constexpr static auto call() noexcept {
        // we don't add anything because a constant is not an input
        return inputs_t<NodesAlive...>::template call<Inputs...>();
    }
};

template <std::size_t N, typename... NodesAlive>
struct inputs_t<double_t<N>, NodesAlive...> {
    template <typename... Inputs> constexpr static auto call() noexcept {
        using this_type = double_t<N>;

        constexpr bool already_included = has_type<this_type, Inputs...>();

        if constexpr (already_included) {
            // It's already been added, so we don't add the leaf
            return inputs_t<NodesAlive...>::template call<Inputs...>();
        } else {
            // we add the leaf
            return inputs_t<NodesAlive...>::template call<Inputs...,
                                                          this_type>();
        }
    }
};

template <template <class...> class Xvariate, class... Node,
          typename... NodesAlive>
struct inputs_t<Xvariate<Node...>, NodesAlive...> {
    template <typename... Inputs> constexpr static auto call() noexcept {
        return inputs_t<Node..., NodesAlive...>::template call<Inputs...>();
    }
};

template <typename... Derivatives> struct orders_t {
    template <std::size_t... Orders> constexpr static auto call() noexcept;
};

template <> struct orders_t<> {
    template <std::size_t... Orders> constexpr static auto call() noexcept {
        return std::index_sequence<Orders...>{};
    }
};

template <typename D, typename... Derivatives>
struct orders_t<D, Derivatives...> {
    template <std::size_t... Orders> constexpr static auto call() noexcept {
        auto constexpr this_order = d_order(D{});

        constexpr bool already_included = ((this_order == Orders) || ...);

        if constexpr (already_included) {
            // It's already been added, so we don't add the new order
            return orders_t<Derivatives...>::template call<Orders...>();
        } else {
            // we add the order
            return orders_t<Derivatives...>::template call<Orders...,
                                                           this_order>();
        }
    }
};

template <class Output, std::size_t... Orders> constexpr auto expand() {
    return std::tuple<der::d<Output, Orders>...>{};
}

} // namespace detail

template <class... OutputsAndDerivatives>
constexpr auto Outputs(OutputsAndDerivatives... /* o */) {
    return detail::outputs_t<OutputsAndDerivatives...>::template call<>();
}

template <class... OutputsAndDerivatives>
constexpr auto Derivatives(OutputsAndDerivatives... /* o */) {
    return detail::derivatives_t<OutputsAndDerivatives...>::template call<>();
}

template <class... Outputs>
constexpr auto Inputs(std::tuple<Outputs...> /* o */) {
    return detail::inputs_t<Outputs...>::template call<>();
}

template <class... Derivatives>
constexpr auto Orders(std::tuple<Derivatives...> /* o */) {
    return detail::orders_t<Derivatives...>::template call<>();
}

template <class... Outputs, std::size_t... Orders>
constexpr auto ExpandOutputs(std::tuple<Outputs...> /* o */,
                             std::index_sequence<Orders...> /* orders */) {
    return std::tuple_cat((detail::expand<Outputs, Orders...>())...);
}

} // namespace adhoc3

#endif // ADHOC3_TAPE_ARGS_TUPLE_HPP
