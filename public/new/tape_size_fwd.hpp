#ifndef NEW_TAPE_SIZE_FWD_HPP
#define NEW_TAPE_SIZE_FWD_HPP

#include "adhoc.hpp"
#include "constants_type.hpp"
#include "dependency.hpp"
#include "utils.hpp"

#include <tuple>

namespace adhoc {

namespace detail {

template <typename... TypesAlive>
constexpr auto fwd_calc_order_t(const std::tuple<TypesAlive...> &);

template <> constexpr auto fwd_calc_order_t(const std::tuple<> &) {
    return std::tuple<>{};
}

#if __cplusplus >= 202002L

template <constants::detail::AsTemplateArg<double> D, typename... TypesAlive>
constexpr auto
fwd_calc_order_t(const std::tuple<constants::CD<D> const, TypesAlive...> &) {
    return fwd_calc_order_t(
        std::tuple<TypesAlive...>{}); // we don't add anything because it's on
                                      // the input tape
}

template <constants::detail::AsTemplateArg<double> D, typename... TypesAlive>
constexpr auto
fwd_calc_order_t(const std::tuple<constants::CD<D>, TypesAlive...> &) {
    return fwd_calc_order_t(
        std::tuple<TypesAlive...>{}); // we don't add anything because it's on
                                      // the input tape
}

#else

template <std::uint64_t D, typename... TypesAlive>
constexpr auto
fwd_calc_order_t(const std::tuple<constants::CD<D> const, TypesAlive...> &) {
    return fwd_calc_order_t(
        std::tuple<TypesAlive...>{}); // we don't add anything because it's on
                                      // the input tape
}

template <std::uint64_t D, typename... TypesAlive>
constexpr auto
fwd_calc_order_t(const std::tuple<constants::CD<D>, TypesAlive...> &) {
    return fwd_calc_order_t(
        std::tuple<TypesAlive...>{}); // we don't add anything because it's on
                                      // the input tape
}

#endif

template <std::size_t N, typename... TypesAlive>
constexpr auto
fwd_calc_order_t(const std::tuple<adouble_aux<N> const, TypesAlive...> &) {
    return fwd_calc_order_t(
        std::tuple<TypesAlive...>{}); // we don't add anything because it's on
                                      // the input tape
}

template <std::size_t N, typename... TypesAlive>
constexpr auto
fwd_calc_order_t(const std::tuple<adouble_aux<N>, TypesAlive...> &) {
    return fwd_calc_order_t(
        std::tuple<TypesAlive...>{}); // we don't add anything because it's on
                                      // the input tape
}

template <typename CurrentType, typename NextType, typename... TypesAlive>
constexpr auto
skip_type(const std::tuple<CurrentType, NextType, TypesAlive...> &) {
    return fwd_calc_order_t(
        std::tuple<NextType, TypesAlive...>{}); // current type will come up
                                                // again because it is included
                                                // on one of the other types
}

template <template <class...> class Xvariate, class... Input,
          typename... TypesAlive>
constexpr auto
fwd_calc_order_t(const std::tuple<Xvariate<Input...> const, TypesAlive...> &) {
    using this_type = Xvariate<Input...>;
    static_assert(!has_type2<this_type, TypesAlive...>());

    constexpr bool other_types_depend_on_this =
        (depends<TypesAlive, this_type>::call() || ...);

    if constexpr (other_types_depend_on_this) {
        return std::tuple<decltype(skip_type(
            std::tuple<this_type const, TypesAlive...>{}))>{};
    } else {
        return std::tuple<decltype(fwd_calc_order_t(
                              std::tuple<Input..., TypesAlive...>{})),
                          this_type>{};
    }
}

template <typename T> struct flatten_tuple { using type = std::tuple<T>; };

template <typename T> using flatten_tuple_t = typename flatten_tuple<T>::type;

template <typename... Ts> struct flatten_tuple<std::tuple<Ts...>> {
    using type =
        decltype(std::tuple_cat(std::declval<flatten_tuple_t<Ts>>()...));
};

} // namespace detail

template <class Output> constexpr auto fwd_calc_order_t(Output const &) {
    return detail::flatten_tuple_t<decltype(detail::fwd_calc_order_t(
        std::tuple<Output const>{}))>{};
}

} // namespace adhoc

#endif // NEW_TAPE_SIZE_FWD_HPP