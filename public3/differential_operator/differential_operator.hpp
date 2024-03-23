#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_DIFFERENTIAL_OPERATOR_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_DIFFERENTIAL_OPERATOR_HPP

#include "../dependency.hpp"
#include "../utils.hpp"

#include <tuple>

namespace adhoc3 {

namespace der2 {

template <std::size_t Order, class Id> struct d {};
template <std::size_t Power, class Id> struct p {};

} // namespace der2

template <std::size_t Order = 1UL, class Id> constexpr auto d(Id /* id */) {
    return std::tuple<der2::p<1UL, der2::d<Order, Id>>>{};
}

template <std::size_t Power, class... Ids, std::size_t... Orders,
          std::size_t... Powers>
constexpr auto
pow(std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...> /* id */) {
    if constexpr (Power == 0) {
        return std::tuple<>{};
    } else {
        return std::tuple<der2::p<Powers * Power, der2::d<Orders, Ids>>...>{};
    }
}

namespace detail {

template <class Id, std::size_t Order, std::size_t Power1, std::size_t Power2>
constexpr auto add_duplicate(der2::p<Power1, der2::d<Order, Id>> /* in1 */,
                             der2::p<Power2, der2::d<Order, Id>> /* in2 */) {
    return der2::p<Power1 + Power2, der2::d<Order, Id>>{};
}

template <class Id1, std::size_t Order1, std::size_t Power1, class... Ids2,
          std::size_t... Orders2, std::size_t... Powers2>
constexpr auto
add_duplicates(der2::p<Power1, der2::d<Order1, Id1>> in,
               std::tuple<der2::p<Powers2, der2::d<Orders2, Ids2>>...> id2) {
    if constexpr (has_type<der2::d<Order1, Id1>, der2::d<Orders2, Ids2>...>()) {
        constexpr auto idx =
            idx_type2<der2::d<Order1, Id1>, der2::d<Orders2, Ids2>...>();
        constexpr auto v = std::get<idx>(id2);
        constexpr auto added = add_duplicate(in, v);
        return std::make_tuple(added);
    } else {
        return std::make_tuple(in);
    }
};

template <class Id1, std::size_t Order1, std::size_t Power1, class... Ids2,
          std::size_t... Orders2, std::size_t... Powers2>
constexpr auto remove_duplicates(
    der2::p<Power1, der2::d<Order1, Id1>> in,
    std::tuple<der2::p<Powers2, der2::d<Orders2, Ids2>>...> /* id2 */) {
    if constexpr (has_type<der2::d<Order1, Id1>, der2::d<Orders2, Ids2>...>()) {
        return std::tuple{};
    } else {
        return std::make_tuple(in);
    }
};

} // namespace detail

template <class... Ids1, std::size_t... Orders1, std::size_t... Powers1,
          class... Ids2, std::size_t... Orders2, std::size_t... Powers2>
constexpr auto
operator*(std::tuple<der2::p<Powers1, der2::d<Orders1, Ids1>>...> id1,
          std::tuple<der2::p<Powers2, der2::d<Orders2, Ids2>>...> id2) {
    auto augmented_in1 = std::apply(
        [id2](auto... id1s) {
            return std::tuple_cat(detail::add_duplicates(id1s, id2)...);
        },
        id1);

    auto reduced_in2 = std::apply(
        [id1](auto... id2s) {
            return std::tuple_cat(detail::remove_duplicates(id2s, id1)...);
        },
        id2);

    return std::tuple_cat(augmented_in1, reduced_in2);
}

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_DIFFERENTIAL_OPERATOR_HPP
