#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_UTIL_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_UTIL_HPP

#include "adhoc.hpp"
#include "differential_operator_2.hpp"

namespace adhoc3 {

template <typename T> struct is_root_class : std::false_type {};

template <std::size_t... Ns, std::size_t... Orders, std::size_t... Powers>
struct is_root_class<
    std::tuple<der2::p<Powers, der2::d<Orders, double_t<Ns>>>...>>
    : std::true_type {};

template <class T> constexpr bool is_root_class_v = is_root_class<T>::value;

// tuple version
template <class... Ders>
constexpr auto select_root_derivatives(std::tuple<Ders...> /* ids */) {
    return std::tuple_cat(
        std::conditional_t<is_root_class_v<Ders>, std::tuple<Ders>,
                           std::tuple<>>{}...);
}

// tuple version
template <class... Ders>
constexpr auto select_non_root_derivatives(std::tuple<Ders...> /* ids */) {
    return std::tuple_cat(
        std::conditional_t<!is_root_class_v<Ders>, std::tuple<Ders>,
                           std::tuple<>>{}...);
}

// variadic version
template <class... Ders> constexpr auto select_root_derivatives2() {
    return std::tuple_cat(
        std::conditional_t<is_root_class_v<Ders>, std::tuple<Ders>,
                           std::tuple<>>{}...);
}

// variadic version
template <class... Ders> constexpr auto select_non_root_derivatives2() {
    return std::tuple_cat(
        std::conditional_t<!is_root_class_v<Ders>, std::tuple<Ders>,
                           std::tuple<>>{}...);
}

// template <class... Ids, std::size_t... Orders1, std::size_t... Powers1,
//           std::size_t... Ns, std::size_t... Orders2, std::size_t... Powers2>
// constexpr auto derivative_non_null(
//     std::tuple<der2::p<Powers1, der2::d<Orders1, Ids>>...> /* ids1 */,
//     std::tuple<
//         der2::p<Powers2, der2::d<Orders2, double_t<Ns>>>...> /* ids2 */) {
//     return true;
// }

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_UTIL_HPP