#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_UTIL_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_UTIL_HPP

#include "adhoc.hpp"
#include "dependency.hpp"
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

namespace detail {

template <std::size_t... Args> auto constexpr sum() -> std::size_t {
    return (Args + ...);
}

template <std::size_t Input1, std::size_t Input2, std::size_t... Inputs>
auto constexpr max() -> std::size_t {
    auto constexpr thismax = Input1 > Input2 ? Input1 : Input2;
    if constexpr (sizeof...(Inputs) == 0) {
        return thismax;
    } else {
        return max<thismax, Inputs...>();
    }
}

template <std::size_t Input1> auto constexpr max() -> std::size_t {
    return Input1;
}

// once we hit max, which represents infinity, we stay there
template <std::size_t Input1, std::size_t... Inputs>
auto constexpr sum_no_overflow2() -> std::size_t {
    if constexpr (sizeof...(Inputs) == 0) {
        return Input1;
    } else {
        if (std::numeric_limits<std::size_t>::max() - Input1 <
            sum_no_overflow<Inputs...>()) {
            return std::numeric_limits<std::size_t>::max();
        }

        if (std::numeric_limits<std::size_t>::max() -
                sum_no_overflow<Inputs...>() <
            Input1) {
            return std::numeric_limits<std::size_t>::max();
        }

        return Input1 + sum_no_overflow<Inputs...>();
    }
}

template <class IdNum, std::size_t OrderNum, std::size_t PowerNum,
          std::size_t N, std::size_t PowerDen>
auto constexpr derivative_non_null_simple(
    std::tuple<der2::p<PowerNum, der2::d<OrderNum, IdNum>>> /* ids1 */,
    der2::p<PowerDen, der2::d<1, double_t<N>>> /* id2 */) -> bool {
    auto constexpr summation = order<IdNum, double_t<N>>::call() * PowerNum;
    return summation >= PowerDen;
}

} // namespace detail

template <class... Ids, std::size_t... Orders1, std::size_t... Powers1,
          std::size_t... Ns, std::size_t... Orders2, std::size_t... Powers2>
// constexpr auto derivative_non_null(
auto derivative_non_null(
    std::tuple<der2::p<Powers1, der2::d<Orders1, Ids>>...> /* ids1 */,
    std::tuple<
        der2::p<Powers2, der2::d<Orders2, double_t<Ns>>>...> /* ids2 */) {

    auto constexpr numer_pow = detail::sum<Orders1 * Powers1...>();
    auto constexpr denom_pow = detail::sum<Orders2 * Powers2...>();

    if constexpr (denom_pow != numer_pow) {
        return false;
    } else {
        auto constexpr numer_max_d = detail::max<Orders1...>();
        auto constexpr denom_max_d = detail::max<Orders2...>();
        if constexpr (denom_max_d > numer_max_d) {
            return false;
        } else {
            if (denom_max_d == 1 && sizeof...(Ids) == 1) {
                return (
                    detail::derivative_non_null_simple(
                        std::tuple<
                            der2::p<Powers1, der2::d<Orders1, Ids>>...>{},
                        der2::p<Powers2, der2::d<Orders2, double_t<Ns>>>{}) &&
                    ...);
            } else {
                return true;
            }
        }
    }
}

// namespace detail {

// template <class... Ids, std::size_t... Orders, std::size_t... Powers>
// // constexpr auto expand_aux(
// auto expand_aux(
//     std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...> /* ids1 */) {
//     return std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...>{};
// }

// template <class Id1, class Id2, std::size_t Order1, std::size_t Power1,
//           class... Ids, std::size_t... Orders, std::size_t... Powers>
// // constexpr auto expand_aux(
// auto expand_aux(
//     std::tuple<der2::p<Power1, der2::d<Order1, add_t<Id1, Id2>>>,
//                der2::p<Powers, der2::d<Orders, Ids>>...> /* ids1 */) {
//     return std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...>{};
// }

// } // namespace detail

template <class... Ids, std::size_t... Orders, std::size_t... Powers>
// constexpr auto expand(
auto expand(std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...> /* ids1 */) {
    if constexpr (is_root_class_v<
                      std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...>>) {
        return std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...>{};
    } else {
        auto intermediate =
            expand_aux(std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...>{});
        return expand(intermediate);
    }
}

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_UTIL_HPP