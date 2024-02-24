#ifndef ADHOC3_DIFFERENTIAL_OPERATOR_UTIL_HPP
#define ADHOC3_DIFFERENTIAL_OPERATOR_UTIL_HPP

#include "adhoc.hpp"
#include "differential_operator_2.hpp"

namespace adhoc3 {

template <std::size_t... Ns, std::size_t... Orders, std::size_t... Powers>
constexpr auto is_root_derivative(
    std::tuple<der2::p<Powers, der2::d<Orders, double_t<Ns>>>...> /* id */) {
    return true;
}

template <class... Ids, std::size_t... Orders, std::size_t... Powers>
constexpr auto is_root_derivative(
    std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...> /* id */) {
    return false;
}

} // namespace adhoc3

#endif // ADHOC3_DIFFERENTIAL_OPERATOR_UTIL_HPP