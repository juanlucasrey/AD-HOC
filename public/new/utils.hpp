#ifndef ADHOC_UTILS_HPP
#define ADHOC_UTILS_HPP

#include <type_traits>

namespace adhoc2 {

namespace detail {
template <typename T, typename... Ts> struct get_index;

template <typename T, typename... Ts>
struct get_index<T, T, Ts...> : std::integral_constant<std::size_t, 0> {};

template <typename T, typename... Ts>
struct get_index<T, const T, Ts...> : std::integral_constant<std::size_t, 0> {};

template <typename T, typename... Ts>
struct get_index<const T, T, Ts...> : std::integral_constant<std::size_t, 0> {};

template <typename T, typename Tail, typename... Ts>
struct get_index<T, Tail, Ts...>
    : std::integral_constant<std::size_t, 1 + get_index<T, Ts...>::value> {};
} // namespace detail

template <typename T, typename... Ts> constexpr auto has_type() -> bool {
    return (std::is_same_v<T, Ts> || ...) ||
           (std::is_same_v<T, Ts const> || ...);
};

template <typename T, typename... Ts> constexpr auto idx_type() -> std::size_t {
    static_assert(has_type<T, Ts...>(), "variable not on this tape");
    return detail::get_index<T, Ts...>::value;
};

} // namespace adhoc2

#endif // ADHOC_UTILS_HPP
