#ifndef UTILS_HPP
#define UTILS_HPP

#include <type_traits>

namespace adhoc2 {

namespace detail {
template <typename T, typename... Ts> struct get_index;

template <typename T, typename... Ts>
struct get_index<T, T, Ts...> : std::integral_constant<std::size_t, 0> {};

template <typename T, typename Tail, typename... Ts>
struct get_index<T, Tail, Ts...>
    : std::integral_constant<std::size_t, 1 + get_index<T, Ts...>::value> {};
} // namespace detail

template <typename... Ts, typename T>
constexpr auto has_type(T const & /* in */) -> bool {
    return (std::is_same_v<T, Ts> || ...);
};

template <typename... Ts, typename T>
constexpr auto idx_type(T const & /* in */) -> std::size_t {
    static_assert((std::is_same_v<T, Ts> || ...));
    return detail::get_index<T, Ts...>::value;
};

} // namespace adhoc2

#endif // UTILS_HPP
