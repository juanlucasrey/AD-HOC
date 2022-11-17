#ifndef NEW_UTILS_HPP
#define NEW_UTILS_HPP

#include <type_traits>

namespace adhoc {

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

template <typename... Ts, typename T>
constexpr auto idx_type(T const & /* in */) -> std::size_t {
    static_assert((std::is_same_v<T, Ts> || ...));
    return detail::get_index<T, Ts...>::value;
};

template <typename T, typename... Ts> constexpr auto has_type2() -> bool {
    return (std::is_same_v<T, Ts> || ...) ||
           (std::is_same_v<T, Ts const> || ...);
};

template <typename T, typename... Ts>
constexpr auto idx_type2() -> std::size_t {
    static_assert(has_type2<T, Ts...>());
    return detail::get_index<T, Ts...>::value;
};

template <std::size_t N, std::size_t Head, std::size_t... Tail>
struct Get : Get<N - 1, Tail...> {};

template <std::size_t Head, std::size_t... Tail> struct Get<0, Head, Tail...> {
    static const std::size_t value = Head;
};

} // namespace adhoc

#endif // NEW_UTILS_HPP
