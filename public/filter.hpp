#ifndef FILTER_HPP
#define FILTER_HPP

#include <tuple>
#include <type_traits>

namespace adhoc2 {

// template <typename...> struct filter;

// template <> struct filter<> { using type = std::tuple<>; };

// template <typename Head, typename ...Tail>
// struct filter<Head, Tail...>
// {
//     using type = typename std::conditional<Predicate<Head>::value,
//                                typename Cons<Head, typename
//                                filter<Tail...>::type>::type, typename
//                                filter<Tail...>::type
//                           >::type;
// };

// template <template <class> class, template <class...> class, class...>
// struct filter;
// template <template <class> class Pred, template <class...> class Variadic>
// struct filter<Pred, Variadic> {
//     using type = Variadic<>;
// };

// template <template <class> class Pred, template <class...> class Variadic,
//           class T, class... Ts>
// struct filter<Pred, Variadic, T, Ts...> {
//     template <class, class> struct Cons;
//     template <class Head, class... Tail> struct Cons<Head, Variadic<Tail...>>
//     {
//         using type = Variadic<Head, Tail...>;
//     };

//     using type = typename std::conditional<
//         Pred<T>::value,
//         typename Cons<T, typename filter<Pred, Variadic, Ts...>::type>::type,
//         typename filter<Pred, Variadic, Ts...>::type>::type;
// };

} // namespace adhoc2

#endif // FILTER_HPP
