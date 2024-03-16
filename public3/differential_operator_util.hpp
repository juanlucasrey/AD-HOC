// #ifndef ADHOC3_DIFFERENTIAL_OPERATOR_UTIL_HPP
// #define ADHOC3_DIFFERENTIAL_OPERATOR_UTIL_HPP

// #include "adhoc.hpp"
// #include "dependency.hpp"
// #include "differential_operator/differential_operator.hpp"

// namespace adhoc3 {

// template <typename T> struct is_root_class : std::false_type {};

// template <std::size_t... Ns, std::size_t... Orders, std::size_t... Powers>
// struct is_root_class<
//     std::tuple<der2::p<Powers, der2::d<Orders, double_t<Ns>>>...>>
//     : std::true_type {};

// template <class T> constexpr bool is_root_class_v = is_root_class<T>::value;

// // tuple version
// template <class... Ders>
// constexpr auto select_root_derivatives(std::tuple<Ders...> /* ids */) {
//     return std::tuple_cat(
//         std::conditional_t<is_root_class_v<Ders>, std::tuple<Ders>,
//                            std::tuple<>>{}...);
// }

// // tuple version
// template <class... Ders>
// constexpr auto select_non_root_derivatives(std::tuple<Ders...> /* ids */) {
//     return std::tuple_cat(
//         std::conditional_t<!is_root_class_v<Ders>, std::tuple<Ders>,
//                            std::tuple<>>{}...);
// }

// // variadic version
// template <class... Ders> constexpr auto select_root_derivatives2() {
//     return std::tuple_cat(
//         std::conditional_t<is_root_class_v<Ders>, std::tuple<Ders>,
//                            std::tuple<>>{}...);
// }

// // variadic version
// template <class... Ders> constexpr auto select_non_root_derivatives2() {
//     return std::tuple_cat(
//         std::conditional_t<!is_root_class_v<Ders>, std::tuple<Ders>,
//                            std::tuple<>>{}...);
// }

// namespace detail {

// template <std::size_t... Args> auto constexpr sum() -> std::size_t {
//     return (Args + ...);
// }

// template <std::size_t Input1, std::size_t Input2, std::size_t... Inputs>
// auto constexpr max() -> std::size_t {
//     auto constexpr thismax = Input1 > Input2 ? Input1 : Input2;
//     if constexpr (sizeof...(Inputs) == 0) {
//         return thismax;
//     } else {
//         return max<thismax, Inputs...>();
//     }
// }

// template <std::size_t Input1> auto constexpr max() -> std::size_t {
//     return Input1;
// }

// // once we hit max, which represents infinity, we stay there
// template <std::size_t Input1, std::size_t... Inputs>
// auto constexpr sum_no_overflow2() -> std::size_t {
//     if constexpr (sizeof...(Inputs) == 0) {
//         return Input1;
//     } else {
//         if (std::numeric_limits<std::size_t>::max() - Input1 <
//             sum_no_overflow<Inputs...>()) {
//             return std::numeric_limits<std::size_t>::max();
//         }

//         if (std::numeric_limits<std::size_t>::max() -
//                 sum_no_overflow<Inputs...>() <
//             Input1) {
//             return std::numeric_limits<std::size_t>::max();
//         }

//         return Input1 + sum_no_overflow<Inputs...>();
//     }
// }

// template <class IdNum, std::size_t OrderNum, std::size_t PowerNum,
//           std::size_t N, std::size_t PowerDen>
// auto constexpr derivative_non_null_simple(
//     std::tuple<der2::p<PowerNum, der2::d<OrderNum, IdNum>>> /* ids1 */,
//     der2::p<PowerDen, der2::d<1, double_t<N>>> /* id2 */) -> bool {
//     auto constexpr summation = order<IdNum, double_t<N>>::call() * PowerNum;
//     return summation >= PowerDen;
// }

// } // namespace detail

// template <class... Ids, std::size_t... Orders1, std::size_t... Powers1,
//           std::size_t... Ns, std::size_t... Orders2, std::size_t... Powers2>
// // constexpr auto derivative_non_null(
// auto derivative_non_null(
//     std::tuple<der2::p<Powers1, der2::d<Orders1, Ids>>...> /* ids1 */,
//     std::tuple<
//         der2::p<Powers2, der2::d<Orders2, double_t<Ns>>>...> /* ids2 */) {

//     auto constexpr numer_pow = detail::sum<Orders1 * Powers1...>();
//     auto constexpr denom_pow = detail::sum<Orders2 * Powers2...>();

//     if constexpr (denom_pow != numer_pow) {
//         return false;
//     } else {
//         auto constexpr numer_max_d = detail::max<Orders1...>();
//         auto constexpr denom_max_d = detail::max<Orders2...>();
//         if constexpr (denom_max_d > numer_max_d) {
//             return false;
//         } else {
//             if (denom_max_d == 1 && sizeof...(Ids) == 1) {
//                 return (
//                     detail::derivative_non_null_simple(
//                         std::tuple<
//                             der2::p<Powers1, der2::d<Orders1, Ids>>...>{},
//                         der2::p<Powers2, der2::d<Orders2, double_t<Ns>>>{})
//                         &&
//                     ...);
//             } else {
//                 return true;
//             }
//         }
//     }
// }

// // namespace detail {

// // template <class... Ids, std::size_t... Orders, std::size_t... Powers>
// // // constexpr auto expand_aux(
// // auto expand_aux(
// //     std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...> /* ids1 */) {
// //     return std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...>{};
// // }

// // template <class Id1, class Id2, std::size_t Order1, std::size_t Power1,
// //           class... Ids, std::size_t... Orders, std::size_t... Powers>
// // // constexpr auto expand_aux(
// // auto expand_aux(
// //     std::tuple<der2::p<Power1, der2::d<Order1, add_t<Id1, Id2>>>,
// //                der2::p<Powers, der2::d<Orders, Ids>>...> /* ids1 */) {
// //     return std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...>{};
// // }

// // } // namespace detail

// template <class... Ids, std::size_t... Orders, std::size_t... Powers>
// // constexpr auto expand(
// auto expand(std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...> /* ids1 */)
// {
//     if constexpr (is_root_class_v<
//                       std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...>>)
//                       {
//         return std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...>{};
//     } else {
//         auto intermediate =
//             expand_aux(std::tuple<der2::p<Powers, der2::d<Orders,
//             Ids>>...>{});
//         return expand(intermediate);
//     }
// }

// namespace detail {

// template <class In, class Id, std::size_t Power, std::size_t Order>
// constexpr auto
// create_single_tuple2(In /* in */, der2::p<Power, der2::d<Order, Id>> /* o */)
// {
//     if constexpr (std::is_same_v<In, Id>) {
//         return std::tuple<der2::p<Power, der2::d<Order, Id>>>{};
//     } else {
//         return std::tuple<>{};
//     }
// };

// template <std::size_t I, class Id, std::size_t Power, std::size_t Order>
// constexpr auto order_by_order4(der2::p<Power, der2::d<Order, Id>> /* o */) {
//     if constexpr (I == Order) {
//         return std::tuple<der2::p<Power, der2::d<Order, Id>>>{};
//     } else {
//         return std::tuple<>{};
//     }
// }

// template <std::size_t I, class Id, std::size_t... Powers, std::size_t...
// Orders> constexpr auto order_by_order3(std::tuple<der2::p<Powers,
// der2::d<Orders, Id>>...> o) {

//     return std::apply(
//         [](auto &&...args) {
//             return std::tuple_cat(order_by_order4<I>(args)...);
//         },
//         o);
// }

// template <class Id, std::size_t... Powers, std::size_t... Orders,
//           std::size_t... I>
// constexpr auto
// order_by_order2(std::tuple<der2::p<Powers, der2::d<Orders, Id>>...> o,
//                 std::index_sequence<I...> /* i */) {
//     return std::tuple_cat(order_by_order3<I + 1>(o)...);
// }

// template <class In> constexpr auto order_by_order(In /* o */);

// template <class Id, std::size_t... Powers, std::size_t... Orders>
// constexpr auto
// order_by_order(std::tuple<der2::p<Powers, der2::d<Orders, Id>>...> o) {
//     auto constexpr numer_max_d = detail::max<Orders...>();
//     auto seq = std::make_integer_sequence<std::size_t, numer_max_d>{};
//     return order_by_order2(o, seq);
// };

// template <> constexpr auto order_by_order(std::tuple<> /* o */) {
//     return std::tuple<>{};
// };

// template <class In, class... Ids, std::size_t... Powers, std::size_t...
// Orders> constexpr auto create_single_tuple(In in,
//                     std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...> o) {
//     auto tuple_to_order = std::apply(
//         [in](auto &&...args) {
//             return std::tuple_cat(detail::create_single_tuple2(in, args)...);
//         },
//         o);

//     return order_by_order(tuple_to_order);
// };

// } // namespace detail

// template <class... Ids, std::size_t... Orders, std::size_t... Powers,
//           class... CalculationNodes>
// constexpr auto order_differential_operator(
//     std::tuple<der2::p<Powers, der2::d<Orders, Ids>>...> diff_operator,
//     std::tuple<CalculationNodes...> nodes) {
//     return std::apply(
//         [diff_operator](auto &&...args) {
//             return std::tuple_cat(
//                 detail::create_single_tuple(args, diff_operator)...);
//         },
//         nodes);
// }

// template <class... Ids1, std::size_t... Orders1, std::size_t... Powers1,
//           class... Ids2, std::size_t... Orders2, std::size_t... Powers2,
//           class... CalculationNodes>
// constexpr auto
// less_than(std::tuple<der2::p<Powers1, der2::d<Orders1, Ids1>>...> in1,
//           std::tuple<der2::p<Powers2, der2::d<Orders2, Ids2>>...> in2,
//           std::tuple<CalculationNodes...> nodes);

// namespace detail {

// template <class Id1, class... Ids1, std::size_t Order1, std::size_t...
// Orders1,
//           std::size_t Power1, std::size_t... Powers1, class Id2, class...
//           Ids2, std::size_t Order2, std::size_t... Orders2, std::size_t
//           Power2, std::size_t... Powers2, class... CalculationNodes>
// constexpr auto less_than_check_first(
//     std::tuple<der2::p<Power1, der2::d<Order1, Id1>>,
//                der2::p<Powers1, der2::d<Orders1, Ids1>>...>
//     /* in1 */,
//     std::tuple<der2::p<Power2, der2::d<Order2, Id2>>,
//                der2::p<Powers2, der2::d<Orders2, Ids2>>...> /* in2 */,
//     std::tuple<CalculationNodes...> nodes);

// template <class... Ids1, std::size_t... Orders1, std::size_t... Powers1,
//           class... Ids2, std::size_t... Orders2, std::size_t... Powers2,
//           class... CalculationNodes>
// constexpr auto less_than_check_empty(
//     std::tuple<der2::p<Powers1, der2::d<Orders1, Ids1>>...> in1,
//     std::tuple<der2::p<Powers2, der2::d<Orders2, Ids2>>...> in2,
//     std::tuple<CalculationNodes...> nodes) {

//     if constexpr (std::is_same<decltype(in1), decltype(in2)>::value) {
//         return false;
//     } else if constexpr (std::is_same<decltype(in1), std::tuple<>>::value) {
//         return true;
//     } else if constexpr (std::is_same<decltype(in2), std::tuple<>>::value) {
//         return false;
//     } else {
//         // we can now assume that both ids are different and non empty
//         return less_than_check_first(in1, in2, nodes);
//     }
// }

// template <class Id1, class... Ids1, std::size_t Order1, std::size_t...
// Orders1,
//           std::size_t Power1, std::size_t... Powers1, class Id2, class...
//           Ids2, std::size_t Order2, std::size_t... Orders2, std::size_t
//           Power2, std::size_t... Powers2, class... CalculationNodes>
// constexpr auto less_than_check_first(
//     std::tuple<der2::p<Power1, der2::d<Order1, Id1>>,
//                der2::p<Powers1, der2::d<Orders1, Ids1>>...>
//     /* in1 */,
//     std::tuple<der2::p<Power2, der2::d<Order2, Id2>>,
//                der2::p<Powers2, der2::d<Orders2, Ids2>>...> /* in2 */,
//     std::tuple<CalculationNodes...> nodes) {
//     constexpr auto idx1 = get_idx_first2<Id1>(nodes);
//     constexpr auto idx2 = get_idx_first2<Id2>(nodes);
//     if constexpr (idx1 < idx2) {
//         return true;
//     } else if constexpr (idx1 > idx2) {
//         return false;
//     } else if constexpr (Order1 < Order2) {
//         return true;
//     } else if constexpr (Order1 > Order2) {
//         return false;
//     } else if constexpr (Power1 < Power2) {
//         return true;
//     } else if constexpr (Power1 > Power2) {
//         return false;
//     } else {
//         static_assert(
//             std::is_same<der2::p<Power1, der2::d<Order1, Id1>>,
//                          der2::p<Power2, der2::d<Order2, Id2>>>::value);
//         return less_than_check_empty(
//             std::tuple<der2::p<Powers1, der2::d<Orders1, Ids1>>...>{},
//             std::tuple<der2::p<Powers2, der2::d<Orders2, Ids2>>...>{},
//             nodes);
//     }
// }
// } // namespace detail

// template <class... Ids1, std::size_t... Orders1, std::size_t... Powers1,
//           class... Ids2, std::size_t... Orders2, std::size_t... Powers2,
//           class... CalculationNodes>
// constexpr auto
// less_than(std::tuple<der2::p<Powers1, der2::d<Orders1, Ids1>>...> in1,
//           std::tuple<der2::p<Powers2, der2::d<Orders2, Ids2>>...> in2,
//           std::tuple<CalculationNodes...> nodes) {
//     auto constexpr ordered_in1 = order_differential_operator(in1, nodes);
//     auto constexpr ordered_in2 = order_differential_operator(in2, nodes);
//     return detail::less_than_check_empty(ordered_in1, ordered_in2, nodes);
// }

// namespace detail {

// template <class... CalculationNodes>
// constexpr auto is_ordered_aux(std::tuple<> /* in */,
//                               std::tuple<CalculationNodes...> /* nodes */) {

//     return true;
// }

// template <class Id, class... CalculationNodes>
// constexpr auto is_ordered_aux(std::tuple<Id> /* in */,
//                               std::tuple<CalculationNodes...> /* nodes */) {

//     return true;
// }

// template <class Id1, class Id2, class... Ids, class... CalculationNodes>
// constexpr auto is_ordered_aux(std::tuple<Id1, Id2, Ids...> /* in */,
//                               std::tuple<CalculationNodes...> nodes) {

//     if constexpr (less_than(Id2{}, Id1{}, nodes)) {
//         return is_ordered_aux(std::tuple<Id2, Ids...>{}, nodes);
//     } else {
//         return false;
//     }
// }

// } // namespace detail

// template <class... Ids, class... CalculationNodes>
// constexpr auto is_ordered(std::tuple<Ids...> in,
//                           std::tuple<CalculationNodes...> nodes) {

//     return detail::is_ordered_aux(in, nodes);
// }

// namespace detail {

// template <class Id1, class... Ids1, class Id2, class... Ids2, class... Out,
//           class... CalculationNodes>
// constexpr auto merge_ordered_compare(std::tuple<Id1, Ids1...> in1,
//                                      std::tuple<Id2, Ids2...> in2,
//                                      std::tuple<Out...> /* out */,
//                                      std::tuple<CalculationNodes...> nodes);

// template <class... Ids1, class... Ids2, class... Out, class...
// CalculationNodes> constexpr auto merge_ordered_empty(std::tuple<Ids1...> in1,
//                                    std::tuple<Ids2...> in2,
//                                    std::tuple<Out...> out,
//                                    std::tuple<CalculationNodes...> nodes) {

//     if constexpr (sizeof...(Ids2) == 0) {
//         return std::tuple_cat(out, in1);
//     } else if constexpr (sizeof...(Ids1) == 0) {
//         return std::tuple_cat(out, in2);
//     } else {
//         return merge_ordered_compare(in1, in2, out, nodes);
//     }
// }

// template <class Id1, class... Ids1, class Id2, class... Ids2, class... Out,
//           class... CalculationNodes>
// constexpr auto merge_ordered_compare(std::tuple<Id1, Ids1...> in1,
//                                      std::tuple<Id2, Ids2...> in2,
//                                      std::tuple<Out...> /* out */,
//                                      std::tuple<CalculationNodes...> nodes) {
//     if constexpr (std::is_same<Id1, Id2>::value) {
//         return merge_ordered_empty(std::tuple<Ids1...>{},
//         std::tuple<Ids2...>{},
//                                    std::tuple<Out..., Id1>{}, nodes);

//     } else if constexpr (less_than(Id1{}, Id2{}, nodes)) {
//         return merge_ordered_empty(in1, std::tuple<Ids2...>{},
//                                    std::tuple<Out..., Id2>{}, nodes);

//     } else {
//         return merge_ordered_empty(std::tuple<Ids1...>{}, in2,
//                                    std::tuple<Out..., Id1>{}, nodes);
//     }
// }

// } // namespace detail

// template <class... Ids1, class... Ids2, class... CalculationNodes>
// constexpr auto merge_ordered(std::tuple<Ids1...> in1, std::tuple<Ids2...>
// in2,
//                              std::tuple<CalculationNodes...> nodes) {

//     static_assert(is_ordered(in1, nodes));
//     static_assert(is_ordered(in2, nodes));
//     return detail::merge_ordered_empty(in1, in2, std::tuple<>{}, nodes);
// }

// } // namespace adhoc3

// #endif // ADHOC3_DIFFERENTIAL_OPERATOR_UTIL_HPP
