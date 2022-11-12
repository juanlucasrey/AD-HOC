#ifndef EVALUATE_BWD_HPP
#define EVALUATE_BWD_HPP

// #include <adhoc3.hpp>
// #include <dependency3.hpp>
// #include <tape_size.hpp>

#ifdef CODELOGGER
#include <iostream>
#endif

namespace adhoc3 {
namespace detail {

// template <std::size_t I, std::size_t M, std::size_t N, typename...
// LeavesAlive,
//           std::size_t... IdxLeavesAlive>
// static inline void
// evaluate_final(std::array<double, M> & /* temp */,
//                std::array<double, N> const & /* tape */,
//                args<LeavesAlive...> const &, args<> const &,
//                std::index_sequence<IdxLeavesAlive...> const &) {}

// template <std::size_t I, std::size_t M, std::size_t N, typename...
// LeavesAlive,
//           typename CurrentLeaf, typename... Leaves,
//           std::size_t... IdxLeavesAlive>
// static inline void
// evaluate_final(std::array<double, M> &temp, std::array<double, N> const
// &tape,
//                args<LeavesAlive...> const &,
//                args<CurrentLeaf, Leaves...> const &,
//                std::index_sequence<IdxLeavesAlive...> const &) {

//     // what if a leaf is never alive? treat separately
//     constexpr bool is_leaf_active =
//         has_type2<CurrentLeaf const, LeavesAlive...>();

//     if constexpr (is_leaf_active) {
//         if constexpr (I < M) {
//             constexpr auto position =
//                 idx_type2<CurrentLeaf const, LeavesAlive...>();
//             constexpr auto position_on_tape =
//                 Get<position, IdxLeavesAlive...>::value;
//             temp[I] = tape[position_on_tape];
//             evaluate_final<I + 1>(temp, tape, args<LeavesAlive...>{},
//                                   args<Leaves...>{},
//                                   std::index_sequence<IdxLeavesAlive...>{});
//         }
//     } else {
//         evaluate_final<I>(temp, tape, args<LeavesAlive...>{},
//         args<Leaves...>{},
//                           std::index_sequence<IdxLeavesAlive...>{});
//     }
// }

// template <std::size_t N, typename... LeavesAlive, typename... Leaves,
//           std::size_t... IdxLeavesAlive, std::size_t... IdxAvailable>
// static inline void evaluate(args<LeavesAlive...> const &,
//                             args<Leaves...> const &,
//                             std::index_sequence<> const &,
//                             std::index_sequence<IdxLeavesAlive...> const &,
//                             std::index_sequence<IdxAvailable...> const &,
//                             std::array<double, N> &tape) {

//     // final stage, we order tape for final output
//     constexpr std::size_t final_result_size = sizeof...(IdxLeavesAlive);
//     std::array<double, final_result_size> temp;
//     evaluate_final<0>(temp, tape, args<LeavesAlive...>{}, args<Leaves...>{},
//                       std::index_sequence<IdxLeavesAlive...>{});

//     for (std::size_t index = 0; index < final_result_size; ++index) {
//         tape[index] = temp[index];
//     }
// }

// template <std::size_t N, class FirstType, class SecondType,
//           typename... TypesAlive, typename... LeavesAlive, typename...
//           Leaves, std::size_t IdxTypesAliveFirst, std::size_t
//           IdxTypesAliveSecond, std::size_t... IdxTypesAlive, std::size_t...
//           IdxLeavesAlive, std::size_t... IdxAvailable>
// static inline void
// evaluate_skip(args<LeavesAlive...> const &, args<Leaves...> const &,
//               std::index_sequence<IdxTypesAliveFirst, IdxTypesAliveSecond,
//                                   IdxTypesAlive...> const &,
//               std::index_sequence<IdxLeavesAlive...> const &,
//               std::index_sequence<IdxAvailable...> const &,
//               std::array<double, N> &tape, FirstType const &first,
//               SecondType const &second, TypesAlive const &...next) {

//     evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
//              std::index_sequence<
//                  IdxTypesAliveSecond, IdxTypesAlive...,
//                  IdxTypesAliveFirst>{} /* we place first at the end */,
//              std::index_sequence<IdxLeavesAlive...>{},
//              std::index_sequence<IdxAvailable...>{}, tape, second, next...,
//              first /* we place first at the end */);
// }

// template <class Input, std::size_t N, class this_type, typename...
// TypesAlive,
//           typename... LeavesAlive, typename... Leaves,
//           std::size_t IdxTypesAliveCurrent, std::size_t... IdxTypesAlive,
//           std::size_t... IdxLeavesAlive, std::size_t... IdxAvailable>
// static inline void evaluate_univariate_noskip(
//     args<LeavesAlive...> const &, args<Leaves...> const &,
//     std::index_sequence<IdxTypesAliveCurrent, IdxTypesAlive...> const &,
//     std::index_sequence<IdxLeavesAlive...> const &,
//     std::index_sequence<IdxAvailable...> const &, std::array<double, N>
//     &tape, this_type const &in, TypesAlive const &...next) {

//     constexpr bool is_input_leaf = has_type2<Input, Leaves...>();
//     if constexpr (is_input_leaf) {
//         constexpr bool is_input_new_leaf = !has_type2<Input,
//         LeavesAlive...>();

//         if constexpr (is_input_new_leaf) {
// #ifdef CODELOGGER
//             std::cout << "tape[" << IdxTypesAliveCurrent << "] *= ";
// #endif
//             tape[IdxTypesAliveCurrent] *= in.d();
// #ifdef CODELOGGER
//             std::cout << ";" << std::endl;
// #endif
//             evaluate(
//                 args<Input, LeavesAlive...>{}, args<Leaves...>{},
//                 std::index_sequence<IdxTypesAlive...>{},
//                 std::index_sequence<IdxTypesAliveCurrent,
//                 IdxLeavesAlive...>{}, std::index_sequence<IdxAvailable...>{},
//                 tape, next...);
//         } else {
//             constexpr auto position = idx_type2<Input, LeavesAlive...>();
//             constexpr auto position_on_tape =
//                 Get<position, IdxLeavesAlive...>::value;
// #ifdef CODELOGGER
//             std::cout << "tape[" << position_on_tape << "] += tape["
//                       << IdxTypesAliveCurrent << "] * ";
// #endif
//             tape[position_on_tape] += tape[IdxTypesAliveCurrent] * in.d();
// #ifdef CODELOGGER
//             std::cout << ";" << std::endl;
// #endif
//             evaluate(
//                 args<LeavesAlive...>{}, args<Leaves...>{},
//                 std::index_sequence<IdxTypesAlive...>{},
//                 std::index_sequence<IdxLeavesAlive...>{},
//                 std::index_sequence<IdxTypesAliveCurrent, IdxAvailable...>{},
//                 tape, next...);
//         }
//     } else {
//         constexpr bool is_input_new = !has_type2<Input, TypesAlive...>();

//         if constexpr (is_input_new) {
// #ifdef CODELOGGER
//             std::cout << "tape[" << IdxTypesAliveCurrent << "] *= ";
// #endif
//             tape[IdxTypesAliveCurrent] *= in.d();
// #ifdef CODELOGGER
//             std::cout << ";" << std::endl;
// #endif
//             evaluate(
//                 args<LeavesAlive...>{}, args<Leaves...>{},
//                 std::index_sequence<IdxTypesAliveCurrent,
//                 IdxTypesAlive...>{},
//                 std::index_sequence<IdxLeavesAlive...>{},
//                 std::index_sequence<IdxAvailable...>{}, tape, in.input(),
//                 next...);
//         } else {
//             constexpr auto position = idx_type2<Input, TypesAlive...>();
//             constexpr auto position_on_tape =
//                 Get<position, IdxTypesAlive...>::value;
// #ifdef CODELOGGER
//             std::cout << "tape[" << position_on_tape << "] += tape["
//                       << IdxTypesAliveCurrent << "] * ";
// #endif
//             tape[position_on_tape] += tape[IdxTypesAliveCurrent] * in.d();
// #ifdef CODELOGGER
//             std::cout << ";" << std::endl;
// #endif
//             evaluate(
//                 args<LeavesAlive...>{}, args<Leaves...>{},
//                 std::index_sequence<IdxTypesAlive...>{},
//                 std::index_sequence<IdxLeavesAlive...>{},
//                 std::index_sequence<IdxTypesAliveCurrent, IdxAvailable...>{},
//                 tape, next...);
//         }
//     }
// }

// template <class Input, std::size_t N, class this_type, typename...
// TypesAlive,
//           typename... LeavesAlive, typename... Leaves,
//           std::size_t... IdxTypesAlive, std::size_t... IdxLeavesAlive,
//           std::size_t... IdxAvailable>
// static inline void
// evaluate_univariate(args<LeavesAlive...> const &, args<Leaves...> const &,
//                     std::index_sequence<IdxTypesAlive...> const &,
//                     std::index_sequence<IdxLeavesAlive...> const &,
//                     std::index_sequence<IdxAvailable...> const &,
//                     std::array<double, N> &tape, this_type const &in,
//                     TypesAlive const &...next) {

//     static_assert(!has_type2<this_type, Leaves...>());

//     constexpr bool other_types_depend_on_this =
//         (depends<TypesAlive, this_type>::call() || ...);

//     if constexpr (other_types_depend_on_this) {
//         evaluate_skip(args<LeavesAlive...>{}, args<Leaves...>{},
//                       std::index_sequence<IdxTypesAlive...>{},
//                       std::index_sequence<IdxLeavesAlive...>{},
//                       std::index_sequence<IdxAvailable...>{}, tape, in,
//                       next...);
//     } else {
//         evaluate_univariate_noskip<Input>(
//             args<LeavesAlive...>{}, args<Leaves...>{},
//             std::index_sequence<IdxTypesAlive...>{},
//             std::index_sequence<IdxLeavesAlive...>{},
//             std::index_sequence<IdxAvailable...>{}, tape, in, next...);
//     }
// }

// template <class Input1, class Input2, std::size_t N, class this_type,
//           typename... TypesAlive, typename... LeavesAlive, typename...
//           Leaves, std::size_t IdxTypesAliveCurrent, std::size_t...
//           IdxTypesAlive, std::size_t... IdxLeavesAlive, std::size_t
//           IdxAvailableCurrent, std::size_t... IdxAvailable>
// static inline void evaluate_bivariate_noskip_newloc(
//     args<LeavesAlive...> const &, args<Leaves...> const &,
//     std::index_sequence<IdxTypesAliveCurrent, IdxTypesAlive...> const &,
//     std::index_sequence<IdxLeavesAlive...> const &,
//     std::index_sequence<IdxAvailableCurrent, IdxAvailable...> const &,
//     std::array<double, N> &tape, this_type const &in,
//     TypesAlive const &...next) {

//     constexpr bool input1_has_0_deriv =
//         !equal_or_depends_many<Input1, Leaves...>();
//     constexpr bool input2_has_0_deriv =
//         !equal_or_depends_many<Input2, Leaves...>();
//     static_assert(!input1_has_0_deriv);
//     static_assert(!input2_has_0_deriv);
//     constexpr bool is_input1_new =
//         !has_type2<Input1, TypesAlive..., LeavesAlive...>();
//     constexpr bool is_input2_new =
//         !has_type2<Input2, TypesAlive..., LeavesAlive...>();
//     static_assert(is_input1_new);
//     static_assert(is_input2_new);

//     constexpr bool is_input1_leaf = has_type2<Input1, Leaves...>();
//     constexpr bool is_input2_leaf = has_type2<Input2, Leaves...>();

//     if constexpr (is_input1_leaf && is_input2_leaf) {
// #ifdef CODELOGGER
//         std::cout << "tape[" << IdxAvailableCurrent << "] = tape["
//                   << IdxTypesAliveCurrent << "] * ";
// #endif
//         tape[IdxAvailableCurrent] = tape[IdxTypesAliveCurrent] * in.d2();
// #ifdef CODELOGGER
//         std::cout << ";" << std::endl;
//         std::cout << "tape[" << IdxTypesAliveCurrent << "] *= ";
// #endif
//         tape[IdxTypesAliveCurrent] *= in.d1();
// #ifdef CODELOGGER
//         std::cout << ";" << std::endl;
// #endif
//         evaluate(args<Input1, Input2, LeavesAlive...>{}, args<Leaves...>{},
//                  std::index_sequence<IdxTypesAlive...>{},
//                  std::index_sequence<IdxTypesAliveCurrent,
//                  IdxAvailableCurrent,
//                                      IdxLeavesAlive...>{},
//                  std::index_sequence<IdxAvailable...>{}, tape, next...);
//     } else if constexpr (is_input1_leaf) {
// #ifdef CODELOGGER
//         std::cout << "tape[" << IdxAvailableCurrent << "] = tape["
//                   << IdxTypesAliveCurrent << "] * ";
// #endif
//         tape[IdxAvailableCurrent] = tape[IdxTypesAliveCurrent] * in.d1();
// #ifdef CODELOGGER
//         std::cout << ";" << std::endl;
//         std::cout << "tape[" << IdxTypesAliveCurrent << "] *= ";
// #endif
//         tape[IdxTypesAliveCurrent] *= in.d2();
// #ifdef CODELOGGER
//         std::cout << ";" << std::endl;
// #endif
//         evaluate(args<Input1, LeavesAlive...>{}, args<Leaves...>{},
//                  std::index_sequence<IdxTypesAliveCurrent,
//                  IdxTypesAlive...>{},
//                  std::index_sequence<IdxAvailableCurrent,
//                  IdxLeavesAlive...>{},
//                  std::index_sequence<IdxAvailable...>{}, tape, in.input2(),
//                  next...);
//     } else if constexpr (is_input2_leaf) {
// #ifdef CODELOGGER
//         std::cout << "tape[" << IdxAvailableCurrent << "] = tape["
//                   << IdxTypesAliveCurrent << "] * ";
// #endif
//         tape[IdxAvailableCurrent] = tape[IdxTypesAliveCurrent] * in.d2();
// #ifdef CODELOGGER
//         std::cout << ";" << std::endl;
//         std::cout << "tape[" << IdxTypesAliveCurrent << "] *= ";
// #endif
//         tape[IdxTypesAliveCurrent] *= in.d1();
// #ifdef CODELOGGER
//         std::cout << ";" << std::endl;
// #endif
//         evaluate(args<Input2, LeavesAlive...>{}, args<Leaves...>{},
//                  std::index_sequence<IdxTypesAliveCurrent,
//                  IdxTypesAlive...>{},
//                  std::index_sequence<IdxAvailableCurrent,
//                  IdxLeavesAlive...>{},
//                  std::index_sequence<IdxAvailable...>{}, tape, in.input1(),
//                  next...);
//     } else {
// #ifdef CODELOGGER
//         std::cout << "tape[" << IdxAvailableCurrent << "] = tape["
//                   << IdxTypesAliveCurrent << "] * ";
// #endif
//         tape[IdxAvailableCurrent] = tape[IdxTypesAliveCurrent] * in.d2();
// #ifdef CODELOGGER
//         std::cout << ";" << std::endl;
//         std::cout << "tape[" << IdxTypesAliveCurrent << "] *= ";
// #endif
//         tape[IdxTypesAliveCurrent] *= in.d1();
// #ifdef CODELOGGER
//         std::cout << ";" << std::endl;
// #endif
//         evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
//                  std::index_sequence<IdxTypesAliveCurrent,
//                  IdxAvailableCurrent,
//                                      IdxTypesAlive...>{},
//                  std::index_sequence<IdxLeavesAlive...>{},
//                  std::index_sequence<IdxAvailable...>{}, tape, in.input1(),
//                  in.input2(), next...);
//     }
// }

// template <class Input1, class Input2, std::size_t N, class this_type,
//           typename... TypesAlive, typename... LeavesAlive, typename...
//           Leaves, std::size_t IdxTypesAliveCurrent, std::size_t...
//           IdxTypesAlive, std::size_t... IdxLeavesAlive, std::size_t...
//           IdxAvailable>
// static inline void evaluate_bivariate_noskip(
//     args<LeavesAlive...> const &, args<Leaves...> const &,
//     std::index_sequence<IdxTypesAliveCurrent, IdxTypesAlive...> const &,
//     std::index_sequence<IdxLeavesAlive...> const &,
//     std::index_sequence<IdxAvailable...> const &, std::array<double, N>
//     &tape, this_type const &in, TypesAlive const &...next) {

//     constexpr bool input1_has_0_deriv =
//         !equal_or_depends_many<Input1, Leaves...>();
//     constexpr bool input2_has_0_deriv =
//         !equal_or_depends_many<Input2, Leaves...>();
//     static_assert(!input1_has_0_deriv || !input2_has_0_deriv);

//     if constexpr (!input1_has_0_deriv && !input2_has_0_deriv) {
//         constexpr bool is_input1_leaf = has_type2<Input1, Leaves...>();
//         constexpr bool is_input2_leaf = has_type2<Input2, Leaves...>();
//         if constexpr (is_input1_leaf && is_input2_leaf) {
//             constexpr bool is_input1_new_leaf =
//                 !has_type2<Input1, LeavesAlive...>();
//             constexpr bool is_input2_new_leaf =
//                 !has_type2<Input2, LeavesAlive...>();

//             if constexpr (is_input1_new_leaf && is_input2_new_leaf) {
//                 if constexpr (std::is_same_v<Input1, Input2>) {
// #ifdef CODELOGGER
//                     std::cout << "tape[" << IdxTypesAliveCurrent << "] *= 2 *
//                     ";
// #endif
//                     tape[IdxTypesAliveCurrent] *= 2 * in.d1();
// #ifdef CODELOGGER
//                     std::cout << ";" << std::endl;
// #endif
//                     evaluate(args<Input1, LeavesAlive...>{},
//                     args<Leaves...>{},
//                              std::index_sequence<IdxTypesAlive...>{},
//                              std::index_sequence<IdxTypesAliveCurrent,
//                                                  IdxLeavesAlive...>{},
//                              std::index_sequence<IdxAvailable...>{}, tape,
//                              next...);
//                 } else {
//                     static_assert(sizeof...(IdxAvailable) > 0);
//                     evaluate_bivariate_noskip_newloc<Input1, Input2>(
//                         args<LeavesAlive...>{}, args<Leaves...>{},
//                         std::index_sequence<IdxTypesAliveCurrent,
//                                             IdxTypesAlive...>{},
//                         std::index_sequence<IdxLeavesAlive...>{},
//                         std::index_sequence<IdxAvailable...>{}, tape, in,
//                         next...);
//                 }
//             } else if constexpr (is_input1_new_leaf) {
//                 constexpr auto position = idx_type2<Input2,
//                 LeavesAlive...>(); constexpr auto position_on_tape =
//                     Get<position, IdxLeavesAlive...>::value;
// #ifdef CODELOGGER
//                 std::cout << "tape[" << position_on_tape << "] += tape["
//                           << IdxTypesAliveCurrent << "] * ";
// #endif
//                 tape[position_on_tape] += tape[IdxTypesAliveCurrent] *
//                 in.d2();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
//                 std::cout << "tape[" << IdxTypesAliveCurrent << "] *= ";
// #endif
//                 tape[IdxTypesAliveCurrent] *= in.d1();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
// #endif
//                 evaluate(args<Input1, LeavesAlive...>{}, args<Leaves...>{},
//                          std::index_sequence<IdxTypesAlive...>{},
//                          std::index_sequence<IdxTypesAliveCurrent,
//                                              IdxLeavesAlive...>{},
//                          std::index_sequence<IdxAvailable...>{}, tape,
//                          next...);
//             } else if constexpr (is_input2_new_leaf) {
//                 constexpr auto position = idx_type2<Input1,
//                 LeavesAlive...>(); constexpr auto position_on_tape =
//                     Get<position, IdxLeavesAlive...>::value;
// #ifdef CODELOGGER
//                 std::cout << "tape[" << position_on_tape << "] += tape["
//                           << IdxTypesAliveCurrent << "] * ";
// #endif
//                 tape[position_on_tape] += tape[IdxTypesAliveCurrent] *
//                 in.d1();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
//                 std::cout << "tape[" << IdxTypesAliveCurrent << "] *= ";
// #endif
//                 tape[IdxTypesAliveCurrent] *= in.d2();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
// #endif
//                 evaluate(args<Input2, LeavesAlive...>{}, args<Leaves...>{},
//                          std::index_sequence<IdxTypesAlive...>{},
//                          std::index_sequence<IdxTypesAliveCurrent,
//                                              IdxLeavesAlive...>{},
//                          std::index_sequence<IdxAvailable...>{}, tape,
//                          next...);
//             } else {
//                 if constexpr (std::is_same_v<Input1, Input2>) {
//                     constexpr auto position =
//                         idx_type2<Input1, LeavesAlive...>();
//                     static_assert(position ==
//                                   idx_type2<Input2, LeavesAlive...>());
//                     constexpr auto position_on_tape =
//                         Get<position, IdxLeavesAlive...>::value;
// #ifdef CODELOGGER
//                     std::cout << "tape[" << position_on_tape << "] += 2 *
//                     tape["
//                               << IdxTypesAliveCurrent << "] * ";
// #endif
//                     tape[position_on_tape] +=
//                         2 * tape[IdxTypesAliveCurrent] * in.d1();
// #ifdef CODELOGGER
//                     std::cout << ";" << std::endl;
// #endif
//                     evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
//                              std::index_sequence<IdxTypesAlive...>{},
//                              std::index_sequence<IdxLeavesAlive...>{},
//                              std::index_sequence<IdxTypesAliveCurrent,
//                                                  IdxAvailable...>{},
//                              tape, next...);
//                 } else {
//                     constexpr auto position1 =
//                         idx_type2<Input1, LeavesAlive...>();
//                     constexpr auto position2 =
//                         idx_type2<Input2, LeavesAlive...>();
//                     constexpr auto position1_on_tape =
//                         Get<position1, IdxLeavesAlive...>::value;
//                     constexpr auto position2_on_tape =
//                         Get<position2, IdxLeavesAlive...>::value;
// #ifdef CODELOGGER
//                     std::cout << "tape[" << position1_on_tape << "] += tape["
//                               << IdxTypesAliveCurrent << "] * ";
// #endif
//                     tape[position1_on_tape] +=
//                         tape[IdxTypesAliveCurrent] * in.d1();
// #ifdef CODELOGGER
//                     std::cout << ";" << std::endl;
//                     std::cout << "tape[" << position2_on_tape << "] += tape["
//                               << IdxTypesAliveCurrent << "] * ";
// #endif
//                     tape[position2_on_tape] +=
//                         tape[IdxTypesAliveCurrent] * in.d2();
// #ifdef CODELOGGER
//                     std::cout << ";" << std::endl;
// #endif
//                     evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
//                              std::index_sequence<IdxTypesAlive...>{},
//                              std::index_sequence<IdxLeavesAlive...>{},
//                              std::index_sequence<IdxTypesAliveCurrent,
//                                                  IdxAvailable...>{},
//                              tape, next...);
//                 }
//             }
//         } else if constexpr (is_input1_leaf) {
//             constexpr bool is_input1_new_leaf =
//                 !has_type2<Input1, LeavesAlive...>();
//             constexpr bool is_input2_new = !has_type2<Input2,
//             TypesAlive...>();

//             if constexpr (is_input1_new_leaf && is_input2_new) {
//                 static_assert(sizeof...(IdxAvailable) > 0);
//                 evaluate_bivariate_noskip_newloc<Input1, Input2>(
//                     args<LeavesAlive...>{}, args<Leaves...>{},
//                     std::index_sequence<IdxTypesAliveCurrent,
//                                         IdxTypesAlive...>{},
//                     std::index_sequence<IdxLeavesAlive...>{},
//                     std::index_sequence<IdxAvailable...>{}, tape, in,
//                     next...);
//             } else if constexpr (is_input1_new_leaf) {
//                 constexpr auto position = idx_type2<Input2, TypesAlive...>();
//                 constexpr auto position_on_tape =
//                     Get<position, IdxTypesAlive...>::value;
// #ifdef CODELOGGER
//                 std::cout << "tape[" << position_on_tape << "] += tape["
//                           << IdxTypesAliveCurrent << "] * ";
//                 std::cout << "tape[" << IdxTypesAliveCurrent << "] *=
//                 in.d1();"
//                           << std::endl;
// #endif
//                 tape[position_on_tape] += tape[IdxTypesAliveCurrent] *
//                 in.d2();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
//                 std::cout << "tape[" << IdxTypesAliveCurrent << "] *= ";
// #endif
//                 tape[IdxTypesAliveCurrent] *= in.d1();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
// #endif
//                 evaluate(args<Input1, LeavesAlive...>{}, args<Leaves...>{},
//                          std::index_sequence<IdxTypesAlive...>{},
//                          std::index_sequence<IdxTypesAliveCurrent,
//                                              IdxLeavesAlive...>{},
//                          std::index_sequence<IdxAvailable...>{}, tape,
//                          next...);
//             } else if constexpr (is_input2_new) {
//                 constexpr auto position = idx_type2<Input1,
//                 LeavesAlive...>(); constexpr auto position_on_tape =
//                     Get<position, IdxLeavesAlive...>::value;
// #ifdef CODELOGGER
//                 std::cout << "tape[" << position_on_tape << "] += tape["
//                           << IdxTypesAliveCurrent << "] * ";
// #endif
//                 tape[position_on_tape] += tape[IdxTypesAliveCurrent] *
//                 in.d1();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
//                 std::cout << "tape[" << IdxTypesAliveCurrent << "] *= ";
// #endif
//                 tape[IdxTypesAliveCurrent] *= in.d2();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
// #endif
//                 evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
//                          std::index_sequence<IdxTypesAliveCurrent,
//                                              IdxTypesAlive...>{},
//                          std::index_sequence<IdxLeavesAlive...>{},
//                          std::index_sequence<IdxAvailable...>{}, tape,
//                          in.input2(), next...);
//             } else {
//                 constexpr auto position1 = idx_type2<Input1,
//                 LeavesAlive...>(); constexpr auto position2 =
//                 idx_type2<Input2, TypesAlive...>(); constexpr auto
//                 position1_on_tape =
//                     Get<position1, IdxLeavesAlive...>::value;
//                 constexpr auto position2_on_tape =
//                     Get<position2, IdxTypesAlive...>::value;
// #ifdef CODELOGGER
//                 std::cout << "tape[" << position1_on_tape << "] += tape["
//                           << IdxTypesAliveCurrent << "] * ";
// #endif
//                 tape[position1_on_tape] += tape[IdxTypesAliveCurrent] *
//                 in.d1();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
//                 std::cout << "tape[" << position2_on_tape << "] += tape["
//                           << IdxTypesAliveCurrent << "] * ";
// #endif
//                 tape[position2_on_tape] += tape[IdxTypesAliveCurrent] *
//                 in.d2();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
// #endif
//                 evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
//                          std::index_sequence<IdxTypesAlive...>{},
//                          std::index_sequence<IdxLeavesAlive...>{},
//                          std::index_sequence<IdxTypesAliveCurrent,
//                                              IdxAvailable...>{},
//                          tape, next...);
//             }
//         } else if constexpr (is_input2_leaf) {
//             constexpr bool is_input1_new = !has_type2<Input1,
//             TypesAlive...>(); constexpr bool is_input2_new_leaf =
//                 !has_type2<Input2, LeavesAlive...>();

//             if constexpr (is_input1_new && is_input2_new_leaf) {
//                 static_assert(sizeof...(IdxAvailable) > 0);
//                 evaluate_bivariate_noskip_newloc<Input1, Input2>(
//                     args<LeavesAlive...>{}, args<Leaves...>{},
//                     std::index_sequence<IdxTypesAliveCurrent,
//                                         IdxTypesAlive...>{},
//                     std::index_sequence<IdxLeavesAlive...>{},
//                     std::index_sequence<IdxAvailable...>{}, tape, in,
//                     next...);
//             } else if constexpr (is_input1_new) {
//                 constexpr auto position = idx_type2<Input2,
//                 LeavesAlive...>(); constexpr auto position_on_tape =
//                     Get<position, IdxLeavesAlive...>::value;
// #ifdef CODELOGGER
//                 std::cout << "tape[" << position_on_tape << "] += tape["
//                           << IdxTypesAliveCurrent << "] * ";
// #endif
//                 tape[position_on_tape] += tape[IdxTypesAliveCurrent] *
//                 in.d2();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
//                 std::cout << "tape[" << IdxTypesAliveCurrent << "] *= ";
// #endif
//                 tape[IdxTypesAliveCurrent] *= in.d1();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
// #endif
//                 evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
//                          std::index_sequence<IdxTypesAliveCurrent,
//                                              IdxTypesAlive...>{},
//                          std::index_sequence<IdxLeavesAlive...>{},
//                          std::index_sequence<IdxAvailable...>{}, tape,
//                          in.input1(), next...);
//             } else if constexpr (is_input2_new_leaf) {
//                 constexpr auto position = idx_type2<Input1, TypesAlive...>();
//                 constexpr auto position_on_tape =
//                     Get<position, IdxTypesAlive...>::value;
// #ifdef CODELOGGER
//                 std::cout << "tape[" << position_on_tape << "] += tape["
//                           << IdxTypesAliveCurrent << "] * ";
// #endif
//                 tape[position_on_tape] += tape[IdxTypesAliveCurrent] *
//                 in.d1();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
//                 std::cout << "tape[" << IdxTypesAliveCurrent << "] *= ";
// #endif
//                 tape[IdxTypesAliveCurrent] *= in.d2();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
// #endif
//                 evaluate(args<Input2, LeavesAlive...>{}, args<Leaves...>{},
//                          std::index_sequence<IdxTypesAlive...>{},
//                          std::index_sequence<IdxTypesAliveCurrent,
//                                              IdxLeavesAlive...>{},
//                          std::index_sequence<IdxAvailable...>{}, tape,
//                          next...);
//             } else {
//                 constexpr auto position1 = idx_type2<Input1,
//                 TypesAlive...>(); constexpr auto position2 =
//                 idx_type2<Input2, LeavesAlive...>(); constexpr auto
//                 position1_on_tape =
//                     Get<position1, IdxTypesAlive...>::value;
//                 constexpr auto position2_on_tape =
//                     Get<position2, IdxLeavesAlive...>::value;
// #ifdef CODELOGGER
//                 std::cout << "tape[" << position1_on_tape << "] += tape["
//                           << IdxTypesAliveCurrent << "] * ";
// #endif
//                 tape[position1_on_tape] += tape[IdxTypesAliveCurrent] *
//                 in.d1();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
//                 std::cout << "tape[" << position2_on_tape << "] += tape["
//                           << IdxTypesAliveCurrent << "] * ";
// #endif
//                 tape[position2_on_tape] += tape[IdxTypesAliveCurrent] *
//                 in.d2();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
// #endif
//                 evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
//                          std::index_sequence<IdxTypesAlive...>{},
//                          std::index_sequence<IdxLeavesAlive...>{},
//                          std::index_sequence<IdxTypesAliveCurrent,
//                                              IdxAvailable...>{},
//                          tape, next...);
//             }
//         } else {
//             constexpr bool is_input1_new = !has_type2<Input1,
//             TypesAlive...>(); constexpr bool is_input2_new =
//             !has_type2<Input2, TypesAlive...>();

//             if constexpr (is_input1_new && is_input2_new) {
//                 if constexpr (std::is_same_v<Input1, Input2>) {
// #ifdef CODELOGGER
//                     std::cout << "tape[" << IdxTypesAliveCurrent << "] *= 2 *
//                     ";
// #endif
//                     tape[IdxTypesAliveCurrent] *= 2 * in.d1();
// #ifdef CODELOGGER
//                     std::cout << ";" << std::endl;
// #endif
//                     evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
//                              std::index_sequence<IdxTypesAliveCurrent,
//                                                  IdxTypesAlive...>{},
//                              std::index_sequence<IdxLeavesAlive...>{},
//                              std::index_sequence<IdxAvailable...>{}, tape,
//                              in.input1(), next...);
//                 } else {
//                     static_assert(sizeof...(IdxAvailable) > 0);
//                     evaluate_bivariate_noskip_newloc<Input1, Input2>(
//                         args<LeavesAlive...>{}, args<Leaves...>{},
//                         std::index_sequence<IdxTypesAliveCurrent,
//                                             IdxTypesAlive...>{},
//                         std::index_sequence<IdxLeavesAlive...>{},
//                         std::index_sequence<IdxAvailable...>{}, tape, in,
//                         next...);
//                 }

//             } else if constexpr (is_input1_new) {
//                 constexpr auto position = idx_type2<Input2, TypesAlive...>();
//                 constexpr auto position_on_tape =
//                     Get<position, IdxTypesAlive...>::value;
// #ifdef CODELOGGER
//                 std::cout << "tape[" << position_on_tape << "] += tape["
//                           << IdxTypesAliveCurrent << "] * ";
// #endif
//                 tape[position_on_tape] += tape[IdxTypesAliveCurrent] *
//                 in.d2();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
//                 std::cout << "tape[" << IdxTypesAliveCurrent << "] *= ";
// #endif
//                 tape[IdxTypesAliveCurrent] *= in.d1();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
// #endif
//                 evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
//                          std::index_sequence<IdxTypesAliveCurrent,
//                                              IdxTypesAlive...>{},
//                          std::index_sequence<IdxLeavesAlive...>{},
//                          std::index_sequence<IdxAvailable...>{}, tape,
//                          in.input1(), next...);
//             } else if constexpr (is_input2_new) {
//                 constexpr auto position = idx_type2<Input1, TypesAlive...>();
//                 constexpr auto position_on_tape =
//                     Get<position, IdxTypesAlive...>::value;
// #ifdef CODELOGGER
//                 std::cout << "tape[" << position_on_tape << "] += tape["
//                           << IdxTypesAliveCurrent << "] * ";
// #endif
//                 tape[position_on_tape] += tape[IdxTypesAliveCurrent] *
//                 in.d1();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
//                 std::cout << "tape[" << IdxTypesAliveCurrent << "] *= ";
// #endif
//                 tape[IdxTypesAliveCurrent] *= in.d2();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
// #endif
//                 evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
//                          std::index_sequence<IdxTypesAliveCurrent,
//                                              IdxTypesAlive...>{},
//                          std::index_sequence<IdxLeavesAlive...>{},
//                          std::index_sequence<IdxAvailable...>{}, tape,
//                          in.input2(), next...);
//             } else {
//                 if constexpr (std::is_same_v<Input1, Input2>) {
//                     constexpr auto position =
//                         idx_type2<Input1, TypesAlive...>();
//                     static_assert(position ==
//                                   idx_type2<Input2, TypesAlive...>());
//                     constexpr auto position_on_tape =
//                         Get<position, IdxTypesAlive...>::value;
// #ifdef CODELOGGER
//                     std::cout << "tape[" << position_on_tape << "] += 2 *
//                     tape["
//                               << IdxTypesAliveCurrent << "] * ";
// #endif
//                     tape[position_on_tape] +=
//                         2 * tape[IdxTypesAliveCurrent] * in.d1();
// #ifdef CODELOGGER
//                     std::cout << ";" << std::endl;
// #endif
//                     evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
//                              std::index_sequence<IdxTypesAlive...>{},
//                              std::index_sequence<IdxLeavesAlive...>{},
//                              std::index_sequence<IdxTypesAliveCurrent,
//                                                  IdxAvailable...>{},
//                              tape, next...);
//                 } else {
//                     constexpr auto position1 =
//                         idx_type2<Input1, TypesAlive...>();
//                     constexpr auto position2 =
//                         idx_type2<Input2, TypesAlive...>();
//                     constexpr auto position1_on_tape =
//                         Get<position1, IdxTypesAlive...>::value;
//                     constexpr auto position2_on_tape =
//                         Get<position2, IdxTypesAlive...>::value;
// #ifdef CODELOGGER
//                     std::cout << "tape[" << position1_on_tape << "] += tape["
//                               << IdxTypesAliveCurrent << "] * ";
// #endif
//                     tape[position1_on_tape] +=
//                         tape[IdxTypesAliveCurrent] * in.d1();
// #ifdef CODELOGGER
//                     std::cout << ";" << std::endl;
//                     std::cout << "tape[" << position2_on_tape << "] += tape["
//                               << IdxTypesAliveCurrent << "] * ";
// #endif
//                     tape[position2_on_tape] +=
//                         tape[IdxTypesAliveCurrent] * in.d2();
// #ifdef CODELOGGER
//                     std::cout << ";" << std::endl;
// #endif
//                     evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
//                              std::index_sequence<IdxTypesAlive...>{},
//                              std::index_sequence<IdxLeavesAlive...>{},
//                              std::index_sequence<IdxTypesAliveCurrent,
//                                                  IdxAvailable...>{},
//                              tape, next...);
//                 }
//             }
//         }
//     } else if constexpr (input1_has_0_deriv) {
//         constexpr bool is_input_leaf = has_type2<Input2, Leaves...>();
//         if constexpr (is_input_leaf) {
//             constexpr bool is_input_new_leaf =
//                 !has_type2<Input2, LeavesAlive...>();

//             if constexpr (is_input_new_leaf) {
// #ifdef CODELOGGER
//                 std::cout << "tape[" << IdxTypesAliveCurrent << "] *= ";
// #endif
//                 tape[IdxTypesAliveCurrent] *= in.d2();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
// #endif
//                 evaluate(args<Input2, LeavesAlive...>{}, args<Leaves...>{},
//                          std::index_sequence<IdxTypesAlive...>{},
//                          std::index_sequence<IdxTypesAliveCurrent,
//                                              IdxLeavesAlive...>{},
//                          std::index_sequence<IdxAvailable...>{}, tape,
//                          next...);
//             } else {
//                 constexpr auto position = idx_type2<Input2,
//                 LeavesAlive...>(); constexpr auto position_on_tape =
//                     Get<position, IdxLeavesAlive...>::value;
// #ifdef CODELOGGER
//                 std::cout << "tape[" << position_on_tape << "] += tape["
//                           << IdxTypesAliveCurrent << "] * ";
// #endif
//                 tape[position_on_tape] += tape[IdxTypesAliveCurrent] *
//                 in.d2();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
// #endif
//                 evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
//                          std::index_sequence<IdxTypesAlive...>{},
//                          std::index_sequence<IdxLeavesAlive...>{},
//                          std::index_sequence<IdxTypesAliveCurrent,
//                                              IdxAvailable...>{},
//                          tape, next...);
//             }
//         } else {
//             constexpr bool is_input_new = !has_type2<Input2,
//             TypesAlive...>();

//             if constexpr (is_input_new) {
// #ifdef CODELOGGER
//                 std::cout << "tape[" << IdxTypesAliveCurrent << "] *= ";
// #endif
//                 tape[IdxTypesAliveCurrent] *= in.d2();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
// #endif
//                 evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
//                          std::index_sequence<IdxTypesAliveCurrent,
//                                              IdxTypesAlive...>{},
//                          std::index_sequence<IdxLeavesAlive...>{},
//                          std::index_sequence<IdxAvailable...>{}, tape,
//                          in.input2(), next...);
//             } else {
//                 constexpr auto position = idx_type2<Input2, TypesAlive...>();
//                 constexpr auto position_on_tape =
//                     Get<position, IdxTypesAlive...>::value;
// #ifdef CODELOGGER
//                 std::cout << "tape[" << position_on_tape << "] *= ";
// #endif
//                 tape[position_on_tape] += tape[IdxTypesAliveCurrent] *
//                 in.d2();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
// #endif
//                 evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
//                          std::index_sequence<IdxTypesAlive...>{},
//                          std::index_sequence<IdxLeavesAlive...>{},
//                          std::index_sequence<IdxTypesAliveCurrent,
//                                              IdxAvailable...>{},
//                          tape, next...);
//             }
//         }
//     } else if constexpr (input2_has_0_deriv) {
//         constexpr bool is_input_leaf = has_type2<Input1, Leaves...>();
//         if constexpr (is_input_leaf) {
//             constexpr bool is_input_new_leaf =
//                 !has_type2<Input1, LeavesAlive...>();

//             if constexpr (is_input_new_leaf) {
// #ifdef CODELOGGER
//                 std::cout << "tape[" << IdxTypesAliveCurrent << "] *= ";
// #endif
//                 tape[IdxTypesAliveCurrent] *= in.d1();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
// #endif
//                 evaluate(args<Input1, LeavesAlive...>{}, args<Leaves...>{},
//                          std::index_sequence<IdxTypesAlive...>{},
//                          std::index_sequence<IdxTypesAliveCurrent,
//                                              IdxLeavesAlive...>{},
//                          std::index_sequence<IdxAvailable...>{}, tape,
//                          next...);
//             } else {
//                 constexpr auto position = idx_type2<Input1,
//                 LeavesAlive...>(); constexpr auto position_on_tape =
//                     Get<position, IdxLeavesAlive...>::value;
// #ifdef CODELOGGER
//                 std::cout << "tape[" << position_on_tape << "] += tape["
//                           << IdxTypesAliveCurrent << "] * ";
// #endif
//                 tape[position_on_tape] += tape[IdxTypesAliveCurrent] *
//                 in.d1();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
// #endif
//                 evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
//                          std::index_sequence<IdxTypesAlive...>{},
//                          std::index_sequence<IdxLeavesAlive...>{},
//                          std::index_sequence<IdxTypesAliveCurrent,
//                                              IdxAvailable...>{},
//                          tape, next...);
//             }
//         } else {
//             constexpr bool is_input_new = !has_type2<Input1,
//             TypesAlive...>();

//             if constexpr (is_input_new) {
// #ifdef CODELOGGER
//                 std::cout << "tape[" << IdxTypesAliveCurrent << "] *= ";
// #endif
//                 tape[IdxTypesAliveCurrent] *= in.d1();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
// #endif
//                 evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
//                          std::index_sequence<IdxTypesAliveCurrent,
//                                              IdxTypesAlive...>{},
//                          std::index_sequence<IdxLeavesAlive...>{},
//                          std::index_sequence<IdxAvailable...>{}, tape,
//                          in.input1(), next...);
//             } else {
//                 constexpr auto position = idx_type2<Input1, TypesAlive...>();
//                 constexpr auto position_on_tape =
//                     Get<position, IdxTypesAlive...>::value;
// #ifdef CODELOGGER
//                 std::cout << "tape[" << position_on_tape << "] += tape["
//                           << IdxTypesAliveCurrent << "] * ";
// #endif
//                 tape[position_on_tape] += tape[IdxTypesAliveCurrent] *
//                 in.d1();
// #ifdef CODELOGGER
//                 std::cout << ";" << std::endl;
// #endif
//                 evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
//                          std::index_sequence<IdxTypesAlive...>{},
//                          std::index_sequence<IdxLeavesAlive...>{},
//                          std::index_sequence<IdxTypesAliveCurrent,
//                                              IdxAvailable...>{},
//                          tape, next...);
//             }
//         }
//     }
// }

// template <class Input1, class Input2, std::size_t N, class this_type,
//           typename... TypesAlive, typename... LeavesAlive, typename...
//           Leaves, std::size_t... IdxTypesAlive, std::size_t...
//           IdxLeavesAlive, std::size_t... IdxAvailable>
// static inline void
// evaluate_bivariate(args<LeavesAlive...> const &, args<Leaves...> const &,
//                    std::index_sequence<IdxTypesAlive...> const &,
//                    std::index_sequence<IdxLeavesAlive...> const &,
//                    std::index_sequence<IdxAvailable...> const &,
//                    std::array<double, N> &tape, this_type const &in,
//                    TypesAlive const &...next) {

//     static_assert(!has_type2<this_type, Leaves...>());

//     constexpr bool other_types_depend_on_this =
//         (depends<TypesAlive, this_type>::call() || ...);

//     if constexpr (other_types_depend_on_this) {
//         evaluate_skip(args<LeavesAlive...>{}, args<Leaves...>{},
//                       std::index_sequence<IdxTypesAlive...>{},
//                       std::index_sequence<IdxLeavesAlive...>{},
//                       std::index_sequence<IdxAvailable...>{}, tape, in,
//                       next...);
//     } else {
//         evaluate_bivariate_noskip<Input1, Input2>(
//             args<LeavesAlive...>{}, args<Leaves...>{},
//             std::index_sequence<IdxTypesAlive...>{},
//             std::index_sequence<IdxLeavesAlive...>{},
//             std::index_sequence<IdxAvailable...>{}, tape, in, next...);
//     }
// }

// template <std::size_t N, template <class> class Univariate,
//           typename... TypesAlive, typename... LeavesAlive, typename...
//           Leaves, std::size_t... IdxTypesAlive, std::size_t...
//           IdxLeavesAlive, std::size_t... IdxAvailable, class Input>
// static inline void
// evaluate(args<LeavesAlive...> const &, args<Leaves...> const &,
//          std::index_sequence<IdxTypesAlive...> const &,
//          std::index_sequence<IdxLeavesAlive...> const &,
//          std::index_sequence<IdxAvailable...> const &,
//          std::array<double, N> &tape, Univariate<Input> const &in,
//          TypesAlive const &...next) {
//     static_assert((sizeof...(TypesAlive) + 1) == sizeof...(IdxTypesAlive));
//     static_assert(sizeof...(LeavesAlive) == sizeof...(IdxLeavesAlive));
//     static_assert(sizeof...(IdxTypesAlive) + sizeof...(IdxLeavesAlive) +
//                       sizeof...(IdxAvailable) ==
//                   N);
//     evaluate_univariate<Input>(args<LeavesAlive...>{}, args<Leaves...>{},
//                                std::index_sequence<IdxTypesAlive...>{},
//                                std::index_sequence<IdxLeavesAlive...>{},
//                                std::index_sequence<IdxAvailable...>{}, tape,
//                                in, next...);
// }

// template <std::size_t N, template <class, class> class Bivariate,
//           typename... TypesAlive, typename... LeavesAlive, typename...
//           Leaves, std::size_t... IdxTypesAlive, std::size_t...
//           IdxLeavesAlive, std::size_t... IdxAvailable, class Input1, class
//           Input2>
// static inline void
// evaluate(args<LeavesAlive...> const &, args<Leaves...> const &,
//          std::index_sequence<IdxTypesAlive...> const &,
//          std::index_sequence<IdxLeavesAlive...> const &,
//          std::index_sequence<IdxAvailable...> const &,
//          std::array<double, N> &tape, Bivariate<Input1, Input2> const &in,
//          TypesAlive const &...next) {
//     static_assert((sizeof...(TypesAlive) + 1) == sizeof...(IdxTypesAlive));
//     static_assert(sizeof...(LeavesAlive) == sizeof...(IdxLeavesAlive));
//     static_assert(sizeof...(IdxTypesAlive) + sizeof...(IdxLeavesAlive) +
//                       sizeof...(IdxAvailable) ==
//                   N);
//     evaluate_bivariate<Input1, Input2>(
//         args<LeavesAlive...>{}, args<Leaves...>{},
//         std::index_sequence<IdxTypesAlive...>{},
//         std::index_sequence<IdxLeavesAlive...>{},
//         std::index_sequence<IdxAvailable...>{}, tape, in, next...);
// }

// template <std::size_t N, class this_type, typename... Leaves,
//           std::size_t... IdxTypesAlive, std::size_t... IdxLeavesAlive,
//           std::size_t FirstIdxAvailable, std::size_t... IdxAvailable>
// static inline void evaluate_first(
//     std::array<double, N> &tape, this_type const &in, double init,
//     args<Leaves...> const &,
//     std::index_sequence<FirstIdxAvailable, IdxAvailable...> const &) {
//     tape[FirstIdxAvailable] = init;
// #ifdef CODELOGGER
//     std::cout << "tape[" << FirstIdxAvailable << "] = " << init << ";"
//               << std::endl;
// #endif
//     constexpr bool is_input_leaf = has_type2<this_type, Leaves...>();
//     if constexpr (!is_input_leaf) {
//         evaluate(args<>{}, args<Leaves...>{},
//                  std::index_sequence<FirstIdxAvailable>{},
//                  std::index_sequence<>{},
//                  std::index_sequence<IdxAvailable...>{}, tape, in);
//     }
// }

} // namespace detail

// template <typename... Leaves, class Output>
// constexpr static inline auto evaluate(Output const &in, double init = 1.0)
//     -> std::array<double, sizeof...(Leaves)> {
//     static_assert((depends<Output, Leaves const>::call() && ...));
//     std::array<double, sizeof...(Leaves)> results = {};
//     if constexpr (results.empty()) {
//         return results;
//     } else {
//         constexpr std::size_t size = detail::tape_size(
//             args<Output const>{}, args<>{}, args<Leaves const...>{});
//         std::array<double, size> tape{};
// #ifdef CODELOGGER
//         std::cout << "std::array<double, " << size << "> tape{};" <<
//         std::endl;
// #endif
//         detail::evaluate_first(tape, in, init, args<Leaves const...>{},
//                                std::make_index_sequence<size>{});

//         for (std::size_t index = 0; index < sizeof...(Leaves); ++index) {
//             results[index] = tape[index];
//         }
//         return results;
//     }
// }

// // in case you have the variable around and you don't want to type decltype
// template <typename... Leaves, class Output>
// constexpr static inline auto evaluate(Output const &in, double init,
//                                       Leaves const &.../* leaves */)
//     -> std::array<double, sizeof...(Leaves)> {
//     return evaluate<Leaves...>(in, init);
// }

} // namespace adhoc3

#endif // EVALUATE_BWD_HPP