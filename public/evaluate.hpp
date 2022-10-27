#ifndef EVALUATE_HPP
#define EVALUATE_HPP

#include <adhoc2.hpp>
#include <tape_size.hpp>

namespace adhoc2 {
namespace detail {

template <std::size_t I, std::size_t M, std::size_t N, typename... LeavesAlive,
          std::size_t... IdxLeavesAlive>
static inline void
evaluate_final(std::array<double, M> &temp, std::array<double, N> const &tape,
               args<LeavesAlive...> const &, args<> const &,
               std::index_sequence<IdxLeavesAlive...> const &) {}

template <std::size_t I, std::size_t M, std::size_t N, typename... LeavesAlive,
          typename CurrentLeaf, typename... Leaves,
          std::size_t... IdxLeavesAlive>
static inline void
evaluate_final(std::array<double, M> &temp, std::array<double, N> const &tape,
               args<LeavesAlive...> const &,
               args<CurrentLeaf, Leaves...> const &,
               std::index_sequence<IdxLeavesAlive...> const &) {

    // what if a leaf is never alive? treat separately
    constexpr bool is_leaf_active =
        has_type2<CurrentLeaf const, LeavesAlive...>();

    if constexpr (is_leaf_active) {
        if constexpr (I < M) {
            constexpr auto position =
                idx_type2<CurrentLeaf const, LeavesAlive...>();
            constexpr auto position_on_tape =
                Get<position, IdxLeavesAlive...>::value;
            temp[I] = tape[position_on_tape];
            evaluate_final<I + 1>(temp, tape, args<LeavesAlive...>{},
                                  args<Leaves...>{},
                                  std::index_sequence<IdxLeavesAlive...>{});
        }
    } else {
        evaluate_final<I>(temp, tape, args<LeavesAlive...>{}, args<Leaves...>{},
                          std::index_sequence<IdxLeavesAlive...>{});
    }
}

template <std::size_t N, typename... LeavesAlive, typename... Leaves,
          std::size_t... IdxLeavesAlive, std::size_t... IdxAvailable>
static inline void evaluate(args<LeavesAlive...> const &,
                            args<Leaves...> const &,
                            std::index_sequence<> const &,
                            std::index_sequence<IdxLeavesAlive...> const &,
                            std::index_sequence<IdxAvailable...> const &,
                            std::array<double, N> &tape) {

    // final stage, we order tape for final output
    constexpr std::size_t final_result_size = sizeof...(IdxLeavesAlive);
    std::array<double, final_result_size> temp;
    evaluate_final<0>(temp, tape, args<LeavesAlive...>{}, args<Leaves...>{},
                      std::index_sequence<IdxLeavesAlive...>{});

    for (std::size_t index = 0; index < final_result_size; ++index) {
        tape[index] = temp[index];
    }
}

template <std::size_t N, class FirstType, class SecondType,
          typename... TypesAlive, typename... LeavesAlive, typename... Leaves,
          std::size_t IdxTypesAliveFirst, std::size_t IdxTypesAliveSecond,
          std::size_t... IdxTypesAlive, std::size_t... IdxLeavesAlive,
          std::size_t... IdxAvailable>
static inline void
evaluate_skip(args<LeavesAlive...> const &, args<Leaves...> const &,
              std::index_sequence<IdxTypesAliveFirst, IdxTypesAliveSecond,
                                  IdxTypesAlive...> const &,
              std::index_sequence<IdxLeavesAlive...> const &,
              std::index_sequence<IdxAvailable...> const &,
              std::array<double, N> &tape, FirstType const &first,
              SecondType const &second, TypesAlive const &...next) {

    evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
             std::index_sequence<
                 IdxTypesAliveSecond, IdxTypesAlive...,
                 IdxTypesAliveFirst>{} /* we place first at the end */,
             std::index_sequence<IdxLeavesAlive...>{},
             std::index_sequence<IdxAvailable...>{}, tape, second, next...,
             first /* we place first at the end */);
}

template <class Input, std::size_t N, class this_type, typename... TypesAlive,
          typename... LeavesAlive, typename... Leaves,
          std::size_t IdxTypesAliveCurrent, std::size_t... IdxTypesAlive,
          std::size_t... IdxLeavesAlive, std::size_t... IdxAvailable>
static inline void evaluate_univariate_noskip(
    args<LeavesAlive...> const &, args<Leaves...> const &,
    std::index_sequence<IdxTypesAliveCurrent, IdxTypesAlive...> const &,
    std::index_sequence<IdxLeavesAlive...> const &,
    std::index_sequence<IdxAvailable...> const &, std::array<double, N> &tape,
    this_type const &in, TypesAlive const &...next) {

    constexpr bool is_input_leaf = has_type2<Input, Leaves...>();
    if constexpr (is_input_leaf) {
        constexpr bool is_input_new_leaf = !has_type2<Input, LeavesAlive...>();

        if constexpr (is_input_new_leaf) {
            tape[IdxTypesAliveCurrent] *= in.d();
            evaluate(
                args<Input, LeavesAlive...>{}, args<Leaves...>{},
                std::index_sequence<IdxTypesAlive...>{},
                std::index_sequence<IdxTypesAliveCurrent, IdxLeavesAlive...>{},
                std::index_sequence<IdxAvailable...>{}, tape, next...);
        } else {
            constexpr auto position = idx_type2<Input, LeavesAlive...>();
            tape[position] += tape[IdxTypesAliveCurrent] * in.d();
            evaluate(
                args<LeavesAlive...>{}, args<Leaves...>{},
                std::index_sequence<IdxTypesAlive...>{},
                std::index_sequence<IdxLeavesAlive...>{},
                std::index_sequence<IdxTypesAliveCurrent, IdxAvailable...>{},
                tape, next...);
        }
    } else {
        constexpr bool is_input_new = !has_type2<Input, TypesAlive...>();

        if constexpr (is_input_new) {
            tape[IdxTypesAliveCurrent] *= in.d();
            evaluate(
                args<LeavesAlive...>{}, args<Leaves...>{},
                std::index_sequence<IdxTypesAliveCurrent, IdxTypesAlive...>{},
                std::index_sequence<IdxLeavesAlive...>{},
                std::index_sequence<IdxAvailable...>{}, tape, in.input(),
                next...);
        } else {
            constexpr auto position = idx_type2<Input, TypesAlive...>();
            tape[position] += tape[IdxTypesAliveCurrent] * in.d();
            evaluate(
                args<LeavesAlive...>{}, args<Leaves...>{},
                std::index_sequence<IdxTypesAlive...>{},
                std::index_sequence<IdxLeavesAlive...>{},
                std::index_sequence<IdxTypesAliveCurrent, IdxAvailable...>{},
                tape, next...);
        }
    }
}

template <class Input, std::size_t N, class this_type, typename... TypesAlive,
          typename... LeavesAlive, typename... Leaves,
          std::size_t... IdxTypesAlive, std::size_t... IdxLeavesAlive,
          std::size_t... IdxAvailable>
static inline void
evaluate_univariate(args<LeavesAlive...> const &, args<Leaves...> const &,
                    std::index_sequence<IdxTypesAlive...> const &,
                    std::index_sequence<IdxLeavesAlive...> const &,
                    std::index_sequence<IdxAvailable...> const &,
                    std::array<double, N> &tape, this_type const &in,
                    TypesAlive const &...next) {

    static_assert(!has_type2<this_type, Leaves...>());

    constexpr bool other_types_depend_on_this =
        (TypesAlive::template depends_in<this_type>() || ...);

    if constexpr (other_types_depend_on_this) {
        evaluate_skip(args<LeavesAlive...>{}, args<Leaves...>{},
                      std::index_sequence<IdxTypesAlive...>{},
                      std::index_sequence<IdxLeavesAlive...>{},
                      std::index_sequence<IdxAvailable...>{}, tape, in,
                      next...);
    } else {
        evaluate_univariate_noskip<Input>(
            args<LeavesAlive...>{}, args<Leaves...>{},
            std::index_sequence<IdxTypesAlive...>{},
            std::index_sequence<IdxLeavesAlive...>{},
            std::index_sequence<IdxAvailable...>{}, tape, in, next...);
    }
}

template <class Input1, class Input2, std::size_t N, class this_type,
          typename... TypesAlive, typename... LeavesAlive, typename... Leaves,
          std::size_t IdxTypesAliveCurrent, std::size_t... IdxTypesAlive,
          std::size_t... IdxLeavesAlive, std::size_t IdxAvailableCurrent,
          std::size_t... IdxAvailable>
static inline void evaluate_bivariate_noskip_newloc(
    args<LeavesAlive...> const &, args<Leaves...> const &,
    std::index_sequence<IdxTypesAliveCurrent, IdxTypesAlive...> const &,
    std::index_sequence<IdxLeavesAlive...> const &,
    std::index_sequence<IdxAvailableCurrent, IdxAvailable...> const &,
    std::array<double, N> &tape, this_type const &in,
    TypesAlive const &...next) {

    constexpr bool input1_has_0_deriv =
        !static_cast<bool>(Input1::template depends3<Leaves...>());
    constexpr bool input2_has_0_deriv =
        !static_cast<bool>(Input2::template depends3<Leaves...>());
    static_assert(!input1_has_0_deriv);
    static_assert(!input2_has_0_deriv);
    constexpr bool is_input1_new =
        !has_type2<Input1, TypesAlive..., LeavesAlive...>();
    constexpr bool is_input2_new =
        !has_type2<Input2, TypesAlive..., LeavesAlive...>();
    static_assert(is_input1_new);
    static_assert(is_input2_new);

    constexpr bool is_input1_leaf = has_type2<Input1, Leaves...>();
    constexpr bool is_input2_leaf = has_type2<Input2, Leaves...>();

    if constexpr (is_input1_leaf && is_input2_leaf) {
        tape[IdxAvailableCurrent] = tape[IdxTypesAliveCurrent] * in.d2();
        tape[IdxTypesAliveCurrent] *= in.d1();
        evaluate(args<Input1, Input2, LeavesAlive...>{}, args<Leaves...>{},
                 std::index_sequence<IdxTypesAlive...>{},
                 std::index_sequence<IdxTypesAliveCurrent, IdxAvailableCurrent,
                                     IdxLeavesAlive...>{},
                 std::index_sequence<IdxAvailable...>{}, tape, next...);
    } else if constexpr (is_input1_leaf) {
        tape[IdxAvailableCurrent] = tape[IdxTypesAliveCurrent] * in.d1();
        tape[IdxTypesAliveCurrent] *= in.d2();
        evaluate(args<Input1, LeavesAlive...>{}, args<Leaves...>{},
                 std::index_sequence<IdxTypesAliveCurrent, IdxTypesAlive...>{},
                 std::index_sequence<IdxAvailableCurrent, IdxLeavesAlive...>{},
                 std::index_sequence<IdxAvailable...>{}, tape, in.input2(),
                 next...);
    } else if constexpr (is_input2_leaf) {
        tape[IdxAvailableCurrent] = tape[IdxTypesAliveCurrent] * in.d2();
        tape[IdxTypesAliveCurrent] *= in.d1();
        evaluate(args<Input2, LeavesAlive...>{}, args<Leaves...>{},
                 std::index_sequence<IdxTypesAliveCurrent, IdxTypesAlive...>{},
                 std::index_sequence<IdxAvailableCurrent, IdxLeavesAlive...>{},
                 std::index_sequence<IdxAvailable...>{}, tape, in.input1(),
                 next...);
    } else {
        tape[IdxAvailableCurrent] = tape[IdxTypesAliveCurrent] * in.d2();
        tape[IdxTypesAliveCurrent] *= in.d1();
        evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
                 std::index_sequence<IdxTypesAliveCurrent, IdxAvailableCurrent,
                                     IdxTypesAlive...>{},
                 std::index_sequence<IdxLeavesAlive...>{},
                 std::index_sequence<IdxAvailable...>{}, tape, in.input1(),
                 in.input2(), next...);
    }
}

template <class Input1, class Input2, std::size_t N, class this_type,
          typename... TypesAlive, typename... LeavesAlive, typename... Leaves,
          std::size_t IdxTypesAliveCurrent, std::size_t... IdxTypesAlive,
          std::size_t... IdxLeavesAlive, std::size_t... IdxAvailable>
static inline void evaluate_bivariate_noskip(
    args<LeavesAlive...> const &, args<Leaves...> const &,
    std::index_sequence<IdxTypesAliveCurrent, IdxTypesAlive...> const &,
    std::index_sequence<IdxLeavesAlive...> const &,
    std::index_sequence<IdxAvailable...> const &, std::array<double, N> &tape,
    this_type const &in, TypesAlive const &...next) {

    constexpr bool input1_has_0_deriv =
        !static_cast<bool>(Input1::template depends3<Leaves...>());
    constexpr bool input2_has_0_deriv =
        !static_cast<bool>(Input2::template depends3<Leaves...>());
    static_assert(!input1_has_0_deriv || !input2_has_0_deriv);

    if constexpr (!input1_has_0_deriv && !input2_has_0_deriv) {
        constexpr bool is_input1_leaf = has_type2<Input1, Leaves...>();
        constexpr bool is_input2_leaf = has_type2<Input2, Leaves...>();
        if constexpr (is_input1_leaf && is_input2_leaf) {
            constexpr bool is_input1_new_leaf =
                !has_type2<Input1, LeavesAlive...>();
            constexpr bool is_input2_new_leaf =
                !has_type2<Input2, LeavesAlive...>();

            if constexpr (is_input1_new_leaf && is_input2_new_leaf) {
                if constexpr (std::is_same_v<Input1, Input2>) {
                    tape[IdxTypesAliveCurrent] *= 2 * in.d1();
                    evaluate(args<Input1, LeavesAlive...>{}, args<Leaves...>{},
                             std::index_sequence<IdxTypesAlive...>{},
                             std::index_sequence<IdxTypesAliveCurrent,
                                                 IdxLeavesAlive...>{},
                             std::index_sequence<IdxAvailable...>{}, tape,
                             next...);
                } else {
                    static_assert(sizeof...(IdxAvailable) > 0);
                    evaluate_bivariate_noskip_newloc<Input1, Input2>(
                        args<LeavesAlive...>{}, args<Leaves...>{},
                        std::index_sequence<IdxTypesAliveCurrent,
                                            IdxTypesAlive...>{},
                        std::index_sequence<IdxLeavesAlive...>{},
                        std::index_sequence<IdxAvailable...>{}, tape, in,
                        next...);
                }
            } else if constexpr (is_input1_new_leaf) {
                constexpr auto position = idx_type2<Input2, LeavesAlive...>();
                tape[position] += tape[IdxTypesAliveCurrent] * in.d2();
                tape[IdxTypesAliveCurrent] *= in.d1();
                evaluate(args<Input1, LeavesAlive...>{}, args<Leaves...>{},
                         std::index_sequence<IdxTypesAlive...>{},
                         std::index_sequence<IdxTypesAliveCurrent,
                                             IdxLeavesAlive...>{},
                         std::index_sequence<IdxAvailable...>{}, tape, next...);
            } else if constexpr (is_input2_new_leaf) {
                constexpr auto position = idx_type2<Input1, LeavesAlive...>();
                tape[position] += tape[IdxTypesAliveCurrent] * in.d1();
                tape[IdxTypesAliveCurrent] *= in.d2();
                evaluate(args<Input2, LeavesAlive...>{}, args<Leaves...>{},
                         std::index_sequence<IdxTypesAlive...>{},
                         std::index_sequence<IdxTypesAliveCurrent,
                                             IdxLeavesAlive...>{},
                         std::index_sequence<IdxAvailable...>{}, tape, next...);
            } else {
                constexpr auto position1 = idx_type2<Input1, LeavesAlive...>();
                constexpr auto position2 = idx_type2<Input2, LeavesAlive...>();
                tape[position1] += tape[IdxTypesAliveCurrent] * in.d1();
                tape[position2] += tape[IdxTypesAliveCurrent] * in.d2();
                evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
                         std::index_sequence<IdxTypesAlive...>{},
                         std::index_sequence<IdxLeavesAlive...>{},
                         std::index_sequence<IdxTypesAliveCurrent,
                                             IdxAvailable...>{},
                         tape, next...);
            }
        } else if constexpr (is_input1_leaf) {
            constexpr bool is_input1_new_leaf =
                !has_type2<Input1, LeavesAlive...>();
            constexpr bool is_input2_new = !has_type2<Input2, TypesAlive...>();

            if constexpr (is_input1_new_leaf && is_input2_new) {
                static_assert(sizeof...(IdxAvailable) > 0);
                evaluate_bivariate_noskip_newloc<Input1, Input2>(
                    args<LeavesAlive...>{}, args<Leaves...>{},
                    std::index_sequence<IdxTypesAliveCurrent,
                                        IdxTypesAlive...>{},
                    std::index_sequence<IdxLeavesAlive...>{},
                    std::index_sequence<IdxAvailable...>{}, tape, in, next...);
            } else if constexpr (is_input1_new_leaf) {
                constexpr auto position = idx_type2<Input2, TypesAlive...>();
                tape[position] += tape[IdxTypesAliveCurrent] * in.d2();
                tape[IdxTypesAliveCurrent] *= in.d1();
                evaluate(args<Input1, LeavesAlive...>{}, args<Leaves...>{},
                         std::index_sequence<IdxTypesAlive...>{},
                         std::index_sequence<IdxTypesAliveCurrent,
                                             IdxLeavesAlive...>{},
                         std::index_sequence<IdxAvailable...>{}, tape, next...);
            } else if constexpr (is_input2_new) {
                constexpr auto position = idx_type2<Input1, LeavesAlive...>();
                tape[position] += tape[IdxTypesAliveCurrent] * in.d1();
                tape[IdxTypesAliveCurrent] *= in.d2();
                evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
                         std::index_sequence<IdxTypesAliveCurrent,
                                             IdxTypesAlive...>{},
                         std::index_sequence<IdxLeavesAlive...>{},
                         std::index_sequence<IdxAvailable...>{}, tape,
                         in.input2(), next...);
            } else {
                constexpr auto position1 = idx_type2<Input1, LeavesAlive...>();
                constexpr auto position2 = idx_type2<Input2, TypesAlive...>();
                tape[position1] += tape[IdxTypesAliveCurrent] * in.d1();
                tape[position2] += tape[IdxTypesAliveCurrent] * in.d2();
                evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
                         std::index_sequence<IdxTypesAlive...>{},
                         std::index_sequence<IdxLeavesAlive...>{},
                         std::index_sequence<IdxTypesAliveCurrent,
                                             IdxAvailable...>{},
                         tape, next...);
            }

        } else if constexpr (is_input2_leaf) {
            constexpr bool is_input1_new = !has_type2<Input1, TypesAlive...>();
            constexpr bool is_input2_new_leaf =
                !has_type2<Input2, LeavesAlive...>();

            if constexpr (is_input1_new && is_input2_new_leaf) {
                static_assert(sizeof...(IdxAvailable) > 0);
                evaluate_bivariate_noskip_newloc<Input1, Input2>(
                    args<LeavesAlive...>{}, args<Leaves...>{},
                    std::index_sequence<IdxTypesAliveCurrent,
                                        IdxTypesAlive...>{},
                    std::index_sequence<IdxLeavesAlive...>{},
                    std::index_sequence<IdxAvailable...>{}, tape, in, next...);
            } else if constexpr (is_input1_new) {
                constexpr auto position = idx_type2<Input2, LeavesAlive...>();
                tape[position] += tape[IdxTypesAliveCurrent] * in.d2();
                tape[IdxTypesAliveCurrent] *= in.d1();
                evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
                         std::index_sequence<IdxTypesAliveCurrent,
                                             IdxTypesAlive...>{},
                         std::index_sequence<IdxLeavesAlive...>{},
                         std::index_sequence<IdxAvailable...>{}, tape,
                         in.input1(), next...);
            } else if constexpr (is_input2_new_leaf) {
                constexpr auto position = idx_type2<Input1, TypesAlive...>();
                tape[position] += tape[IdxTypesAliveCurrent] * in.d1();
                tape[IdxTypesAliveCurrent] *= in.d2();
                evaluate(args<Input2, LeavesAlive...>{}, args<Leaves...>{},
                         std::index_sequence<IdxTypesAlive...>{},
                         std::index_sequence<IdxTypesAliveCurrent,
                                             IdxLeavesAlive...>{},
                         std::index_sequence<IdxAvailable...>{}, tape, next...);
            } else {
                constexpr auto position1 = idx_type2<Input1, TypesAlive...>();
                constexpr auto position2 = idx_type2<Input2, LeavesAlive...>();
                tape[position1] += tape[IdxTypesAliveCurrent] * in.d1();
                tape[position2] += tape[IdxTypesAliveCurrent] * in.d2();
                evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
                         std::index_sequence<IdxTypesAlive...>{},
                         std::index_sequence<IdxLeavesAlive...>{},
                         std::index_sequence<IdxTypesAliveCurrent,
                                             IdxAvailable...>{},
                         tape, next...);
            }
        } else {
            constexpr bool is_input1_new = !has_type2<Input1, TypesAlive...>();
            constexpr bool is_input2_new = !has_type2<Input2, TypesAlive...>();

            if constexpr (is_input1_new && is_input2_new) {
                if constexpr (std::is_same_v<Input1, Input2>) {
                    tape[IdxTypesAliveCurrent] *= 2 * in.d1();
                    evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
                             std::index_sequence<IdxTypesAliveCurrent,
                                                 IdxTypesAlive...>{},
                             std::index_sequence<IdxLeavesAlive...>{},
                             std::index_sequence<IdxAvailable...>{}, tape,
                             in.input1(), next...);
                } else {
                    static_assert(sizeof...(IdxAvailable) > 0);
                    evaluate_bivariate_noskip_newloc<Input1, Input2>(
                        args<LeavesAlive...>{}, args<Leaves...>{},
                        std::index_sequence<IdxTypesAliveCurrent,
                                            IdxTypesAlive...>{},
                        std::index_sequence<IdxLeavesAlive...>{},
                        std::index_sequence<IdxAvailable...>{}, tape, in,
                        next...);
                }

            } else if constexpr (is_input1_new) {
                constexpr auto position = idx_type2<Input2, LeavesAlive...>();
                tape[position] += tape[IdxTypesAliveCurrent] * in.d2();
                tape[IdxTypesAliveCurrent] *= in.d1();
                evaluate(args<TypesAlive...>{}, args<Leaves...>{},
                         std::index_sequence<IdxTypesAliveCurrent,
                                             IdxTypesAlive...>{},
                         std::index_sequence<IdxLeavesAlive...>{},
                         std::index_sequence<IdxAvailable...>{}, tape,
                         in.input1(), next...);
            } else if constexpr (is_input2_new) {
                constexpr auto position = idx_type2<Input1, TypesAlive...>();
                tape[position] += tape[IdxTypesAliveCurrent] * in.d1();
                tape[IdxTypesAliveCurrent] *= in.d2();
                evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
                         std::index_sequence<IdxTypesAliveCurrent,
                                             IdxTypesAlive...>{},
                         std::index_sequence<IdxLeavesAlive...>{},
                         std::index_sequence<IdxAvailable...>{}, tape,
                         in.input2(), next...);
            } else {
                constexpr auto position1 = idx_type2<Input1, TypesAlive...>();
                constexpr auto position2 = idx_type2<Input2, TypesAlive...>();
                tape[position1] += tape[IdxTypesAliveCurrent] * in.d1();
                tape[position2] += tape[IdxTypesAliveCurrent] * in.d2();
                evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
                         std::index_sequence<IdxTypesAlive...>{},
                         std::index_sequence<IdxLeavesAlive...>{},
                         std::index_sequence<IdxTypesAliveCurrent,
                                             IdxAvailable...>{},
                         tape, next...);
            }
        }
    } else if constexpr (input1_has_0_deriv) {
        constexpr bool is_input_leaf = has_type2<Input2, Leaves...>();
        if constexpr (is_input_leaf) {
            constexpr bool is_input_new_leaf =
                !has_type2<Input2, LeavesAlive...>();

            if constexpr (is_input_new_leaf) {
                tape[IdxTypesAliveCurrent] *= in.d2();
                evaluate(args<Input2, LeavesAlive...>{}, args<Leaves...>{},
                         std::index_sequence<IdxTypesAlive...>{},
                         std::index_sequence<IdxTypesAliveCurrent,
                                             IdxLeavesAlive...>{},
                         std::index_sequence<IdxAvailable...>{}, tape, next...);
            } else {
                constexpr auto position = idx_type2<Input2, LeavesAlive...>();
                tape[position] += tape[IdxTypesAliveCurrent] * in.d2();
                evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
                         std::index_sequence<IdxTypesAlive...>{},
                         std::index_sequence<IdxLeavesAlive...>{},
                         std::index_sequence<IdxTypesAliveCurrent,
                                             IdxAvailable...>{},
                         tape, next...);
            }
        } else {
            constexpr bool is_input_new = !has_type2<Input2, TypesAlive...>();

            if constexpr (is_input_new) {
                tape[IdxTypesAliveCurrent] *= in.d2();
                evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
                         std::index_sequence<IdxTypesAliveCurrent,
                                             IdxTypesAlive...>{},
                         std::index_sequence<IdxLeavesAlive...>{},
                         std::index_sequence<IdxAvailable...>{}, tape,
                         in.input2(), next...);
            } else {
                constexpr auto position = idx_type2<Input2, TypesAlive...>();
                tape[position] += tape[IdxTypesAliveCurrent] * in.d2();
                evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
                         std::index_sequence<IdxTypesAlive...>{},
                         std::index_sequence<IdxLeavesAlive...>{},
                         std::index_sequence<IdxTypesAliveCurrent,
                                             IdxAvailable...>{},
                         tape, next...);
            }
        }
    } else if constexpr (input2_has_0_deriv) {
        constexpr bool is_input_leaf = has_type2<Input1, Leaves...>();
        if constexpr (is_input_leaf) {
            constexpr bool is_input_new_leaf =
                !has_type2<Input1, LeavesAlive...>();

            if constexpr (is_input_new_leaf) {
                tape[IdxTypesAliveCurrent] *= in.d1();
                evaluate(args<Input1, LeavesAlive...>{}, args<Leaves...>{},
                         std::index_sequence<IdxTypesAlive...>{},
                         std::index_sequence<IdxTypesAliveCurrent,
                                             IdxLeavesAlive...>{},
                         std::index_sequence<IdxAvailable...>{}, tape, next...);
            } else {
                constexpr auto position = idx_type2<Input1, LeavesAlive...>();
                tape[position] += tape[IdxTypesAliveCurrent] * in.d1();
                evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
                         std::index_sequence<IdxTypesAlive...>{},
                         std::index_sequence<IdxLeavesAlive...>{},
                         std::index_sequence<IdxTypesAliveCurrent,
                                             IdxAvailable...>{},
                         tape, next...);
            }
        } else {
            constexpr bool is_input_new = !has_type2<Input1, TypesAlive...>();

            if constexpr (is_input_new) {
                tape[IdxTypesAliveCurrent] *= in.d1();
                evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
                         std::index_sequence<IdxTypesAliveCurrent,
                                             IdxTypesAlive...>{},
                         std::index_sequence<IdxLeavesAlive...>{},
                         std::index_sequence<IdxAvailable...>{}, tape,
                         in.input1(), next...);
            } else {
                constexpr auto position = idx_type2<Input1, TypesAlive...>();
                tape[position] += tape[IdxTypesAliveCurrent] * in.d1();
                evaluate(args<LeavesAlive...>{}, args<Leaves...>{},
                         std::index_sequence<IdxTypesAlive...>{},
                         std::index_sequence<IdxLeavesAlive...>{},
                         std::index_sequence<IdxTypesAliveCurrent,
                                             IdxAvailable...>{},
                         tape, next...);
            }
        }
    }
}

template <class Input1, class Input2, std::size_t N, class this_type,
          typename... TypesAlive, typename... LeavesAlive, typename... Leaves,
          std::size_t... IdxTypesAlive, std::size_t... IdxLeavesAlive,
          std::size_t... IdxAvailable>
static inline void
evaluate_bivariate(args<LeavesAlive...> const &, args<Leaves...> const &,
                   std::index_sequence<IdxTypesAlive...> const &,
                   std::index_sequence<IdxLeavesAlive...> const &,
                   std::index_sequence<IdxAvailable...> const &,
                   std::array<double, N> &tape, this_type const &in,
                   TypesAlive const &...next) {

    static_assert(!has_type2<this_type, Leaves...>());

    constexpr bool other_types_depend_on_this =
        (TypesAlive::template depends_in<this_type>() || ...);

    if constexpr (other_types_depend_on_this) {
        evaluate_skip(args<LeavesAlive...>{}, args<Leaves...>{},
                      std::index_sequence<IdxTypesAlive...>{},
                      std::index_sequence<IdxLeavesAlive...>{},
                      std::index_sequence<IdxAvailable...>{}, tape, in,
                      next...);
    } else {
        evaluate_bivariate_noskip<Input1, Input2>(
            args<LeavesAlive...>{}, args<Leaves...>{},
            std::index_sequence<IdxTypesAlive...>{},
            std::index_sequence<IdxLeavesAlive...>{},
            std::index_sequence<IdxAvailable...>{}, tape, in, next...);
    }
}

template <std::size_t N, class this_type, typename... TypesAlive,
          typename... LeavesAlive, typename... Leaves,
          std::size_t... IdxTypesAlive, std::size_t... IdxLeavesAlive,
          std::size_t... IdxAvailable, class Input>
static inline void
evaluate(args<LeavesAlive...> const &, args<Leaves...> const &,
         std::index_sequence<IdxTypesAlive...> const &,
         std::index_sequence<IdxLeavesAlive...> const &,
         std::index_sequence<IdxAvailable...> const &,
         std::array<double, N> &tape, Univariate<Input, this_type> const &in,
         TypesAlive const &...next) {
    static_assert((sizeof...(TypesAlive) + 1) == sizeof...(IdxTypesAlive));
    static_assert(sizeof...(LeavesAlive) == sizeof...(IdxLeavesAlive));
    static_assert(sizeof...(IdxTypesAlive) + sizeof...(IdxLeavesAlive) +
                      sizeof...(IdxAvailable) ==
                  N);
    evaluate_univariate<Input>(args<LeavesAlive...>{}, args<Leaves...>{},
                               std::index_sequence<IdxTypesAlive...>{},
                               std::index_sequence<IdxLeavesAlive...>{},
                               std::index_sequence<IdxAvailable...>{}, tape,
                               *static_cast<this_type const *>(&in), next...);
}

// template <class Input1, class Input2, std::size_t N, class this_type,
//           typename... TypesAlive, typename... LeavesAlive, typename...
//           Leaves, std::size_t... IdxTypesAlive, std::size_t...
//           IdxLeavesAlive, std::size_t... IdxAvailable>
// static inline void
// evaluate(std::array<double, N> &tape,
//          Bivariate<Input1, Input2, this_type> const &in,
//          TypesAlive const &...next, args<LeavesAlive...> const &,
//          args<Leaves...> const &, std::index_sequence<IdxTypesAlive...> const
//          &, std::index_sequence<IdxLeavesAlive...> const &,
//          std::index_sequence<IdxAvailable...> const &) {
//     static_assert((sizeof...(TypesAlive) + 1) == sizeof...(IdxTypesAlive));
//     static_assert(sizeof...(LeavesAlive) == sizeof...(IdxLeavesAlive));
//     static_assert(sizeof...(IdxTypesAlive) + sizeof...(IdxLeavesAlive) +
//                       sizeof...(IdxAvailable) ==
//                   N);
//     evaluate_bivariate<Input1, Input2>(
//         tape, *static_cast<this_type const *>(&in), next...,
//         args<LeavesAlive...>{}, args<Leaves...>{},
//         std::index_sequence<IdxTypesAlive...>{},
//         std::index_sequence<IdxLeavesAlive...>{},
//         std::index_sequence<IdxAvailable...>{});
// }

template <std::size_t N, typename... TypesAlive, typename... LeavesAlive,
          typename... Leaves, std::size_t... IdxTypesAlive,
          std::size_t... IdxLeavesAlive, std::size_t... IdxAvailable,
          class Input1, class Input2>
static inline void
evaluate(args<LeavesAlive...> const &, args<Leaves...> const &,
         std::index_sequence<IdxTypesAlive...> const &,
         std::index_sequence<IdxLeavesAlive...> const &,
         std::index_sequence<IdxAvailable...> const &,
         std::array<double, N> &tape, mul<Input1, Input2> const &in,
         TypesAlive const &...next) {
    static_assert((sizeof...(TypesAlive) + 1) == sizeof...(IdxTypesAlive));
    static_assert(sizeof...(LeavesAlive) == sizeof...(IdxLeavesAlive));
    static_assert(sizeof...(IdxTypesAlive) + sizeof...(IdxLeavesAlive) +
                      sizeof...(IdxAvailable) ==
                  N);
    evaluate_bivariate<Input1, Input2>(
        args<LeavesAlive...>{}, args<Leaves...>{},
        std::index_sequence<IdxTypesAlive...>{},
        std::index_sequence<IdxLeavesAlive...>{},
        std::index_sequence<IdxAvailable...>{}, tape, in, next...);
}

template <std::size_t N, typename... TypesAlive, typename... LeavesAlive,
          typename... Leaves, std::size_t... IdxTypesAlive,
          std::size_t... IdxLeavesAlive, std::size_t... IdxAvailable,
          class Input1, class Input2>
static inline void
evaluate(args<LeavesAlive...> const &, args<Leaves...> const &,
         std::index_sequence<IdxTypesAlive...> const &,
         std::index_sequence<IdxLeavesAlive...> const &,
         std::index_sequence<IdxAvailable...> const &,
         std::array<double, N> &tape, add<Input1, Input2> const &in,
         TypesAlive const &...next) {
    static_assert((sizeof...(TypesAlive) + 1) == sizeof...(IdxTypesAlive));
    static_assert(sizeof...(LeavesAlive) == sizeof...(IdxLeavesAlive));
    static_assert(sizeof...(IdxTypesAlive) + sizeof...(IdxLeavesAlive) +
                      sizeof...(IdxAvailable) ==
                  N);
    evaluate_bivariate<Input1, Input2>(
        args<LeavesAlive...>{}, args<Leaves...>{},
        std::index_sequence<IdxTypesAlive...>{},
        std::index_sequence<IdxLeavesAlive...>{},
        std::index_sequence<IdxAvailable...>{}, tape, in, next...);
}

template <std::size_t N, class this_type, typename... Leaves,
          std::size_t... IdxTypesAlive, std::size_t... IdxLeavesAlive,
          std::size_t FirstIdxAvailable, std::size_t... IdxAvailable>
static inline void evaluate_first(
    std::array<double, N> &tape, this_type const &in, args<Leaves...> const &,
    std::index_sequence<FirstIdxAvailable, IdxAvailable...> const &) {
    tape[FirstIdxAvailable] = 1.0;
    evaluate(args<>{}, args<Leaves...>{},
             std::index_sequence<FirstIdxAvailable>{}, std::index_sequence<>{},
             std::index_sequence<IdxAvailable...>{}, tape, in);
}

} // namespace detail

template <typename... Leaves, class Output>
constexpr static inline auto evaluate(Output const &in)
    -> std::array<double, (Output::template depends<Leaves>() + ...)> {

    std::array<double, (Output::template depends<Leaves>() + ...)> results = {};
    if constexpr (results.empty()) {
        return results;
    } else {
        constexpr std::size_t size = detail::tape_size(
            args<Output const>{}, args<>{}, args<Leaves...>{});
        std::array<double, size> tape{};
        detail::evaluate_first(tape, in, args<Leaves...>{},
                               std::make_index_sequence<size>{});

        for (std::size_t index = 0;
             index < (Output::template depends<Leaves>() + ...); ++index) {
            results[index] = tape[index];
        }
        return results;
    }
}

// in case you have the variable around and you don't want to type decltype
template <typename... Leaves, class Output>
constexpr static inline auto evaluate(Output const &in,
                                      Leaves const &.../* leaves */)
    -> std::array<double, (Output::template depends<Leaves>() + ...)> {
    return evaluate<Leaves...>(in);
}

} // namespace adhoc2

#endif // EVALUATE_HPP