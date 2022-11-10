#ifndef NEW_EVALUATE_BWD_HPP
#define NEW_EVALUATE_BWD_HPP

#include "dependency.hpp"
#include "tape_size_bwd.hpp"
#include "tape_static.hpp"

namespace adhoc {
namespace detail {

template <class Val, class... InputTypes, class... IntermediateTypes>
inline auto get2(TapeIntermediate<InputTypes...> const &in,
                 TapeIntermediate<IntermediateTypes...> const &intermediate)
    -> double;

template <std::uint64_t D, class... InputTypes, class... IntermediateTypes>
inline auto get2(TapeIntermediate<InputTypes...> const &,
                 TapeIntermediate<IntermediateTypes...> const &,
                 constants::CD<D> &) -> double {
    return constants::CD<D>::v();
}

template <std::uint64_t D, class... InputTypes, class... IntermediateTypes>
inline auto get2(TapeIntermediate<InputTypes...> const &,
                 TapeIntermediate<IntermediateTypes...> const &,
                 constants::CD<D> const &) -> double {
    return constants::CD<D>::v();
}

template <class Input, class... InputTypes, class... IntermediateTypes>
inline auto get2(TapeIntermediate<InputTypes...> const &in,
                 TapeIntermediate<IntermediateTypes...> const &intermediate,
                 Input &) -> double {
    constexpr bool is_input = has_type2<Input, InputTypes...>();
    if constexpr (is_input) {
        return in.get(Input{});
    } else {
        return intermediate.get(Input{});
    }
}

template <std::size_t I, std::size_t M, std::size_t N, typename... LeavesAlive,
          std::size_t... IdxLeavesAlive>
static inline void
evaluate_final(std::array<double, M> & /* temp */,
               std::array<double, N> const & /* tape */,
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

template <class... InputTypes, class... IntermediateTypes, std::size_t N,
          typename... LeavesAlive, typename... Leaves,
          std::size_t... IdxLeavesAlive, std::size_t... IdxAvailable>
inline void
evaluate_bwd(std::array<double, N> &tape,
             TapeIntermediate<InputTypes...> const & /* in */,
             TapeIntermediate<IntermediateTypes...> const & /* intermediate */,
             args<> const &, args<LeavesAlive...> const &,
             args<Leaves...> const &, std::index_sequence<> const &,
             std::index_sequence<IdxLeavesAlive...> const &,
             std::index_sequence<IdxAvailable...> const &) {

    // final stage, we order tape for final output
    constexpr std::size_t final_result_size = sizeof...(IdxLeavesAlive);
    std::array<double, final_result_size> temp;
    evaluate_final<0>(temp, tape, args<LeavesAlive...>{}, args<Leaves...>{},
                      std::index_sequence<IdxLeavesAlive...>{});

    for (std::size_t index = 0; index < final_result_size; ++index) {
        tape[index] = temp[index];
    }
}

template <class... InputTypes, class... IntermediateTypes, std::size_t N,
          template <class, class> class Bivariate, class Input1, class Input2,
          typename... TypesAlive, typename... LeavesAlive, typename... Leaves,
          std::size_t IdxTypesAliveCurrent, std::size_t... IdxTypesAlive,
          std::size_t... IdxLeavesAlive, std::size_t IdxAvailableCurrent,
          std::size_t... IdxAvailable>
inline void evaluate_bwd_bivariate_noskip_new_loc(
    std::array<double, N> &tape, TapeIntermediate<InputTypes...> const &in,
    TapeIntermediate<IntermediateTypes...> const &intermediate,
    args<Bivariate<Input1, Input2> const, TypesAlive...> const &,
    args<LeavesAlive...> const &, args<Leaves...> const &,
    std::index_sequence<IdxTypesAliveCurrent, IdxTypesAlive...> const &,
    std::index_sequence<IdxLeavesAlive...> const &,
    std::index_sequence<IdxAvailableCurrent, IdxAvailable...> const &) {

    using this_type = Bivariate<Input1, Input2> const;
    constexpr bool input1_has_0_deriv =
        !equal_or_depends_many<Input1, Leaves...>();
    constexpr bool input2_has_0_deriv =
        !equal_or_depends_many<Input2, Leaves...>();
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

    double const d1 = this_type::d1(get2(in, intermediate, this_type{}),
                                    get2(in, intermediate, Input1{}),
                                    get2(in, intermediate, Input2{}));
    double const d2 = this_type::d2(get2(in, intermediate, this_type{}),
                                    get2(in, intermediate, Input1{}),
                                    get2(in, intermediate, Input2{}));

    if constexpr (is_input1_leaf && is_input2_leaf) {
        tape[IdxAvailableCurrent] = tape[IdxTypesAliveCurrent] * d2;
        tape[IdxTypesAliveCurrent] *= d1;
        evaluate_bwd(
            tape, in, intermediate, args<TypesAlive...>{},
            args<Input1, Input2, LeavesAlive...>{}, args<Leaves...>{},
            std::index_sequence<IdxTypesAlive...>{},
            std::index_sequence<IdxTypesAliveCurrent, IdxAvailableCurrent,
                                IdxLeavesAlive...>{},
            std::index_sequence<IdxAvailable...>{});
    } else if constexpr (is_input1_leaf) {
        tape[IdxAvailableCurrent] = tape[IdxTypesAliveCurrent] * d1;
        tape[IdxTypesAliveCurrent] *= d2;
        evaluate_bwd(
            tape, in, intermediate, args<Input2, TypesAlive...>{},
            args<Input1, LeavesAlive...>{}, args<Leaves...>{},
            std::index_sequence<IdxTypesAliveCurrent, IdxTypesAlive...>{},
            std::index_sequence<IdxAvailableCurrent, IdxLeavesAlive...>{},
            std::index_sequence<IdxAvailable...>{});
    } else if constexpr (is_input2_leaf) {
        tape[IdxAvailableCurrent] = tape[IdxTypesAliveCurrent] * d2;
        tape[IdxTypesAliveCurrent] *= d1;
        evaluate_bwd(
            tape, in, intermediate, args<Input1, TypesAlive...>{},
            args<Input2, LeavesAlive...>{}, args<Leaves...>{},
            std::index_sequence<IdxTypesAliveCurrent, IdxTypesAlive...>{},
            std::index_sequence<IdxAvailableCurrent, IdxLeavesAlive...>{},
            std::index_sequence<IdxAvailable...>{});
    } else {
        tape[IdxAvailableCurrent] = tape[IdxTypesAliveCurrent] * d2;
        tape[IdxTypesAliveCurrent] *= d1;
        evaluate_bwd(
            tape, in, intermediate, args<Input1, Input2, TypesAlive...>{},
            args<LeavesAlive...>{}, args<Leaves...>{},
            std::index_sequence<IdxTypesAliveCurrent, IdxAvailableCurrent,
                                IdxTypesAlive...>{},
            std::index_sequence<IdxLeavesAlive...>{},
            std::index_sequence<IdxAvailable...>{});
    }
}

template <class... InputTypes, class... IntermediateTypes, std::size_t N,
          template <class, class> class Bivariate, class Input1, class Input2,
          typename... TypesAlive, typename... LeavesAlive, typename... Leaves,
          std::size_t IdxTypesAliveCurrent, std::size_t... IdxTypesAlive,
          std::size_t... IdxLeavesAlive, std::size_t... IdxAvailable>
inline void evaluate_bwd_bivariate_noskip(
    std::array<double, N> &tape, TapeIntermediate<InputTypes...> const &in,
    TapeIntermediate<IntermediateTypes...> const &intermediate,
    args<Bivariate<Input1, Input2> const, TypesAlive...> const &,
    args<LeavesAlive...> const &, args<Leaves...> const &,
    std::index_sequence<IdxTypesAliveCurrent, IdxTypesAlive...> const &,
    std::index_sequence<IdxLeavesAlive...> const &,
    std::index_sequence<IdxAvailable...> const &) {

    using this_type = Bivariate<Input1, Input2> const;
    constexpr bool input1_has_0_deriv =
        !equal_or_depends_many<Input1, Leaves...>();
    constexpr bool input2_has_0_deriv =
        !equal_or_depends_many<Input2, Leaves...>();
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
                    double const d1 =
                        this_type::d1(get2(in, intermediate, this_type{}),
                                      get2(in, intermediate, Input1{}),
                                      get2(in, intermediate, Input2{}));
                    tape[IdxTypesAliveCurrent] *= 2 * d1;
                    evaluate_bwd(tape, in, intermediate, args<TypesAlive...>{},
                                 args<Input1, LeavesAlive...>{},
                                 args<Leaves...>{},
                                 std::index_sequence<IdxTypesAlive...>{},
                                 std::index_sequence<IdxTypesAliveCurrent,
                                                     IdxLeavesAlive...>{},
                                 std::index_sequence<IdxAvailable...>{});
                } else {
                    static_assert(sizeof...(IdxAvailable) > 0);
                    evaluate_bwd_bivariate_noskip_new_loc(
                        tape, in, intermediate,
                        args<this_type, TypesAlive...>{},
                        args<LeavesAlive...>{}, args<Leaves...>{},
                        std::index_sequence<IdxTypesAliveCurrent,
                                            IdxTypesAlive...>{},
                        std::index_sequence<IdxLeavesAlive...>{},
                        std::index_sequence<IdxAvailable...>{});
                }
            } else if constexpr (is_input1_new_leaf) {
                constexpr auto position = idx_type2<Input2, LeavesAlive...>();
                constexpr auto position_on_tape =
                    Get<position, IdxLeavesAlive...>::value;
                double const d1 =
                    this_type::d1(get2(in, intermediate, this_type{}),
                                  get2(in, intermediate, Input1{}),
                                  get2(in, intermediate, Input2{}));
                double const d2 =
                    this_type::d2(get2(in, intermediate, this_type{}),
                                  get2(in, intermediate, Input1{}),
                                  get2(in, intermediate, Input2{}));
                tape[position_on_tape] += tape[IdxTypesAliveCurrent] * d2;
                tape[IdxTypesAliveCurrent] *= d1;
                evaluate_bwd(tape, in, intermediate, args<TypesAlive...>{},
                             args<Input1, LeavesAlive...>{}, args<Leaves...>{},
                             std::index_sequence<IdxTypesAlive...>{},
                             std::index_sequence<IdxTypesAliveCurrent,
                                                 IdxLeavesAlive...>{},
                             std::index_sequence<IdxAvailable...>{});
            } else if constexpr (is_input2_new_leaf) {
                constexpr auto position = idx_type2<Input1, LeavesAlive...>();
                constexpr auto position_on_tape =
                    Get<position, IdxLeavesAlive...>::value;
                tape[position_on_tape] += tape[IdxTypesAliveCurrent] * in.d1();
                tape[IdxTypesAliveCurrent] *= in.d2();
                evaluate_bwd(tape, in, intermediate, args<TypesAlive...>{},
                             args<Input2, LeavesAlive...>{}, args<Leaves...>{},
                             std::index_sequence<IdxTypesAlive...>{},
                             std::index_sequence<IdxTypesAliveCurrent,
                                                 IdxLeavesAlive...>{},
                             std::index_sequence<IdxAvailable...>{});
            } else {
                if constexpr (std::is_same_v<Input1, Input2>) {
                    constexpr auto position =
                        idx_type2<Input1, LeavesAlive...>();
                    static_assert(position ==
                                  idx_type2<Input2, LeavesAlive...>());
                    constexpr auto position_on_tape =
                        Get<position, IdxLeavesAlive...>::value;
                    double const d1 =
                        this_type::d1(get2(in, intermediate, this_type{}),
                                      get2(in, intermediate, Input1{}),
                                      get2(in, intermediate, Input2{}));
                    tape[position_on_tape] +=
                        2 * tape[IdxTypesAliveCurrent] * d1;
                } else {
                    constexpr auto position1 =
                        idx_type2<Input1, LeavesAlive...>();
                    constexpr auto position2 =
                        idx_type2<Input2, LeavesAlive...>();
                    constexpr auto position1_on_tape =
                        Get<position1, IdxLeavesAlive...>::value;
                    constexpr auto position2_on_tape =
                        Get<position2, IdxLeavesAlive...>::value;
                    double const d1 =
                        this_type::d1(get2(in, intermediate, this_type{}),
                                      get2(in, intermediate, Input1{}),
                                      get2(in, intermediate, Input2{}));
                    double const d2 =
                        this_type::d2(get2(in, intermediate, this_type{}),
                                      get2(in, intermediate, Input1{}),
                                      get2(in, intermediate, Input2{}));
                    tape[position1_on_tape] += tape[IdxTypesAliveCurrent] * d1;
                    tape[position2_on_tape] += tape[IdxTypesAliveCurrent] * d2;
                }
                evaluate_bwd(tape, in, intermediate, args<TypesAlive...>{},
                             args<LeavesAlive...>{}, args<Leaves...>{},
                             std::index_sequence<IdxTypesAlive...>{},
                             std::index_sequence<IdxLeavesAlive...>{},
                             std::index_sequence<IdxTypesAliveCurrent,
                                                 IdxAvailable...>{});
            }
        } else if constexpr (is_input1_leaf) {
            constexpr bool is_input1_new_leaf =
                !has_type2<Input1, LeavesAlive...>();
            constexpr bool is_input2_new = !has_type2<Input2, TypesAlive...>();

            if constexpr (is_input1_new_leaf && is_input2_new) {
                static_assert(sizeof...(IdxAvailable) > 0);
                evaluate_bwd_bivariate_noskip_new_loc(
                    tape, in, intermediate, args<this_type, TypesAlive...>{},
                    args<LeavesAlive...>{}, args<Leaves...>{},
                    std::index_sequence<IdxTypesAliveCurrent,
                                        IdxTypesAlive...>{},
                    std::index_sequence<IdxLeavesAlive...>{},
                    std::index_sequence<IdxAvailable...>{});
            } else if constexpr (is_input1_new_leaf) {
                constexpr auto position = idx_type2<Input2, TypesAlive...>();
                constexpr auto position_on_tape =
                    Get<position, IdxTypesAlive...>::value;
                double const d1 =
                    this_type::d1(get2(in, intermediate, this_type{}),
                                  get2(in, intermediate, Input1{}),
                                  get2(in, intermediate, Input2{}));
                double const d2 =
                    this_type::d2(get2(in, intermediate, this_type{}),
                                  get2(in, intermediate, Input1{}),
                                  get2(in, intermediate, Input2{}));
                tape[position_on_tape] += tape[IdxTypesAliveCurrent] * d2;
                tape[IdxTypesAliveCurrent] *= d1;
                evaluate_bwd(tape, in, intermediate, args<TypesAlive...>{},
                             args<Input1, LeavesAlive...>{}, args<Leaves...>{},
                             std::index_sequence<IdxTypesAlive...>{},
                             std::index_sequence<IdxTypesAliveCurrent,
                                                 IdxLeavesAlive...>{},
                             std::index_sequence<IdxAvailable...>{});
            } else if constexpr (is_input2_new) {
                constexpr auto position = idx_type2<Input1, LeavesAlive...>();
                constexpr auto position_on_tape =
                    Get<position, IdxLeavesAlive...>::value;
                double const d1 =
                    this_type::d1(get2(in, intermediate, this_type{}),
                                  get2(in, intermediate, Input1{}),
                                  get2(in, intermediate, Input2{}));
                double const d2 =
                    this_type::d2(get2(in, intermediate, this_type{}),
                                  get2(in, intermediate, Input1{}),
                                  get2(in, intermediate, Input2{}));
                tape[position_on_tape] += tape[IdxTypesAliveCurrent] * d1;
                tape[IdxTypesAliveCurrent] *= d2;
                evaluate_bwd(tape, in, intermediate,
                             args<Input2, TypesAlive...>{},
                             args<LeavesAlive...>{}, args<Leaves...>{},
                             std::index_sequence<IdxTypesAliveCurrent,
                                                 IdxTypesAlive...>{},
                             std::index_sequence<IdxLeavesAlive...>{},
                             std::index_sequence<IdxAvailable...>{});
            } else {
                constexpr auto position1 = idx_type2<Input1, LeavesAlive...>();
                constexpr auto position2 = idx_type2<Input2, TypesAlive...>();
                constexpr auto position1_on_tape =
                    Get<position1, IdxLeavesAlive...>::value;
                constexpr auto position2_on_tape =
                    Get<position2, IdxTypesAlive...>::value;
                double const d1 =
                    this_type::d1(get2(in, intermediate, this_type{}),
                                  get2(in, intermediate, Input1{}),
                                  get2(in, intermediate, Input2{}));
                double const d2 =
                    this_type::d2(get2(in, intermediate, this_type{}),
                                  get2(in, intermediate, Input1{}),
                                  get2(in, intermediate, Input2{}));
                tape[position1_on_tape] += tape[IdxTypesAliveCurrent] * d1;
                tape[position2_on_tape] += tape[IdxTypesAliveCurrent] * d2;
                evaluate_bwd(tape, in, intermediate, args<TypesAlive...>{},
                             args<LeavesAlive...>{}, args<Leaves...>{},
                             std::index_sequence<IdxTypesAlive...>{},
                             std::index_sequence<IdxLeavesAlive...>{},
                             std::index_sequence<IdxTypesAliveCurrent,
                                                 IdxAvailable...>{});
            }
        } else if constexpr (is_input2_leaf) {
            constexpr bool is_input1_new = !has_type2<Input1, TypesAlive...>();
            constexpr bool is_input2_new_leaf =
                !has_type2<Input2, LeavesAlive...>();

            if constexpr (is_input1_new && is_input2_new_leaf) {
                static_assert(sizeof...(IdxAvailable) > 0);
                evaluate_bwd_bivariate_noskip_new_loc(
                    tape, in, intermediate, args<this_type, TypesAlive...>{},
                    args<LeavesAlive...>{}, args<Leaves...>{},
                    std::index_sequence<IdxTypesAliveCurrent,
                                        IdxTypesAlive...>{},
                    std::index_sequence<IdxLeavesAlive...>{},
                    std::index_sequence<IdxAvailable...>{});
            } else if constexpr (is_input1_new) {
                constexpr auto position = idx_type2<Input2, LeavesAlive...>();
                constexpr auto position_on_tape =
                    Get<position, IdxLeavesAlive...>::value;
                double const d1 =
                    this_type::d1(get2(in, intermediate, this_type{}),
                                  get2(in, intermediate, Input1{}),
                                  get2(in, intermediate, Input2{}));
                double const d2 =
                    this_type::d2(get2(in, intermediate, this_type{}),
                                  get2(in, intermediate, Input1{}),
                                  get2(in, intermediate, Input2{}));
                tape[position_on_tape] += tape[IdxTypesAliveCurrent] * d2;
                tape[IdxTypesAliveCurrent] *= d1;
                evaluate_bwd(tape, in, intermediate,
                             args<Input1, TypesAlive...>{},
                             args<LeavesAlive...>{}, args<Leaves...>{},
                             std::index_sequence<IdxTypesAliveCurrent,
                                                 IdxTypesAlive...>{},
                             std::index_sequence<IdxLeavesAlive...>{},
                             std::index_sequence<IdxAvailable...>{});
            } else if constexpr (is_input2_new_leaf) {
                constexpr auto position = idx_type2<Input1, TypesAlive...>();
                constexpr auto position_on_tape =
                    Get<position, IdxTypesAlive...>::value;
                double const d1 =
                    this_type::d1(get2(in, intermediate, this_type{}),
                                  get2(in, intermediate, Input1{}),
                                  get2(in, intermediate, Input2{}));
                double const d2 =
                    this_type::d2(get2(in, intermediate, this_type{}),
                                  get2(in, intermediate, Input1{}),
                                  get2(in, intermediate, Input2{}));
                tape[position_on_tape] += tape[IdxTypesAliveCurrent] * d1;
                tape[IdxTypesAliveCurrent] *= d2;
                evaluate_bwd(tape, in, intermediate, args<TypesAlive...>{},
                             args<Input2, LeavesAlive...>{}, args<Leaves...>{},
                             std::index_sequence<IdxTypesAlive...>{},
                             std::index_sequence<IdxTypesAliveCurrent,
                                                 IdxLeavesAlive...>{},
                             std::index_sequence<IdxAvailable...>{});
            } else {
                constexpr auto position1 = idx_type2<Input1, TypesAlive...>();
                constexpr auto position2 = idx_type2<Input2, LeavesAlive...>();
                constexpr auto position1_on_tape =
                    Get<position1, IdxTypesAlive...>::value;
                constexpr auto position2_on_tape =
                    Get<position2, IdxLeavesAlive...>::value;
                double const d1 =
                    this_type::d1(get2(in, intermediate, this_type{}),
                                  get2(in, intermediate, Input1{}),
                                  get2(in, intermediate, Input2{}));
                double const d2 =
                    this_type::d2(get2(in, intermediate, this_type{}),
                                  get2(in, intermediate, Input1{}),
                                  get2(in, intermediate, Input2{}));
                tape[position1_on_tape] += tape[IdxTypesAliveCurrent] * d1;
                tape[position2_on_tape] += tape[IdxTypesAliveCurrent] * d2;
                evaluate_bwd(tape, in, intermediate, args<TypesAlive...>{},
                             args<LeavesAlive...>{}, args<Leaves...>{},
                             std::index_sequence<IdxTypesAlive...>{},
                             std::index_sequence<IdxLeavesAlive...>{},
                             std::index_sequence<IdxTypesAliveCurrent,
                                                 IdxAvailable...>{});
            }
        } else {
            constexpr bool is_input1_new = !has_type2<Input1, TypesAlive...>();
            constexpr bool is_input2_new = !has_type2<Input2, TypesAlive...>();

            if constexpr (is_input1_new && is_input2_new) {
                if constexpr (std::is_same_v<Input1, Input2>) {
                    double const d1 =
                        this_type::d1(get2(in, intermediate, this_type{}),
                                      get2(in, intermediate, Input1{}),
                                      get2(in, intermediate, Input2{}));
                    tape[IdxTypesAliveCurrent] *= 2 * d1;
                    evaluate_bwd(tape, in, intermediate,
                                 args<Input1, TypesAlive...>{},
                                 args<LeavesAlive...>{}, args<Leaves...>{},
                                 std::index_sequence<IdxTypesAliveCurrent,
                                                     IdxTypesAlive...>{},
                                 std::index_sequence<IdxLeavesAlive...>{},
                                 std::index_sequence<IdxAvailable...>{});
                } else {
                    static_assert(sizeof...(IdxAvailable) > 0);
                    evaluate_bwd_bivariate_noskip_new_loc(
                        tape, in, intermediate,
                        args<this_type, TypesAlive...>{},
                        args<LeavesAlive...>{}, args<Leaves...>{},
                        std::index_sequence<IdxTypesAliveCurrent,
                                            IdxTypesAlive...>{},
                        std::index_sequence<IdxLeavesAlive...>{},
                        std::index_sequence<IdxAvailable...>{});
                }

            } else if constexpr (is_input1_new) {
                constexpr auto position = idx_type2<Input2, TypesAlive...>();
                constexpr auto position_on_tape =
                    Get<position, IdxTypesAlive...>::value;
                double const d1 =
                    this_type::d1(get2(in, intermediate, this_type{}),
                                  get2(in, intermediate, Input1{}),
                                  get2(in, intermediate, Input2{}));
                double const d2 =
                    this_type::d2(get2(in, intermediate, this_type{}),
                                  get2(in, intermediate, Input1{}),
                                  get2(in, intermediate, Input2{}));
                tape[position_on_tape] += tape[IdxTypesAliveCurrent] * d2;
                tape[IdxTypesAliveCurrent] *= d1;
                evaluate_bwd(tape, in, intermediate,
                             args<Input1, TypesAlive...>{},
                             args<LeavesAlive...>{}, args<Leaves...>{},
                             std::index_sequence<IdxTypesAliveCurrent,
                                                 IdxTypesAlive...>{},
                             std::index_sequence<IdxLeavesAlive...>{},
                             std::index_sequence<IdxAvailable...>{});
            } else if constexpr (is_input2_new) {
                constexpr auto position = idx_type2<Input1, TypesAlive...>();
                constexpr auto position_on_tape =
                    Get<position, IdxTypesAlive...>::value;
                double const d1 =
                    this_type::d1(get2(in, intermediate, this_type{}),
                                  get2(in, intermediate, Input1{}),
                                  get2(in, intermediate, Input2{}));
                double const d2 =
                    this_type::d2(get2(in, intermediate, this_type{}),
                                  get2(in, intermediate, Input1{}),
                                  get2(in, intermediate, Input2{}));
                tape[position_on_tape] += tape[IdxTypesAliveCurrent] * d1;
                tape[IdxTypesAliveCurrent] *= d2;
                evaluate_bwd(tape, in, intermediate,
                             args<Input2, TypesAlive...>{},
                             args<LeavesAlive...>{}, args<Leaves...>{},
                             std::index_sequence<IdxTypesAliveCurrent,
                                                 IdxTypesAlive...>{},
                             std::index_sequence<IdxLeavesAlive...>{},
                             std::index_sequence<IdxAvailable...>{});
            } else {
                if constexpr (std::is_same_v<Input1, Input2>) {
                    constexpr auto position =
                        idx_type2<Input1, TypesAlive...>();
                    static_assert(position ==
                                  idx_type2<Input2, TypesAlive...>());
                    constexpr auto position_on_tape =
                        Get<position, IdxTypesAlive...>::value;
                    double const d1 =
                        this_type::d1(get2(in, intermediate, this_type{}),
                                      get2(in, intermediate, Input1{}),
                                      get2(in, intermediate, Input2{}));
                    tape[position_on_tape] +=
                        2 * tape[IdxTypesAliveCurrent] * d1;
                } else {
                    constexpr auto position1 =
                        idx_type2<Input1, TypesAlive...>();
                    constexpr auto position2 =
                        idx_type2<Input2, TypesAlive...>();
                    constexpr auto position1_on_tape =
                        Get<position1, IdxTypesAlive...>::value;
                    constexpr auto position2_on_tape =
                        Get<position2, IdxTypesAlive...>::value;
                    double const d1 =
                        this_type::d1(get2(in, intermediate, this_type{}),
                                      get2(in, intermediate, Input1{}),
                                      get2(in, intermediate, Input2{}));
                    double const d2 =
                        this_type::d2(get2(in, intermediate, this_type{}),
                                      get2(in, intermediate, Input1{}),
                                      get2(in, intermediate, Input2{}));
                    tape[position1_on_tape] += tape[IdxTypesAliveCurrent] * d1;
                    tape[position2_on_tape] += tape[IdxTypesAliveCurrent] * d2;
                }
                evaluate_bwd(tape, in, intermediate, args<TypesAlive...>{},
                             args<LeavesAlive...>{}, args<Leaves...>{},
                             std::index_sequence<IdxTypesAlive...>{},
                             std::index_sequence<IdxLeavesAlive...>{},
                             std::index_sequence<IdxTypesAliveCurrent,
                                                 IdxAvailable...>{});
            }
        }
    } else if constexpr (input1_has_0_deriv) {
        double const d2 = this_type::d2(get2(in, intermediate, this_type{}),
                                        get2(in, intermediate, Input1{}),
                                        get2(in, intermediate, Input2{}));
        constexpr bool is_input_leaf = has_type2<Input2, Leaves...>();
        if constexpr (is_input_leaf) {
            constexpr bool is_input_new_leaf =
                !has_type2<Input2, LeavesAlive...>();

            if constexpr (is_input_new_leaf) {
                tape[IdxTypesAliveCurrent] *= d2;
                evaluate_bwd(tape, in, intermediate, args<TypesAlive...>{},
                             args<Input2, LeavesAlive...>{}, args<Leaves...>{},
                             std::index_sequence<IdxTypesAlive...>{},
                             std::index_sequence<IdxTypesAliveCurrent,
                                                 IdxLeavesAlive...>{},
                             std::index_sequence<IdxAvailable...>{});
            } else {
                constexpr auto position = idx_type2<Input2, LeavesAlive...>();
                constexpr auto position_on_tape =
                    Get<position, IdxLeavesAlive...>::value;
                tape[position_on_tape] += tape[IdxTypesAliveCurrent] * d2;
                evaluate_bwd(tape, in, intermediate, args<TypesAlive...>{},
                             args<LeavesAlive...>{}, args<Leaves...>{},
                             std::index_sequence<IdxTypesAlive...>{},
                             std::index_sequence<IdxLeavesAlive...>{},
                             std::index_sequence<IdxTypesAliveCurrent,
                                                 IdxAvailable...>{});
            }
        } else {
            constexpr bool is_input_new = !has_type2<Input2, TypesAlive...>();

            if constexpr (is_input_new) {
                tape[IdxTypesAliveCurrent] *= d2;
                evaluate_bwd(tape, in, intermediate,
                             args<Input2, TypesAlive...>{},
                             args<LeavesAlive...>{}, args<Leaves...>{},
                             std::index_sequence<IdxTypesAliveCurrent,
                                                 IdxTypesAlive...>{},
                             std::index_sequence<IdxLeavesAlive...>{},
                             std::index_sequence<IdxAvailable...>{});
            } else {
                constexpr auto position = idx_type2<Input2, TypesAlive...>();
                constexpr auto position_on_tape =
                    Get<position, IdxTypesAlive...>::value;
                tape[position_on_tape] += tape[IdxTypesAliveCurrent] * d2;
                evaluate_bwd(tape, in, intermediate, args<TypesAlive...>{},
                             args<LeavesAlive...>{}, args<Leaves...>{},
                             std::index_sequence<IdxTypesAlive...>{},
                             std::index_sequence<IdxLeavesAlive...>{},
                             std::index_sequence<IdxTypesAliveCurrent,
                                                 IdxAvailable...>{});
            }
        }
    } else if constexpr (input2_has_0_deriv) {
        double const d1 = this_type::d1(get2(in, intermediate, this_type{}),
                                        get2(in, intermediate, Input1{}),
                                        get2(in, intermediate, Input2{}));

        constexpr bool is_input_leaf = has_type2<Input1, Leaves...>();
        if constexpr (is_input_leaf) {
            constexpr bool is_input_new_leaf =
                !has_type2<Input1, LeavesAlive...>();

            if constexpr (is_input_new_leaf) {
                tape[IdxTypesAliveCurrent] *= d1;
                evaluate_bwd(tape, in, intermediate, args<TypesAlive...>{},
                             args<Input1, LeavesAlive...>{}, args<Leaves...>{},
                             std::index_sequence<IdxTypesAlive...>{},
                             std::index_sequence<IdxTypesAliveCurrent,
                                                 IdxLeavesAlive...>{},
                             std::index_sequence<IdxAvailable...>{});
            } else {
                constexpr auto position = idx_type2<Input1, LeavesAlive...>();
                constexpr auto position_on_tape =
                    Get<position, IdxLeavesAlive...>::value;
                tape[position_on_tape] += tape[IdxTypesAliveCurrent] * d1;
                evaluate_bwd(tape, in, intermediate, args<TypesAlive...>{},
                             args<LeavesAlive...>{}, args<Leaves...>{},
                             std::index_sequence<IdxTypesAlive...>{},
                             std::index_sequence<IdxLeavesAlive...>{},
                             std::index_sequence<IdxTypesAliveCurrent,
                                                 IdxAvailable...>{});
            }
        } else {
            constexpr bool is_input_new = !has_type2<Input1, TypesAlive...>();

            if constexpr (is_input_new) {
                tape[IdxTypesAliveCurrent] *= d1;
                evaluate_bwd(tape, in, intermediate,
                             args<Input1, TypesAlive...>{},
                             args<LeavesAlive...>{}, args<Leaves...>{},
                             std::index_sequence<IdxTypesAliveCurrent,
                                                 IdxTypesAlive...>{},
                             std::index_sequence<IdxLeavesAlive...>{},
                             std::index_sequence<IdxAvailable...>{});
            } else {
                constexpr auto position = idx_type2<Input1, TypesAlive...>();
                constexpr auto position_on_tape =
                    Get<position, IdxTypesAlive...>::value;
                tape[position_on_tape] += tape[IdxTypesAliveCurrent] * d1;
                evaluate_bwd(tape, in, intermediate, args<TypesAlive...>{},
                             args<LeavesAlive...>{}, args<Leaves...>{},
                             std::index_sequence<IdxTypesAlive...>{},
                             std::index_sequence<IdxLeavesAlive...>{},
                             std::index_sequence<IdxTypesAliveCurrent,
                                                 IdxAvailable...>{});
            }
        }
    }
}

template <class... InputTypes, class... IntermediateTypes, std::size_t N,
          template <class> class Univariate, class Input,
          typename... TypesAlive, typename... LeavesAlive, typename... Leaves,
          std::size_t IdxTypesAliveCurrent, std::size_t... IdxTypesAlive,
          std::size_t... IdxLeavesAlive, std::size_t... IdxAvailable>
inline void evaluate_bwd_univariate_noskip(
    std::array<double, N> &tape, TapeIntermediate<InputTypes...> const &in,
    TapeIntermediate<IntermediateTypes...> const &intermediate,
    args<Univariate<Input> const, TypesAlive...> const &,
    args<LeavesAlive...> const &, args<Leaves...> const &,
    std::index_sequence<IdxTypesAliveCurrent, IdxTypesAlive...> const &,
    std::index_sequence<IdxLeavesAlive...> const &,
    std::index_sequence<IdxAvailable...> const &) {

    using this_type = Univariate<Input> const;

    double const d = this_type::d(get2(in, intermediate, this_type{}),
                                  get2(in, intermediate, Input{}));

    constexpr bool is_input_leaf = has_type2<Input, Leaves...>();
    if constexpr (is_input_leaf) {
        constexpr bool is_input_new_leaf = !has_type2<Input, LeavesAlive...>();

        if constexpr (is_input_new_leaf) {
            tape[IdxTypesAliveCurrent] *= d;
            evaluate_bwd(
                tape, in, intermediate, args<TypesAlive...>{},
                args<Input, LeavesAlive...>{}, args<Leaves...>{},
                std::index_sequence<IdxTypesAlive...>{},
                std::index_sequence<IdxTypesAliveCurrent, IdxLeavesAlive...>{},
                std::index_sequence<IdxAvailable...>{});
        } else {
            constexpr auto position = idx_type2<Input, LeavesAlive...>();
            constexpr auto position_on_tape =
                Get<position, IdxLeavesAlive...>::value;
            tape[position_on_tape] += tape[IdxTypesAliveCurrent] * d;
            evaluate_bwd(
                tape, in, intermediate, args<TypesAlive...>{},
                args<LeavesAlive...>{}, args<Leaves...>{},
                std::index_sequence<IdxTypesAlive...>{},
                std::index_sequence<IdxLeavesAlive...>{},
                std::index_sequence<IdxTypesAliveCurrent, IdxAvailable...>{});
        }
    } else {
        constexpr bool is_input_new = !has_type2<Input, TypesAlive...>();

        if constexpr (is_input_new) {
            tape[IdxTypesAliveCurrent] *= d;
            evaluate_bwd(
                tape, in, intermediate, args<Input, TypesAlive...>{},
                args<LeavesAlive...>{}, args<Leaves...>{},
                std::index_sequence<IdxTypesAliveCurrent, IdxTypesAlive...>{},
                std::index_sequence<IdxLeavesAlive...>{},
                std::index_sequence<IdxAvailable...>{});
        } else {
            constexpr auto position = idx_type2<Input, TypesAlive...>();
            constexpr auto position_on_tape =
                Get<position, IdxTypesAlive...>::value;
            tape[position_on_tape] += tape[IdxTypesAliveCurrent] * d;
            evaluate_bwd(
                tape, in, intermediate, args<TypesAlive...>{},
                args<LeavesAlive...>{}, args<Leaves...>{},
                std::index_sequence<IdxTypesAlive...>{},
                std::index_sequence<IdxLeavesAlive...>{},
                std::index_sequence<IdxTypesAliveCurrent, IdxAvailable...>{});
        }
    }
}

template <class... InputTypes, class... IntermediateTypes, std::size_t N,
          class FirstType, class SecondType, typename... TypesAlive,
          typename... LeavesAlive, typename... Leaves,
          std::size_t IdxTypesAliveFirst, std::size_t IdxTypesAliveSecond,
          std::size_t... IdxTypesAlive, std::size_t... IdxLeavesAlive,
          std::size_t... IdxAvailable>
inline void
evaluate_skip(std::array<double, N> &tape,
              TapeIntermediate<InputTypes...> const &in,
              TapeIntermediate<IntermediateTypes...> const &intermediate,
              args<FirstType, SecondType, TypesAlive...> const &,
              args<LeavesAlive...> const &, args<Leaves...> const &,
              std::index_sequence<IdxTypesAliveFirst, IdxTypesAliveSecond,
                                  IdxTypesAlive...> const &,
              std::index_sequence<IdxLeavesAlive...> const &,
              std::index_sequence<IdxAvailable...> const &) {

    evaluate_bwd(tape, in, intermediate,
                 args<SecondType, TypesAlive...,
                      FirstType>{} /* we place first at the end */,
                 args<LeavesAlive...>{}, args<Leaves...>{},
                 std::index_sequence<
                     IdxTypesAliveSecond, IdxTypesAlive...,
                     IdxTypesAliveFirst>{} /* we place first at the end */,
                 std::index_sequence<IdxLeavesAlive...>{},
                 std::index_sequence<IdxAvailable...>{});
}

template <class... InputTypes, class... IntermediateTypes, std::size_t N,
          template <class> class Univariate, class Input,
          typename... TypesAlive, typename... LeavesAlive, typename... Leaves,
          std::size_t... IdxTypesAlive, std::size_t... IdxLeavesAlive,
          std::size_t... IdxAvailable>
inline void
evaluate_bwd(std::array<double, N> &tape,
             TapeIntermediate<InputTypes...> const &in,
             TapeIntermediate<IntermediateTypes...> const &intermediate,
             args<Univariate<Input> const, TypesAlive...> const &,
             args<LeavesAlive...> const &, args<Leaves...> const &,
             std::index_sequence<IdxTypesAlive...> const &,
             std::index_sequence<IdxLeavesAlive...> const &,
             std::index_sequence<IdxAvailable...> const &) {

    using this_type = Univariate<Input> const;
    static_assert(!has_type2<this_type, Leaves...>());

    constexpr bool other_types_depend_on_this =
        (depends<TypesAlive, this_type>::call() || ...);

    if constexpr (other_types_depend_on_this) {
        evaluate_skip(tape, in, intermediate, args<this_type, TypesAlive...>{},
                      args<LeavesAlive...>{}, args<Leaves...>{},
                      std::index_sequence<IdxTypesAlive...>{},
                      std::index_sequence<IdxLeavesAlive...>{},
                      std::index_sequence<IdxAvailable...>{});
    } else {
        evaluate_bwd_univariate_noskip(
            tape, in, intermediate, args<this_type, TypesAlive...>{},
            args<LeavesAlive...>{}, args<Leaves...>{},
            std::index_sequence<IdxTypesAlive...>{},
            std::index_sequence<IdxLeavesAlive...>{},
            std::index_sequence<IdxAvailable...>{});
    }
}

template <class... InputTypes, class... IntermediateTypes, std::size_t N,
          template <class, class> class Bivariate, class Input1, class Input2,
          typename... TypesAlive, typename... LeavesAlive, typename... Leaves,
          std::size_t... IdxTypesAlive, std::size_t... IdxLeavesAlive,
          std::size_t... IdxAvailable>
inline void
evaluate_bwd(std::array<double, N> &tape,
             TapeIntermediate<InputTypes...> const &in,
             TapeIntermediate<IntermediateTypes...> const &intermediate,
             args<Bivariate<Input1, Input2> const, TypesAlive...> const &,
             args<LeavesAlive...> const &, args<Leaves...> const &,
             std::index_sequence<IdxTypesAlive...> const &,
             std::index_sequence<IdxLeavesAlive...> const &,
             std::index_sequence<IdxAvailable...> const &) {

    using this_type = Bivariate<Input1, Input2> const;
    static_assert(!has_type2<this_type, Leaves...>());

    constexpr bool other_types_depend_on_this =
        (depends<TypesAlive, this_type>::call() || ...);

    if constexpr (other_types_depend_on_this) {
        evaluate_skip(tape, in, intermediate, args<this_type, TypesAlive...>{},
                      args<LeavesAlive...>{}, args<Leaves...>{},
                      std::index_sequence<IdxTypesAlive...>{},
                      std::index_sequence<IdxLeavesAlive...>{},
                      std::index_sequence<IdxAvailable...>{});
    } else {
        evaluate_bwd_bivariate_noskip(tape, in, intermediate,
                                      args<this_type, TypesAlive...>{},
                                      args<LeavesAlive...>{}, args<Leaves...>{},
                                      std::index_sequence<IdxTypesAlive...>{},
                                      std::index_sequence<IdxLeavesAlive...>{},
                                      std::index_sequence<IdxAvailable...>{});
    }
}

template <class... InputTypes, class... IntermediateTypes, std::size_t N,
          class this_type, typename... Leaves, std::size_t... IdxTypesAlive,
          std::size_t... IdxLeavesAlive, std::size_t FirstIdxAvailable,
          std::size_t... IdxAvailable>
inline void evaluate_bwd_first(
    std::array<double, N> &tape, TapeIntermediate<InputTypes...> const &in,
    TapeIntermediate<IntermediateTypes...> const &intermediate, double init,
    args<this_type const> const &, args<Leaves...> const &,
    std::index_sequence<FirstIdxAvailable, IdxAvailable...> const &) {
    tape[FirstIdxAvailable] = init;
    constexpr bool is_input_leaf = has_type2<this_type, Leaves...>();
    if constexpr (!is_input_leaf) {
        evaluate_bwd(
            tape, in, intermediate, args<this_type const>{}, args<>{},
            args<Leaves...>{}, std::index_sequence<FirstIdxAvailable>{},
            std::index_sequence<>{}, std::index_sequence<IdxAvailable...>{});
    }
}

} // namespace detail

template <class... InputTypes, class... IntermediateTypes, class Output,
          class... Leaves>
inline auto
evaluate_bwd(detail::TapeIntermediate<InputTypes...> const &in,
             detail::TapeIntermediate<IntermediateTypes...> const &intermediate,
             Output const & /* out */, double init,
             Leaves const &.../* leaves */)
    -> std::array<double, sizeof...(Leaves)> {

    static_assert((depends<Output, Leaves const>::call() && ...));
    std::array<double, sizeof...(Leaves)> results = {};
    if constexpr (results.empty()) {
        return results;
    } else {
        constexpr std::size_t size = detail::tape_size(
            args<Output const>{}, args<>{}, args<Leaves const...>{});
        std::array<double, size> tape{};
        detail::evaluate_bwd_first(
            tape, in, intermediate, init, args<Output const>{},
            args<Leaves const...>{}, std::make_index_sequence<size>{});
        for (std::size_t index = 0; index < sizeof...(Leaves); ++index) {
            results[index] = tape[index];
        }
        return results;
    }
}

} // namespace adhoc

#endif // NEW_EVALUATE_BWD_HPP