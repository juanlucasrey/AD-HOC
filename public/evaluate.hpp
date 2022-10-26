#ifndef EVALUATE_HPP
#define EVALUATE_HPP

#include <adhoc2.hpp>
#include <tape_size.hpp>

namespace adhoc2 {
namespace detail {

// template <std::size_t N, typename... TypesAlive, typename... LeavesAlive,
//           typename... Leaves, std::size_t... IdxTypesAlive,
//           std::size_t... IdxLeavesAlive, std::size_t... IdxAvailable>
// static inline void evaluate(std::array<double, N> &, TypesAlive const &...,
//                             args<LeavesAlive...> const &,
//                             args<Leaves...> const &,
//                             std::index_sequence<IdxTypesAlive...> const &,
//                             std::index_sequence<IdxLeavesAlive...> const &,
//                             std::index_sequence<IdxAvailable...> const &);

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
static inline void
evaluate(std::array<double, N> &tape, args<LeavesAlive...> const &,
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
// template <std::size_t N, typename ThisType, typename... TypesAlive,
//           typename... LeavesAlive, typename... Leaves>
// static void evaluate(std::array<double, N> & /* tape */,
//                      ThisType const & /* in */,
//                      const args<ThisType, TypesAlive...> &,
//                      const args<LeavesAlive...> &, const args<Leaves...> &);

// template <std::size_t N, int ID, typename... TypesAlive,
//           typename... LeavesAlive, typename... Leaves>
// static void evaluate(std::array<double, N> & /* tape */,
//                      adouble_aux<ID> const & /* in */,
//                      const args<adouble_aux<ID> const, TypesAlive...> &,
//                      const args<LeavesAlive...> &, const args<Leaves...> &) {

//     //     using this_type = adouble_aux<N>;

//     //     // because we are on a leave, this_type should be included
//     //     static_assert(has_type2<this_type, Leaves...>());

//     //     // this_type should only be once in TypesAlive
//     //     static_assert(!has_type2<this_type, TypesAlive...>());

//     //     constexpr std::size_t curent_tape_size =
//     //         1 + sizeof...(TypesAlive) + sizeof...(LeavesAlive);
//     //     return curent_tape_size;
// }

// template <typename NextType, typename... TypesAlive, typename... LeavesAlive,
//           typename... Leaves>
// constexpr static auto tape_size_aux(const args<NextType, TypesAlive...> &,
//                                     const args<LeavesAlive...> &,
//                                     const args<Leaves...> &) -> std::size_t {
//     return tape_size(args<NextType, TypesAlive...>{}, args<LeavesAlive...>{},
//                      args<Leaves...>{});
// }

// template <typename... LeavesAlive, typename... Leaves>
// constexpr static auto tape_size_aux(const args<> &,
//                                     const args<LeavesAlive...> &,
//                                     const args<Leaves...> &) -> std::size_t {
//     return sizeof...(LeavesAlive);
// }

// template <class Input, class this_type, typename... TypesAlive,
//           typename... LeavesAlive, typename... Leaves>
// constexpr static auto
// tape_size_next_univariate_noskip(args<this_type, TypesAlive...> const &,
//                                  args<LeavesAlive...> const &,
//                                  args<Leaves...> const &) -> std::size_t {
//     constexpr bool is_input_leaf = has_type2<Input, Leaves...>();
//     if constexpr (is_input_leaf) {
//         constexpr bool is_input_new_leaf = !has_type2<Input,
//         LeavesAlive...>();

//         if constexpr (is_input_new_leaf) {
//             return tape_size_aux(args<TypesAlive...>{},
//                                  args<Input, LeavesAlive...>{},
//                                  args<Leaves...>{});
//         } else {
//             return tape_size_aux(args<TypesAlive...>{},
//             args<LeavesAlive...>{},
//                                  args<Leaves...>{});
//         }
//     } else {
//         constexpr bool is_input_new = !has_type2<Input, TypesAlive...>();

//         if constexpr (is_input_new) {
//             return tape_size_aux(args<Input, TypesAlive...>{},
//                                  args<LeavesAlive...>{}, args<Leaves...>{});
//         } else {
//             return tape_size_aux(args<TypesAlive...>{},
//             args<LeavesAlive...>{},
//                                  args<Leaves...>{});
//         }
//     }
// }

// template <class Input, class this_type, typename... TypesAlive,
//           typename... LeavesAlive, typename... Leaves>
// constexpr static auto
// tape_size_next_univariate(args<this_type, TypesAlive...> const &,
//                           args<LeavesAlive...> const &, args<Leaves...> const
//                           &)
//     -> std::size_t {
//     constexpr bool other_types_depend_on_this =
//         (TypesAlive::template depends_in<this_type>() || ...);

//     if constexpr (other_types_depend_on_this) {
//         return skip_type(args<this_type const, TypesAlive...>{},
//                          args<LeavesAlive...>{}, args<Leaves...>{});
//     } else {
//         return tape_size_next_univariate_noskip<Input>(
//             args<this_type, TypesAlive...>{}, args<LeavesAlive...>{},
//             args<Leaves...>{});
//     }
// }

// template <class Input1, class Input2, class this_type, typename...
// TypesAlive,
//           typename... LeavesAlive, typename... Leaves>
// constexpr static auto
// tape_size_next_bivariate(args<this_type, TypesAlive...> const &,
//                          args<LeavesAlive...> const &, args<Leaves...> const
//                          &)
//     -> std::size_t {

//     constexpr bool other_types_depend_on_this =
//         (TypesAlive::template depends_in<this_type>() || ...);

//     if constexpr (other_types_depend_on_this) {
//         return skip_type(args<this_type const, TypesAlive...>{},
//                          args<LeavesAlive...>{}, args<Leaves...>{});
//     } else {
//         constexpr bool input1_has_0_deriv =
//             !static_cast<bool>(Input1::template depends3<Leaves...>());
//         constexpr bool input2_has_0_deriv =
//             !static_cast<bool>(Input2::template depends3<Leaves...>());
//         static_assert(!input1_has_0_deriv || !input2_has_0_deriv);

//         constexpr bool is_input1_leaf = has_type2<Input1, Leaves...>();
//         constexpr bool is_input2_leaf = has_type2<Input2, Leaves...>();
//         if constexpr (!input1_has_0_deriv && !input2_has_0_deriv) {
//             if constexpr (is_input1_leaf && is_input2_leaf) {
//                 constexpr bool is_input1_new_leaf =
//                     !has_type2<Input1, LeavesAlive...>();
//                 constexpr bool is_input2_new_leaf =
//                     !has_type2<Input2, LeavesAlive...>();

//                 if constexpr (is_input1_new_leaf && is_input2_new_leaf) {
//                     return tape_size_aux(args<TypesAlive...>{},
//                                          args<Input1, Input2,
//                                          LeavesAlive...>{},
//                                          args<Leaves...>{});
//                 } else if constexpr (is_input1_new_leaf) {
//                     return tape_size_aux(args<TypesAlive...>{},
//                                          args<Input1, LeavesAlive...>{},
//                                          args<Leaves...>{});
//                 } else if constexpr (is_input2_new_leaf) {
//                     return tape_size_aux(args<TypesAlive...>{},
//                                          args<Input2, LeavesAlive...>{},
//                                          args<Leaves...>{});
//                 } else {
//                     return tape_size_aux(args<TypesAlive...>{},
//                                          args<LeavesAlive...>{},
//                                          args<Leaves...>{});
//                 }
//             } else if constexpr (is_input1_leaf) {
//                 constexpr bool is_input1_new_leaf =
//                     !has_type2<Input1, LeavesAlive...>();
//                 constexpr bool is_input2_new =
//                     !has_type2<Input2, TypesAlive...>();

//                 if constexpr (is_input1_new_leaf && is_input2_new) {
//                     return tape_size_aux(args<Input2, TypesAlive...>{},
//                                          args<Input1, LeavesAlive...>{},
//                                          args<Leaves...>{});
//                 } else if constexpr (is_input1_new_leaf) {
//                     return tape_size_aux(args<TypesAlive...>{},
//                                          args<Input1, LeavesAlive...>{},
//                                          args<Leaves...>{});
//                 } else if constexpr (is_input2_new) {
//                     return tape_size_aux(args<Input2, TypesAlive...>{},
//                                          args<LeavesAlive...>{},
//                                          args<Leaves...>{});
//                 } else {
//                     return tape_size_aux(args<TypesAlive...>{},
//                                          args<LeavesAlive...>{},
//                                          args<Leaves...>{});
//                 }
//             } else if constexpr (is_input2_leaf) {
//                 constexpr bool is_input1_new =
//                     !has_type2<Input1, TypesAlive...>();
//                 constexpr bool is_input2_new_leaf =
//                     !has_type2<Input2, LeavesAlive...>();

//                 if constexpr (is_input1_new && is_input2_new_leaf) {
//                     return tape_size_aux(args<Input1, TypesAlive...>{},
//                                          args<Input2, LeavesAlive...>{},
//                                          args<Leaves...>{});
//                 } else if constexpr (is_input1_new) {
//                     return tape_size_aux(args<Input1, TypesAlive...>{},
//                                          args<LeavesAlive...>{},
//                                          args<Leaves...>{});
//                 } else if constexpr (is_input2_new_leaf) {
//                     return tape_size_aux(args<TypesAlive...>{},
//                                          args<Input2, LeavesAlive...>{},
//                                          args<Leaves...>{});
//                 } else {
//                     return tape_size_aux(args<TypesAlive...>{},
//                                          args<LeavesAlive...>{},
//                                          args<Leaves...>{});
//                 }
//             } else {
//                 constexpr bool is_input1_new =
//                     !has_type2<Input1, TypesAlive...>();
//                 constexpr bool is_input2_new =
//                     !has_type2<Input2, TypesAlive...>();

//                 if constexpr (is_input1_new && is_input2_new) {
//                     return tape_size_aux(args<Input1, Input2,
//                     TypesAlive...>{},
//                                          args<LeavesAlive...>{},
//                                          args<Leaves...>{});
//                 } else if constexpr (is_input1_new) {
//                     return tape_size_aux(args<Input1, TypesAlive...>{},
//                                          args<LeavesAlive...>{},
//                                          args<Leaves...>{});
//                 } else if constexpr (is_input2_new) {
//                     return tape_size_aux(args<Input2, TypesAlive...>{},
//                                          args<LeavesAlive...>{},
//                                          args<Leaves...>{});
//                 } else {
//                     return tape_size_aux(args<TypesAlive...>{},
//                                          args<LeavesAlive...>{},
//                                          args<Leaves...>{});
//                 }
//             }
//         } else if constexpr (input1_has_0_deriv) {
//             return tape_size_next_univariate_noskip<Input2>(
//                 args<this_type, TypesAlive...>{}, args<LeavesAlive...>{},
//                 args<Leaves...>{});
//         } else if constexpr (input2_has_0_deriv) {
//             return tape_size_next_univariate_noskip<Input1>(
//                 args<this_type, TypesAlive...>{}, args<LeavesAlive...>{},
//                 args<Leaves...>{});
//         }
//     }
// }

template <std::size_t N, class FirstType, class SecondType,
          typename... TypesAlive, typename... LeavesAlive, typename... Leaves,
          std::size_t IdxTypesAliveFirst, std::size_t IdxTypesAliveSecond,
          std::size_t... IdxTypesAlive, std::size_t... IdxLeavesAlive,
          std::size_t... IdxAvailable>
static inline void
evaluate_skip(std::array<double, N> &tape, FirstType const &first,
              SecondType const &second, TypesAlive const &...next,
              args<LeavesAlive...> const &, args<Leaves...> const &,
              std::index_sequence<IdxTypesAliveFirst, IdxTypesAliveSecond,
                                  IdxTypesAlive...> const &,
              std::index_sequence<IdxLeavesAlive...> const &,
              std::index_sequence<IdxAvailable...> const &) {

    evaluate(tape, second, next..., first /* we place first at the end */,
             args<LeavesAlive...>{}, args<Leaves...>{},
             std::index_sequence<
                 IdxTypesAliveSecond, IdxTypesAlive...,
                 IdxTypesAliveFirst>{} /* we place first at the end */,
             std::index_sequence<IdxLeavesAlive...>{},
             std::index_sequence<IdxAvailable...>{});
}

template <class Input1, class Input2, std::size_t N, class this_type,
          typename... TypesAlive, typename... LeavesAlive, typename... Leaves,
          std::size_t... IdxTypesAlive, std::size_t... IdxLeavesAlive,
          std::size_t... IdxAvailable>
static inline void
evaluate_bivariate(std::array<double, N> &tape, this_type const &in,
                   TypesAlive const &...next, args<LeavesAlive...> const &,
                   args<Leaves...> const &,
                   std::index_sequence<IdxTypesAlive...> const &,
                   std::index_sequence<IdxLeavesAlive...> const &,
                   std::index_sequence<IdxAvailable...> const &) {

    static_assert(!has_type2<this_type, Leaves...>());

    constexpr bool other_types_depend_on_this =
        (TypesAlive::template depends_in<this_type>() || ...);

    if constexpr (other_types_depend_on_this) {
        evaluate_skip(tape, in, next..., args<LeavesAlive...>{},
                      args<Leaves...>{},
                      std::index_sequence<IdxTypesAlive...>{},
                      std::index_sequence<IdxLeavesAlive...>{},
                      std::index_sequence<IdxAvailable...>{});
    } else {
        // evaluate_univariate_noskip<Input>(
        //     tape, in, next..., args<LeavesAlive...>{}, args<Leaves...>{},
        //     std::index_sequence<IdxTypesAlive...>{},
        //     std::index_sequence<IdxLeavesAlive...>{},
        //     std::index_sequence<IdxAvailable...>{});
    }
}

template <class Input, std::size_t N, class this_type, typename... TypesAlive,
          typename... LeavesAlive, typename... Leaves,
          std::size_t IdxTypesAliveCurrent, std::size_t... IdxTypesAlive,
          std::size_t... IdxLeavesAlive, std::size_t... IdxAvailable>
static inline void evaluate_univariate_noskip(
    std::array<double, N> &tape, this_type const &in, TypesAlive const &...next,
    args<LeavesAlive...> const &, args<Leaves...> const &,
    std::index_sequence<IdxTypesAliveCurrent, IdxTypesAlive...> const &,
    std::index_sequence<IdxLeavesAlive...> const &,
    std::index_sequence<IdxAvailable...> const &) {

    constexpr bool is_input_leaf = has_type2<Input, Leaves...>();
    if constexpr (is_input_leaf) {
        constexpr bool is_input_new_leaf = !has_type2<Input, LeavesAlive...>();

        if constexpr (is_input_new_leaf) {
            tape[IdxTypesAliveCurrent] *= in.d();
            evaluate(
                tape, next..., args<Input, LeavesAlive...>{}, args<Leaves...>{},
                std::index_sequence<IdxTypesAlive...>{},
                std::index_sequence<IdxTypesAliveCurrent, IdxLeavesAlive...>{},
                std::index_sequence<IdxAvailable...>{});
        } else {
            constexpr auto position = idx_type2<Input, LeavesAlive...>();
            tape[position] += tape[IdxTypesAliveCurrent] * in.d();
            evaluate(
                tape, next..., args<LeavesAlive...>{}, args<Leaves...>{},
                std::index_sequence<IdxTypesAlive...>{},
                std::index_sequence<IdxLeavesAlive...>{},
                std::index_sequence<IdxTypesAliveCurrent, IdxAvailable...>{});
        }
    } else {
        constexpr bool is_input_new = !has_type2<Input, TypesAlive...>();

        if constexpr (is_input_new) {
            tape[IdxTypesAliveCurrent] *= in.d();
            evaluate(
                tape, in.input(), next..., args<LeavesAlive...>{},
                args<Leaves...>{},
                std::index_sequence<IdxTypesAliveCurrent, IdxTypesAlive...>{},
                std::index_sequence<IdxLeavesAlive...>{},
                std::index_sequence<IdxAvailable...>{});
        } else {
            constexpr auto position = idx_type2<Input, TypesAlive...>();
            tape[position] += tape[IdxTypesAliveCurrent] * in.d();
            evaluate(
                tape, next..., args<LeavesAlive...>{}, args<Leaves...>{},
                std::index_sequence<IdxTypesAlive...>{},
                std::index_sequence<IdxLeavesAlive...>{},
                std::index_sequence<IdxTypesAliveCurrent, IdxAvailable...>{});
        }
    }
}

template <class Input, std::size_t N, class this_type, typename... TypesAlive,
          typename... LeavesAlive, typename... Leaves,
          std::size_t... IdxTypesAlive, std::size_t... IdxLeavesAlive,
          std::size_t... IdxAvailable>
static inline void
evaluate_univariate(std::array<double, N> &tape, this_type const &in,
                    TypesAlive const &...next, args<LeavesAlive...> const &,
                    args<Leaves...> const &,
                    std::index_sequence<IdxTypesAlive...> const &,
                    std::index_sequence<IdxLeavesAlive...> const &,
                    std::index_sequence<IdxAvailable...> const &) {

    static_assert(!has_type2<this_type, Leaves...>());

    constexpr bool other_types_depend_on_this =
        (TypesAlive::template depends_in<this_type>() || ...);

    if constexpr (other_types_depend_on_this) {
        evaluate_skip(tape, in, next..., args<LeavesAlive...>{},
                      args<Leaves...>{},
                      std::index_sequence<IdxTypesAlive...>{},
                      std::index_sequence<IdxLeavesAlive...>{},
                      std::index_sequence<IdxAvailable...>{});
    } else {
        evaluate_univariate_noskip<Input>(
            tape, in, next..., args<LeavesAlive...>{}, args<Leaves...>{},
            std::index_sequence<IdxTypesAlive...>{},
            std::index_sequence<IdxLeavesAlive...>{},
            std::index_sequence<IdxAvailable...>{});
    }
}

template <std::size_t N, class Input, class this_type, typename... TypesAlive,
          typename... LeavesAlive, typename... Leaves,
          std::size_t... IdxTypesAlive, std::size_t... IdxLeavesAlive,
          std::size_t... IdxAvailable>
static inline void
evaluate(std::array<double, N> &tape, Univariate<Input, this_type> const &in,
         TypesAlive const &...next, args<LeavesAlive...> const &,
         args<Leaves...> const &, std::index_sequence<IdxTypesAlive...> const &,
         std::index_sequence<IdxLeavesAlive...> const &,
         std::index_sequence<IdxAvailable...> const &) {
    static_assert((sizeof...(TypesAlive) + 1) == sizeof...(IdxTypesAlive));
    static_assert(sizeof...(LeavesAlive) == sizeof...(IdxLeavesAlive));
    static_assert(sizeof...(IdxTypesAlive) + sizeof...(IdxLeavesAlive) +
                      sizeof...(IdxAvailable) ==
                  N);
    evaluate_univariate<Input>(tape, *static_cast<this_type const *>(&in),
                               next..., args<LeavesAlive...>{},
                               args<Leaves...>{},
                               std::index_sequence<IdxTypesAlive...>{},
                               std::index_sequence<IdxLeavesAlive...>{},
                               std::index_sequence<IdxAvailable...>{});
}

template <std::size_t N, class Input1, class Input2, class this_type,
          typename... TypesAlive, typename... LeavesAlive, typename... Leaves,
          std::size_t... IdxTypesAlive, std::size_t... IdxLeavesAlive,
          std::size_t... IdxAvailable>
static inline void
evaluate(std::array<double, N> &tape,
         Bivariate<Input1, Input2, this_type> const &in,
         TypesAlive const &...next, args<LeavesAlive...> const &,
         args<Leaves...> const &, std::index_sequence<IdxTypesAlive...> const &,
         std::index_sequence<IdxLeavesAlive...> const &,
         std::index_sequence<IdxAvailable...> const &) {
    static_assert((sizeof...(TypesAlive) + 1) == sizeof...(IdxTypesAlive));
    static_assert(sizeof...(LeavesAlive) == sizeof...(IdxLeavesAlive));
    static_assert(sizeof...(IdxTypesAlive) + sizeof...(IdxLeavesAlive) +
                      sizeof...(IdxAvailable) ==
                  N);
    evaluate_bivariate<Input1, Input2>(
        tape, *static_cast<this_type const *>(&in), next...,
        args<LeavesAlive...>{}, args<Leaves...>{},
        std::index_sequence<IdxTypesAlive...>{},
        std::index_sequence<IdxLeavesAlive...>{},
        std::index_sequence<IdxAvailable...>{});
}

template <std::size_t N, class this_type, typename... Leaves,
          std::size_t... IdxTypesAlive, std::size_t... IdxLeavesAlive,
          std::size_t FirstIdxAvailable, std::size_t... IdxAvailable>
static inline void evaluate_first(
    std::array<double, N> &tape, this_type const &in, args<Leaves...> const &,
    std::index_sequence<FirstIdxAvailable, IdxAvailable...> const &) {
    tape[FirstIdxAvailable] = 1.0;
    evaluate(tape, in, args<>{}, args<Leaves...>{},
             std::index_sequence<FirstIdxAvailable>{}, std::index_sequence<>{},
             std::index_sequence<IdxAvailable...>{});
}

} // namespace detail

template <class Output, typename... Leaves>
constexpr static inline auto evaluate(Output const &in,
                                      Leaves const &.../* leaves */)
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

} // namespace adhoc2

#endif // EVALUATE_HPP