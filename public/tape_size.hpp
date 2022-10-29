#ifndef TAPE_SIZE_HPP
#define TAPE_SIZE_HPP

#include <adhoc2.hpp>

namespace adhoc2 {
namespace detail {

template <typename ThisType, typename... TypesAlive, typename... LeavesAlive,
          typename... Leaves>
constexpr static auto tape_size(const args<ThisType, TypesAlive...> &,
                                const args<LeavesAlive...> &,
                                const args<Leaves...> &) -> std::size_t;

template <int N, typename... TypesAlive, typename... LeavesAlive,
          typename... Leaves>
constexpr static auto
tape_size(const args<adouble_aux<N> const, TypesAlive...> &,
          const args<LeavesAlive...> &, const args<Leaves...> &)
    -> std::size_t {

    using this_type = adouble_aux<N> const;

    // this_type should only be once in TypesAlive
    static_assert(!has_type2<this_type, TypesAlive...>());

    constexpr std::size_t curent_tape_size = has_type2<this_type, Leaves...>() +
                                             sizeof...(TypesAlive) +
                                             sizeof...(LeavesAlive);
    return curent_tape_size;
}

template <typename NextType, typename... TypesAlive, typename... LeavesAlive,
          typename... Leaves>
constexpr static auto tape_size_aux(const args<NextType, TypesAlive...> &,
                                    const args<LeavesAlive...> &,
                                    const args<Leaves...> &) -> std::size_t {
    return tape_size(args<NextType, TypesAlive...>{}, args<LeavesAlive...>{},
                     args<Leaves...>{});
}

template <typename... LeavesAlive, typename... Leaves>
constexpr static auto tape_size_aux(const args<> &,
                                    const args<LeavesAlive...> &,
                                    const args<Leaves...> &) -> std::size_t {
    return sizeof...(LeavesAlive);
}

template <typename CurrentType, typename NextType, typename... TypesAlive,
          typename... LeavesAlive, typename... Leaves>
constexpr static auto
skip_type(const args<CurrentType, NextType, TypesAlive...> &,
          const args<LeavesAlive...> &, const args<Leaves...> &)
    -> std::size_t {

    return tape_size(args<NextType, TypesAlive..., CurrentType>{},
                     args<LeavesAlive...>{}, args<Leaves...>{});
}

template <class Input, class this_type, typename... TypesAlive,
          typename... LeavesAlive, typename... Leaves>
constexpr static auto
tape_size_next_univariate_noskip(args<this_type, TypesAlive...> const &,
                                 args<LeavesAlive...> const &,
                                 args<Leaves...> const &) -> std::size_t {
    constexpr bool is_input_leaf = has_type2<Input, Leaves...>();
    if constexpr (is_input_leaf) {
        constexpr bool is_input_new_leaf = !has_type2<Input, LeavesAlive...>();

        if constexpr (is_input_new_leaf) {
            return tape_size_aux(args<TypesAlive...>{},
                                 args<Input, LeavesAlive...>{},
                                 args<Leaves...>{});
        } else {
            return tape_size_aux(args<TypesAlive...>{}, args<LeavesAlive...>{},
                                 args<Leaves...>{});
        }
    } else {
        constexpr bool is_input_new = !has_type2<Input, TypesAlive...>();

        if constexpr (is_input_new) {
            return tape_size_aux(args<Input, TypesAlive...>{},
                                 args<LeavesAlive...>{}, args<Leaves...>{});
        } else {
            return tape_size_aux(args<TypesAlive...>{}, args<LeavesAlive...>{},
                                 args<Leaves...>{});
        }
    }
}

template <class Input, class this_type, typename... TypesAlive,
          typename... LeavesAlive, typename... Leaves>
constexpr static auto
tape_size_next_univariate(args<this_type, TypesAlive...> const &,
                          args<LeavesAlive...> const &, args<Leaves...> const &)
    -> std::size_t {
    constexpr bool other_types_depend_on_this =
        (TypesAlive::template depends_in<this_type>() || ...);

    if constexpr (other_types_depend_on_this) {
        return skip_type(args<this_type const, TypesAlive...>{},
                         args<LeavesAlive...>{}, args<Leaves...>{});
    } else {
        return tape_size_next_univariate_noskip<Input>(
            args<this_type, TypesAlive...>{}, args<LeavesAlive...>{},
            args<Leaves...>{});
    }
}

template <class Input1, class Input2, class this_type, typename... TypesAlive,
          typename... LeavesAlive, typename... Leaves>
constexpr static auto
tape_size_next_bivariate(args<this_type, TypesAlive...> const &,
                         args<LeavesAlive...> const &, args<Leaves...> const &)
    -> std::size_t {

    constexpr bool other_types_depend_on_this =
        (TypesAlive::template depends_in<this_type>() || ...);

    if constexpr (other_types_depend_on_this) {
        return skip_type(args<this_type const, TypesAlive...>{},
                         args<LeavesAlive...>{}, args<Leaves...>{});
    } else {
        constexpr bool input1_has_0_deriv =
            !static_cast<bool>(Input1::template depends3<Leaves...>());
        constexpr bool input2_has_0_deriv =
            !static_cast<bool>(Input2::template depends3<Leaves...>());
        static_assert(!input1_has_0_deriv || !input2_has_0_deriv);

        constexpr bool is_input1_leaf = has_type2<Input1, Leaves...>();
        constexpr bool is_input2_leaf = has_type2<Input2, Leaves...>();
        if constexpr (!input1_has_0_deriv && !input2_has_0_deriv) {
            if constexpr (is_input1_leaf && is_input2_leaf) {
                constexpr bool is_input1_new_leaf =
                    !has_type2<Input1, LeavesAlive...>();
                constexpr bool is_input2_new_leaf =
                    !has_type2<Input2, LeavesAlive...>();

                if constexpr (is_input1_new_leaf && is_input2_new_leaf) {
                    if constexpr (std::is_same_v<Input1, Input2>) {
                        return tape_size_aux(args<TypesAlive...>{},
                                             args<Input1, LeavesAlive...>{},
                                             args<Leaves...>{});

                    } else {
                        return tape_size_aux(
                            args<TypesAlive...>{},
                            args<Input1, Input2, LeavesAlive...>{},
                            args<Leaves...>{});
                    }
                } else if constexpr (is_input1_new_leaf) {
                    return tape_size_aux(args<TypesAlive...>{},
                                         args<Input1, LeavesAlive...>{},
                                         args<Leaves...>{});
                } else if constexpr (is_input2_new_leaf) {
                    return tape_size_aux(args<TypesAlive...>{},
                                         args<Input2, LeavesAlive...>{},
                                         args<Leaves...>{});
                } else {
                    return tape_size_aux(args<TypesAlive...>{},
                                         args<LeavesAlive...>{},
                                         args<Leaves...>{});
                }
            } else if constexpr (is_input1_leaf) {
                constexpr bool is_input1_new_leaf =
                    !has_type2<Input1, LeavesAlive...>();
                constexpr bool is_input2_new =
                    !has_type2<Input2, TypesAlive...>();

                if constexpr (is_input1_new_leaf && is_input2_new) {
                    return tape_size_aux(args<Input2, TypesAlive...>{},
                                         args<Input1, LeavesAlive...>{},
                                         args<Leaves...>{});
                } else if constexpr (is_input1_new_leaf) {
                    return tape_size_aux(args<TypesAlive...>{},
                                         args<Input1, LeavesAlive...>{},
                                         args<Leaves...>{});
                } else if constexpr (is_input2_new) {
                    return tape_size_aux(args<Input2, TypesAlive...>{},
                                         args<LeavesAlive...>{},
                                         args<Leaves...>{});
                } else {
                    return tape_size_aux(args<TypesAlive...>{},
                                         args<LeavesAlive...>{},
                                         args<Leaves...>{});
                }
            } else if constexpr (is_input2_leaf) {
                constexpr bool is_input1_new =
                    !has_type2<Input1, TypesAlive...>();
                constexpr bool is_input2_new_leaf =
                    !has_type2<Input2, LeavesAlive...>();

                if constexpr (is_input1_new && is_input2_new_leaf) {
                    return tape_size_aux(args<Input1, TypesAlive...>{},
                                         args<Input2, LeavesAlive...>{},
                                         args<Leaves...>{});
                } else if constexpr (is_input1_new) {
                    return tape_size_aux(args<Input1, TypesAlive...>{},
                                         args<LeavesAlive...>{},
                                         args<Leaves...>{});
                } else if constexpr (is_input2_new_leaf) {
                    return tape_size_aux(args<TypesAlive...>{},
                                         args<Input2, LeavesAlive...>{},
                                         args<Leaves...>{});
                } else {
                    return tape_size_aux(args<TypesAlive...>{},
                                         args<LeavesAlive...>{},
                                         args<Leaves...>{});
                }
            } else {
                constexpr bool is_input1_new =
                    !has_type2<Input1, TypesAlive...>();
                constexpr bool is_input2_new =
                    !has_type2<Input2, TypesAlive...>();

                if constexpr (is_input1_new && is_input2_new) {
                    if constexpr (std::is_same_v<Input1, Input2>) {
                        return tape_size_aux(args<Input1, TypesAlive...>{},
                                             args<LeavesAlive...>{},
                                             args<Leaves...>{});
                    } else {
                        return tape_size_aux(
                            args<Input1, Input2, TypesAlive...>{},
                            args<LeavesAlive...>{}, args<Leaves...>{});
                    }
                } else if constexpr (is_input1_new) {
                    return tape_size_aux(args<Input1, TypesAlive...>{},
                                         args<LeavesAlive...>{},
                                         args<Leaves...>{});
                } else if constexpr (is_input2_new) {
                    return tape_size_aux(args<Input2, TypesAlive...>{},
                                         args<LeavesAlive...>{},
                                         args<Leaves...>{});
                } else {
                    return tape_size_aux(args<TypesAlive...>{},
                                         args<LeavesAlive...>{},
                                         args<Leaves...>{});
                }
            }
        } else if constexpr (input1_has_0_deriv) {
            return tape_size_next_univariate_noskip<Input2>(
                args<this_type, TypesAlive...>{}, args<LeavesAlive...>{},
                args<Leaves...>{});
        } else if constexpr (input2_has_0_deriv) {
            return tape_size_next_univariate_noskip<Input1>(
                args<this_type, TypesAlive...>{}, args<LeavesAlive...>{},
                args<Leaves...>{});
        }
    }
}

template <class Input1, class Input2, class this_type, typename... TypesAlive,
          typename... LeavesAlive, typename... Leaves>
constexpr static auto
tape_size_bivariate(args<this_type, TypesAlive...> const &,
                    args<LeavesAlive...> const &, args<Leaves...> const &)
    -> std::size_t {

    static_assert(!has_type2<this_type, Leaves...>());

    constexpr std::size_t curent_tape_size =
        1 + sizeof...(TypesAlive) + sizeof...(LeavesAlive);

    constexpr std::size_t next_tape_size =
        tape_size_next_bivariate<Input1, Input2>(
            args<this_type const, TypesAlive...>{}, args<LeavesAlive...>{},
            args<Leaves...>{});

    return std::max(curent_tape_size, next_tape_size);
}

template <class Input, class this_type, typename... TypesAlive,
          typename... LeavesAlive, typename... Leaves>
constexpr static auto
tape_size_univariate(args<this_type, TypesAlive...> const &,
                     args<LeavesAlive...> const &, args<Leaves...> const &)
    -> std::size_t {

    static_assert(!has_type2<this_type, Leaves...>());

    constexpr std::size_t curent_tape_size =
        1 + sizeof...(TypesAlive) + sizeof...(LeavesAlive);

    constexpr std::size_t next_tape_size = tape_size_next_univariate<Input>(
        args<this_type const, TypesAlive...>{}, args<LeavesAlive...>{},
        args<Leaves...>{});

    return std::max(curent_tape_size, next_tape_size);
}

// template <class Input, class this_type>
// auto is_my_parent_impl(const Univariate<Input, this_type> *) ->
// std::true_type;

// auto is_my_parent_impl(const void *) -> std::false_type;

// template <class Input, class this_type>
// using is_my_parent =
//     decltype(is_my_parent_impl<Input, this_type>(std::declval<this_type
//     *>()));

// template <class Input, class this_type, typename... TypesAlive,
//           typename... LeavesAlive, typename... Leaves>
// constexpr static auto
// tape_size(args<this_type const, TypesAlive...> const &,
//           args<LeavesAlive...> const &, args<Leaves...> const &)
//     -> std::size_t {
//     return tape_size_univariate<Input>(args<this_type const,
//     TypesAlive...>{},
//                                        args<LeavesAlive...>{},
//                                        args<Leaves...>{});
// }

template <class Input, typename... TypesAlive, typename... LeavesAlive,
          typename... Leaves>
constexpr static auto tape_size(args<exp_t<Input> const, TypesAlive...> const &,
                                args<LeavesAlive...> const &,
                                args<Leaves...> const &) -> std::size_t {

    return tape_size_univariate<Input>(
        args<exp_t<Input> const, TypesAlive...>{}, args<LeavesAlive...>{},
        args<Leaves...>{});
}

template <class Input, typename... TypesAlive, typename... LeavesAlive,
          typename... Leaves>
constexpr static auto
tape_size(args<sqrt_t<Input> const, TypesAlive...> const &,
          args<LeavesAlive...> const &, args<Leaves...> const &)
    -> std::size_t {

    return tape_size_univariate<Input>(
        args<sqrt_t<Input> const, TypesAlive...>{}, args<LeavesAlive...>{},
        args<Leaves...>{});
}

template <class Input, typename... TypesAlive, typename... LeavesAlive,
          typename... Leaves>
constexpr static auto tape_size(args<log_t<Input> const, TypesAlive...> const &,
                                args<LeavesAlive...> const &,
                                args<Leaves...> const &) -> std::size_t {

    return tape_size_univariate<Input>(
        args<log_t<Input> const, TypesAlive...>{}, args<LeavesAlive...>{},
        args<Leaves...>{});
}

template <class Input, typename... TypesAlive, typename... LeavesAlive,
          typename... Leaves>
constexpr static auto tape_size(args<cos_t<Input> const, TypesAlive...> const &,
                                args<LeavesAlive...> const &,
                                args<Leaves...> const &) -> std::size_t {

    return tape_size_univariate<Input>(
        args<cos_t<Input> const, TypesAlive...>{}, args<LeavesAlive...>{},
        args<Leaves...>{});
}

template <class Input, typename... TypesAlive, typename... LeavesAlive,
          typename... Leaves>
constexpr static auto
tape_size(args<erfc_t<Input> const, TypesAlive...> const &,
          args<LeavesAlive...> const &, args<Leaves...> const &)
    -> std::size_t {

    return tape_size_univariate<Input>(
        args<erfc_t<Input> const, TypesAlive...>{}, args<LeavesAlive...>{},
        args<Leaves...>{});
}

template <class Input1, class Input2, typename... TypesAlive,
          typename... LeavesAlive, typename... Leaves>
constexpr static auto
tape_size(args<add<Input1, Input2> const, TypesAlive...> const &,
          args<LeavesAlive...> const &, args<Leaves...> const &)
    -> std::size_t {

    return tape_size_bivariate<Input1, Input2>(
        args<add<Input1, Input2> const, TypesAlive...>{},
        args<LeavesAlive...>{}, args<Leaves...>{});
}

template <class Input1, class Input2, typename... TypesAlive,
          typename... LeavesAlive, typename... Leaves>
constexpr static auto
tape_size(args<mul<Input1, Input2> const, TypesAlive...> const &,
          args<LeavesAlive...> const &, args<Leaves...> const &)
    -> std::size_t {

    return tape_size_bivariate<Input1, Input2>(
        args<mul<Input1, Input2> const, TypesAlive...>{},
        args<LeavesAlive...>{}, args<Leaves...>{});
}

template <class Input1, class Input2, typename... TypesAlive,
          typename... LeavesAlive, typename... Leaves>
constexpr static auto
tape_size(args<div<Input1, Input2> const, TypesAlive...> const &,
          args<LeavesAlive...> const &, args<Leaves...> const &)
    -> std::size_t {

    return tape_size_bivariate<Input1, Input2>(
        args<div<Input1, Input2> const, TypesAlive...>{},
        args<LeavesAlive...>{}, args<Leaves...>{});
}

} // namespace detail

template <class Output, typename... Leaves>
constexpr static auto tape_size(Output const &, Leaves const &...)
    -> std::size_t {
    return detail::tape_size(args<Output const>{}, args<>{},
                             args<Leaves const...>{});
}

} // namespace adhoc2

#endif // TAPE_SIZE_HPP