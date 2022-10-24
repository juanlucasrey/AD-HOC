#ifndef TAPE_SIZE_HPP
#define TAPE_SIZE_HPP

#include <adhoc2.hpp>

namespace adhoc2::detail {

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

    using this_type = adouble_aux<N>;

    // because we are on a leave, thus_type should be included
    static_assert(has_type2<this_type, Leaves...>());

    // this_type should only be once in TypesAlive
    static_assert(!has_type2<this_type, TypesAlive...>());

    constexpr std::size_t curent_tape_size =
        1 + sizeof...(TypesAlive) + sizeof...(LeavesAlive);
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

template <class Input1, class Input2, class this_type, typename... TypesAlive,
          typename... LeavesAlive, typename... Leaves>
constexpr static auto
tape_size_next_bivariate(args<this_type, TypesAlive...> const &,
                         args<LeavesAlive...> const &, args<Leaves...> const &)
    -> std::size_t {

    constexpr bool other_types_depend_on_this =
        static_cast<bool>(this_type::template depends3<TypesAlive...>());

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
                    return tape_size_aux(args<TypesAlive...>{},
                                         args<Input1, Input2, LeavesAlive...>{},
                                         args<Leaves...>{});
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
                    return tape_size_aux(args<Input1, Input2, TypesAlive...>{},
                                         args<LeavesAlive...>{},
                                         args<Leaves...>{});
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
            if constexpr (is_input2_leaf) {
                constexpr bool is_input2_new_leaf =
                    !has_type2<Input2, LeavesAlive...>();

                if constexpr (is_input2_new_leaf) {
                    return tape_size_aux(args<TypesAlive...>{},
                                         args<Input2, LeavesAlive...>{},
                                         args<Leaves...>{});
                } else {
                    return tape_size_aux(args<TypesAlive...>{},
                                         args<LeavesAlive...>{},
                                         args<Leaves...>{});
                }
            } else {
                constexpr bool is_input2_new =
                    !has_type2<Input2, TypesAlive...>();

                if constexpr (is_input2_new) {
                    return tape_size_aux(args<Input2, TypesAlive...>{},
                                         args<LeavesAlive...>{},
                                         args<Leaves...>{});
                } else {
                    return tape_size_aux(args<TypesAlive...>{},
                                         args<LeavesAlive...>{},
                                         args<Leaves...>{});
                }
            }
        } else if constexpr (input2_has_0_deriv) {
            if constexpr (is_input1_leaf) {
                constexpr bool is_input1_new_leaf =
                    !has_type2<Input1, LeavesAlive...>();

                if constexpr (is_input1_new_leaf) {
                    return tape_size_aux(args<TypesAlive...>{},
                                         args<Input1, LeavesAlive...>{},
                                         args<Leaves...>{});
                } else {
                    return tape_size_aux(args<TypesAlive...>{},
                                         args<LeavesAlive...>{},
                                         args<Leaves...>{});
                }
            } else {
                constexpr bool is_input1_new =
                    !has_type2<Input1, TypesAlive...>();

                if constexpr (is_input1_new) {
                    return tape_size_aux(args<Input1, TypesAlive...>{},
                                         args<LeavesAlive...>{},
                                         args<Leaves...>{});
                } else {
                    return tape_size_aux(args<TypesAlive...>{},
                                         args<LeavesAlive...>{},
                                         args<Leaves...>{});
                }
            }
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
            args<mul_active<Input1, Input2> const, TypesAlive...>{},
            args<LeavesAlive...>{}, args<Leaves...>{});

    return std::max(curent_tape_size, next_tape_size);
}

template <class Input1, class Input2, typename... TypesAlive,
          typename... LeavesAlive, typename... Leaves>
constexpr static auto
tape_size(args<mul_active<Input1, Input2> const, TypesAlive...> const &,
          args<LeavesAlive...> const &, args<Leaves...> const &)
    -> std::size_t {

    return tape_size_bivariate<Input1, Input2>(
        args<mul_active<Input1, Input2> const, TypesAlive...>{},
        args<LeavesAlive...>{}, args<Leaves...>{});
}

template <class Input1, class Input2, typename... TypesAlive,
          typename... LeavesAlive, typename... Leaves>
constexpr static auto
tape_size(args<add_active<Input1, Input2> const, TypesAlive...> const &,
          args<LeavesAlive...> const &, args<Leaves...> const &)
    -> std::size_t {

    return tape_size_bivariate<Input1, Input2>(
        args<add_active<Input1, Input2> const, TypesAlive...>{},
        args<LeavesAlive...>{}, args<Leaves...>{});
}

} // namespace adhoc2::detail

#endif // TAPE_SIZE_HPP