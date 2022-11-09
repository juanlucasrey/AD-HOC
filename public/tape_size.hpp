#ifndef TAPE_SIZE_HPP
#define TAPE_SIZE_HPP

#include <adhoc2.hpp>
#include <dependency.hpp>

namespace adhoc2 {
namespace detail {

template <typename... TypesAlive, typename... LeavesAlive, typename... Leaves>
constexpr static auto tape_size(const args<TypesAlive...> &,
                                const args<LeavesAlive...> &,
                                const args<Leaves...> &) -> std::size_t;

template <class Input, class this_type, typename... TypesAlive,
          typename... LeavesAlive, typename... Leaves>
constexpr static auto
tape_size_univariate(args<this_type, TypesAlive...> const &,
                     args<LeavesAlive...> const &, args<Leaves...> const &)
    -> std::size_t {
    constexpr bool is_input_new =
        !has_type2<Input, TypesAlive..., LeavesAlive...>();
    if constexpr (is_input_new) {
        constexpr bool is_input_leaf = has_type2<Input, Leaves...>();
        if constexpr (is_input_leaf) {
            return tape_size(args<TypesAlive...>{},
                             args<Input, LeavesAlive...>{}, args<Leaves...>{});
        } else {
            return tape_size(args<Input, TypesAlive...>{},
                             args<LeavesAlive...>{}, args<Leaves...>{});
        }
    } else {
        return tape_size(args<TypesAlive...>{}, args<LeavesAlive...>{},
                         args<Leaves...>{});
    }
}

template <class Input, template <class> class Univariate,
          typename... TypesAlive, typename... LeavesAlive, typename... Leaves>
constexpr static auto
tape_size_xvariate(args<Univariate<Input>, TypesAlive...> const &,
                   args<LeavesAlive...> const &, args<Leaves...> const &)
    -> std::size_t {
    using this_type = Univariate<Input>;
    return tape_size_univariate<Input>(args<this_type, TypesAlive...>{},
                                       args<LeavesAlive...>{},
                                       args<Leaves...>{});
}

template <class Input1, class Input2, template <class, class> class Bivariate,
          typename... TypesAlive, typename... LeavesAlive, typename... Leaves>
constexpr static auto
tape_size_xvariate(args<Bivariate<Input1, Input2>, TypesAlive...> const &,
                   args<LeavesAlive...> const &, args<Leaves...> const &)
    -> std::size_t {
    using this_type = Bivariate<Input1, Input2>;
    constexpr bool input1_has_0_deriv =
        !equal_or_depends_many<Input1, Leaves...>();
    constexpr bool input2_has_0_deriv =
        !equal_or_depends_many<Input2, Leaves...>();
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
                    return tape_size(args<TypesAlive...>{},
                                     args<Input1, LeavesAlive...>{},
                                     args<Leaves...>{});

                } else {
                    return tape_size(args<TypesAlive...>{},
                                     args<Input1, Input2, LeavesAlive...>{},
                                     args<Leaves...>{});
                }
            } else if constexpr (is_input1_new_leaf) {
                return tape_size(args<TypesAlive...>{},
                                 args<Input1, LeavesAlive...>{},
                                 args<Leaves...>{});
            } else if constexpr (is_input2_new_leaf) {
                return tape_size(args<TypesAlive...>{},
                                 args<Input2, LeavesAlive...>{},
                                 args<Leaves...>{});
            } else {
                return tape_size(args<TypesAlive...>{}, args<LeavesAlive...>{},
                                 args<Leaves...>{});
            }
        } else if constexpr (is_input1_leaf) {
            constexpr bool is_input1_new_leaf =
                !has_type2<Input1, LeavesAlive...>();
            constexpr bool is_input2_new = !has_type2<Input2, TypesAlive...>();

            if constexpr (is_input1_new_leaf && is_input2_new) {
                return tape_size(args<Input2, TypesAlive...>{},
                                 args<Input1, LeavesAlive...>{},
                                 args<Leaves...>{});
            } else if constexpr (is_input1_new_leaf) {
                return tape_size(args<TypesAlive...>{},
                                 args<Input1, LeavesAlive...>{},
                                 args<Leaves...>{});
            } else if constexpr (is_input2_new) {
                return tape_size(args<Input2, TypesAlive...>{},
                                 args<LeavesAlive...>{}, args<Leaves...>{});
            } else {
                return tape_size(args<TypesAlive...>{}, args<LeavesAlive...>{},
                                 args<Leaves...>{});
            }
        } else if constexpr (is_input2_leaf) {
            constexpr bool is_input1_new = !has_type2<Input1, TypesAlive...>();
            constexpr bool is_input2_new_leaf =
                !has_type2<Input2, LeavesAlive...>();

            if constexpr (is_input1_new && is_input2_new_leaf) {
                return tape_size(args<Input1, TypesAlive...>{},
                                 args<Input2, LeavesAlive...>{},
                                 args<Leaves...>{});
            } else if constexpr (is_input1_new) {
                return tape_size(args<Input1, TypesAlive...>{},
                                 args<LeavesAlive...>{}, args<Leaves...>{});
            } else if constexpr (is_input2_new_leaf) {
                return tape_size(args<TypesAlive...>{},
                                 args<Input2, LeavesAlive...>{},
                                 args<Leaves...>{});
            } else {
                return tape_size(args<TypesAlive...>{}, args<LeavesAlive...>{},
                                 args<Leaves...>{});
            }
        } else {
            constexpr bool is_input1_new = !has_type2<Input1, TypesAlive...>();
            constexpr bool is_input2_new = !has_type2<Input2, TypesAlive...>();

            if constexpr (is_input1_new && is_input2_new) {
                if constexpr (std::is_same_v<Input1, Input2>) {
                    return tape_size(args<Input1, TypesAlive...>{},
                                     args<LeavesAlive...>{}, args<Leaves...>{});
                } else {
                    return tape_size(args<Input1, Input2, TypesAlive...>{},
                                     args<LeavesAlive...>{}, args<Leaves...>{});
                }
            } else if constexpr (is_input1_new) {
                return tape_size(args<Input1, TypesAlive...>{},
                                 args<LeavesAlive...>{}, args<Leaves...>{});
            } else if constexpr (is_input2_new) {
                return tape_size(args<Input2, TypesAlive...>{},
                                 args<LeavesAlive...>{}, args<Leaves...>{});
            } else {
                return tape_size(args<TypesAlive...>{}, args<LeavesAlive...>{},
                                 args<Leaves...>{});
            }
        }
    } else if constexpr (input1_has_0_deriv) {
        return tape_size_univariate<Input2>(args<this_type, TypesAlive...>{},
                                            args<LeavesAlive...>{},
                                            args<Leaves...>{});
    } else if constexpr (input2_has_0_deriv) {
        return tape_size_univariate<Input1>(args<this_type, TypesAlive...>{},
                                            args<LeavesAlive...>{},
                                            args<Leaves...>{});
    }
}

template <typename... LeavesAlive, typename... Leaves>
constexpr static auto tape_size(const args<> &, const args<LeavesAlive...> &,
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

template <template <class...> class Xvariate, class... Input,
          typename... TypesAlive, typename... LeavesAlive, typename... Leaves>
constexpr static auto
tape_size(args<Xvariate<Input...> const, TypesAlive...> const &,
          args<LeavesAlive...> const &, args<Leaves...> const &)
    -> std::size_t {

    using this_type = Xvariate<Input...>;
    static_assert(!has_type2<this_type, Leaves...>());

    constexpr std::size_t curent_tape_size =
        1 + sizeof...(TypesAlive) + sizeof...(LeavesAlive);

    constexpr bool other_types_depend_on_this =
        (depends<TypesAlive, this_type>::call() || ...);

    if constexpr (other_types_depend_on_this) {
        constexpr std::size_t next_tape_size =
            skip_type(args<this_type const, TypesAlive...>{},
                      args<LeavesAlive...>{}, args<Leaves...>{});
        return std::max(curent_tape_size, next_tape_size);
    } else {
        constexpr std::size_t next_tape_size =
            tape_size_xvariate(args<this_type, TypesAlive...>{},
                               args<LeavesAlive...>{}, args<Leaves...>{});
        return std::max(curent_tape_size, next_tape_size);
    }
}

// note: this function is only used in the trivial case of a tape size for a
// simple variable output because all other cases are xvariate functions,
// and, if an adhoc type is to be found, it's either an active leaf or
// passive leaf (which are never treated on their own)
template <int N, typename... Leaves>
constexpr static auto tape_size(const args<adouble_aux<N> const> &,
                                const args<> &, const args<Leaves...> &)
    -> std::size_t {
    return has_type2<adouble_aux<N> const, Leaves...>();
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