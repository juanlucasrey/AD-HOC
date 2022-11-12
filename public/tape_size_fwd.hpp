#ifndef TAPE_SIZE_FWD_HPP
#define TAPE_SIZE_FWD_HPP

#include <adhoc2.hpp>
#include <adhoc3.hpp>
#include <constants_type3.hpp>
#include <dependency.hpp>
#include <dependency3.hpp>

namespace adhoc2 {
namespace detail {

template <typename... TypesAlive, typename... LeavesAlive, typename... Leaves>
constexpr static auto tape_size_fwd(const args<TypesAlive...> &,
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
            return tape_size_fwd(args<TypesAlive...>{},
                                 args<Input, LeavesAlive...>{},
                                 args<Leaves...>{});
        } else {
            return tape_size_fwd(args<Input, TypesAlive...>{},
                                 args<LeavesAlive...>{}, args<Leaves...>{});
        }
    } else {
        return tape_size_fwd(args<TypesAlive...>{}, args<LeavesAlive...>{},
                             args<Leaves...>{});
    }
}

template <class Input1, class Input2, class this_type, typename... TypesAlive,
          typename... LeavesAlive, typename... Leaves>
constexpr static auto
tape_size_bivariate(args<this_type, TypesAlive...> const &,
                    args<LeavesAlive...> const &, args<Leaves...> const &)
    -> std::size_t {
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
                    return tape_size_fwd(args<TypesAlive...>{},
                                         args<Input1, LeavesAlive...>{},
                                         args<Leaves...>{});

                } else {
                    return tape_size_fwd(args<TypesAlive...>{},
                                         args<Input1, Input2, LeavesAlive...>{},
                                         args<Leaves...>{});
                }
            } else if constexpr (is_input1_new_leaf) {
                return tape_size_fwd(args<TypesAlive...>{},
                                     args<Input1, LeavesAlive...>{},
                                     args<Leaves...>{});
            } else if constexpr (is_input2_new_leaf) {
                return tape_size_fwd(args<TypesAlive...>{},
                                     args<Input2, LeavesAlive...>{},
                                     args<Leaves...>{});
            } else {
                return tape_size_fwd(args<TypesAlive...>{},
                                     args<LeavesAlive...>{}, args<Leaves...>{});
            }
        } else if constexpr (is_input1_leaf) {
            constexpr bool is_input1_new_leaf =
                !has_type2<Input1, LeavesAlive...>();
            constexpr bool is_input2_new = !has_type2<Input2, TypesAlive...>();

            if constexpr (is_input1_new_leaf && is_input2_new) {
                return tape_size_fwd(args<Input2, TypesAlive...>{},
                                     args<Input1, LeavesAlive...>{},
                                     args<Leaves...>{});
            } else if constexpr (is_input1_new_leaf) {
                return tape_size_fwd(args<TypesAlive...>{},
                                     args<Input1, LeavesAlive...>{},
                                     args<Leaves...>{});
            } else if constexpr (is_input2_new) {
                return tape_size_fwd(args<Input2, TypesAlive...>{},
                                     args<LeavesAlive...>{}, args<Leaves...>{});
            } else {
                return tape_size_fwd(args<TypesAlive...>{},
                                     args<LeavesAlive...>{}, args<Leaves...>{});
            }
        } else if constexpr (is_input2_leaf) {
            constexpr bool is_input1_new = !has_type2<Input1, TypesAlive...>();
            constexpr bool is_input2_new_leaf =
                !has_type2<Input2, LeavesAlive...>();

            if constexpr (is_input1_new && is_input2_new_leaf) {
                return tape_size_fwd(args<Input1, TypesAlive...>{},
                                     args<Input2, LeavesAlive...>{},
                                     args<Leaves...>{});
            } else if constexpr (is_input1_new) {
                return tape_size_fwd(args<Input1, TypesAlive...>{},
                                     args<LeavesAlive...>{}, args<Leaves...>{});
            } else if constexpr (is_input2_new_leaf) {
                return tape_size_fwd(args<TypesAlive...>{},
                                     args<Input2, LeavesAlive...>{},
                                     args<Leaves...>{});
            } else {
                return tape_size_fwd(args<TypesAlive...>{},
                                     args<LeavesAlive...>{}, args<Leaves...>{});
            }
        } else {
            constexpr bool is_input1_new = !has_type2<Input1, TypesAlive...>();
            constexpr bool is_input2_new = !has_type2<Input2, TypesAlive...>();

            if constexpr (is_input1_new && is_input2_new) {
                if constexpr (std::is_same_v<Input1, Input2>) {
                    return tape_size_fwd(args<Input1, TypesAlive...>{},
                                         args<LeavesAlive...>{},
                                         args<Leaves...>{});
                } else {
                    return tape_size_fwd(args<Input1, Input2, TypesAlive...>{},
                                         args<LeavesAlive...>{},
                                         args<Leaves...>{});
                }
            } else if constexpr (is_input1_new) {
                return tape_size_fwd(args<Input1, TypesAlive...>{},
                                     args<LeavesAlive...>{}, args<Leaves...>{});
            } else if constexpr (is_input2_new) {
                return tape_size_fwd(args<Input2, TypesAlive...>{},
                                     args<LeavesAlive...>{}, args<Leaves...>{});
            } else {
                return tape_size_fwd(args<TypesAlive...>{},
                                     args<LeavesAlive...>{}, args<Leaves...>{});
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
constexpr static auto tape_size_fwd(const args<> &,
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

    return tape_size_fwd(args<NextType, TypesAlive..., CurrentType>{},
                         args<LeavesAlive...>{}, args<Leaves...>{});
}

template <template <class...> class Xvariate, class... Input,
          typename... TypesAlive, typename... LeavesAlive, typename... Leaves>
constexpr static auto
tape_size_fwd(args<Xvariate<Input...> const, TypesAlive...> const &,
              args<LeavesAlive...> const &, args<Leaves...> const &)
    -> std::size_t {

    using this_type = Xvariate<Input...>;
    static_assert(!has_type2<this_type, Leaves...>());

    constexpr bool other_types_depend_on_this =
        (depends<TypesAlive, this_type>::call() || ...);

    if constexpr (other_types_depend_on_this) {
        return skip_type(args<this_type const, TypesAlive...>{},
                         args<LeavesAlive...>{}, args<Leaves...>{}) +
               1;
    } else {
        static_assert(sizeof...(Input) == 1 || sizeof...(Input) == 2);
        if constexpr (sizeof...(Input) == 1) {
            return tape_size_univariate<Input...>(
                       args<this_type, TypesAlive...>{}, args<LeavesAlive...>{},
                       args<Leaves...>{}) +
                   1;
        } else {
            return tape_size_bivariate<Input...>(
                       args<this_type, TypesAlive...>{}, args<LeavesAlive...>{},
                       args<Leaves...>{}) +
                   1;
        }
    }
}

// note: this function is only used in the trivial case of a tape size for a
// simple variable output because all other cases are xvariate functions,
// and, if an adhoc type is to be found, it's either an active leaf or
// passive leaf (which are never treated on their own)
template <int N, typename... Leaves>
constexpr static auto tape_size_fwd(const args<adouble_aux<N> const> &,
                                    const args<> &, const args<Leaves...> &)
    -> std::size_t {
    return has_type2<adouble_aux<N> const, Leaves...>();
}

} // namespace detail

template <class Output, typename... Leaves>
constexpr static auto tape_size_fwd(Output const &, Leaves const &...)
    -> std::size_t {
    return detail::tape_size_fwd(args<Output const>{}, args<>{},
                                 args<Leaves const...>{});
}

namespace detail {

template <typename... TypesAlive>
constexpr auto tape_size_fwd2(const args<TypesAlive...> &) -> std::size_t;

template <> constexpr auto tape_size_fwd2(const args<> &) -> std::size_t {
    return 0;
}

template <int N, typename... TypesAlive>
constexpr auto tape_size_fwd2(const args<adouble_aux<N> const, TypesAlive...> &)
    -> std::size_t {
    return tape_size_fwd2(
        args<TypesAlive...>{}); // we don't add anything because it's on the
                                // input tape
}

template <typename CurrentType, typename NextType, typename... TypesAlive>
constexpr auto skip_type2(const args<CurrentType, NextType, TypesAlive...> &)
    -> std::size_t {

    return tape_size_fwd2(
        args<NextType,
             TypesAlive...>{}); // current type will come up again because it is
                                // included on one of the other types
}

template <template <class...> class Xvariate, class... Input,
          typename... TypesAlive>
constexpr auto
tape_size_fwd2(const args<Xvariate<Input...> const, TypesAlive...> &)
    -> std::size_t {
    using this_type = Xvariate<Input...>;
    static_assert(!has_type2<this_type, TypesAlive...>());

    constexpr bool other_types_depend_on_this =
        (depends<TypesAlive, this_type>::call() || ...);

    if constexpr (other_types_depend_on_this) {
        return skip_type2(args<this_type const, TypesAlive...>{});
    } else {
        return tape_size_fwd2(args<Input..., TypesAlive...>{}) + 1;
    }
}

} // namespace detail

template <class Output>
constexpr auto tape_size_fwd2(Output const &) -> std::size_t {
    return detail::tape_size_fwd2(args<Output const>{});
}

} // namespace adhoc2

namespace adhoc3 {

namespace detail {

template <typename... TypesAlive>
constexpr auto fwd_calc_order_t(const std::tuple<TypesAlive...> &);

template <> constexpr auto fwd_calc_order_t(const std::tuple<> &) {
    return std::tuple<>{};
}

#if __cplusplus >= 202002L

template <constants::detail::AsTemplateArg<double> D, typename... TypesAlive>
constexpr auto
fwd_calc_order_t(const std::tuple<constants::CD<D> const, TypesAlive...> &) {
    return fwd_calc_order_t(
        std::tuple<TypesAlive...>{}); // we don't add anything because it's on
                                      // the input tape
}

template <constants::detail::AsTemplateArg<double> D, typename... TypesAlive>
constexpr auto
fwd_calc_order_t(const std::tuple<constants::CD<D>, TypesAlive...> &) {
    return fwd_calc_order_t(
        std::tuple<TypesAlive...>{}); // we don't add anything because it's on
                                      // the input tape
}

#else

template <std::uint64_t D, typename... TypesAlive>
constexpr auto
fwd_calc_order_t(const std::tuple<constants::CD<D> const, TypesAlive...> &) {
    return fwd_calc_order_t(
        std::tuple<TypesAlive...>{}); // we don't add anything because it's on
                                      // the input tape
}

template <std::uint64_t D, typename... TypesAlive>
constexpr auto
fwd_calc_order_t(const std::tuple<constants::CD<D>, TypesAlive...> &) {
    return fwd_calc_order_t(
        std::tuple<TypesAlive...>{}); // we don't add anything because it's on
                                      // the input tape
}

#endif

template <std::size_t N, typename... TypesAlive>
constexpr auto
fwd_calc_order_t(const std::tuple<adouble_aux<N> const, TypesAlive...> &) {
    return fwd_calc_order_t(
        std::tuple<TypesAlive...>{}); // we don't add anything because it's on
                                      // the input tape
}

template <std::size_t N, typename... TypesAlive>
constexpr auto
fwd_calc_order_t(const std::tuple<adouble_aux<N>, TypesAlive...> &) {
    return fwd_calc_order_t(
        std::tuple<TypesAlive...>{}); // we don't add anything because it's on
                                      // the input tape
}

template <typename CurrentType, typename NextType, typename... TypesAlive>
constexpr auto
skip_type(const std::tuple<CurrentType, NextType, TypesAlive...> &) {
    return fwd_calc_order_t(
        std::tuple<NextType, TypesAlive...>{}); // current type will come up
                                                // again because it is included
                                                // on one of the other types
}

template <template <class...> class Xvariate, class... Input,
          typename... TypesAlive>
constexpr auto
fwd_calc_order_t(const std::tuple<Xvariate<Input...> const, TypesAlive...> &) {
    using this_type = Xvariate<Input...>;
    static_assert(!adhoc2::has_type2<this_type, TypesAlive...>());

    constexpr bool other_types_depend_on_this =
        (depends<TypesAlive, this_type>::call() || ...);

    if constexpr (other_types_depend_on_this) {
        return std::tuple<decltype(skip_type(
            std::tuple<this_type const, TypesAlive...>{}))>{};
    } else {
        return std::tuple<decltype(fwd_calc_order_t(
                              std::tuple<Input..., TypesAlive...>{})),
                          this_type>{};
    }
}

template <typename T> struct flatten_tuple { using type = std::tuple<T>; };

template <typename T> using flatten_tuple_t = typename flatten_tuple<T>::type;

template <typename... Ts> struct flatten_tuple<std::tuple<Ts...>> {
    using type =
        decltype(std::tuple_cat(std::declval<flatten_tuple_t<Ts>>()...));
};

} // namespace detail

template <class Output> constexpr auto fwd_calc_order_t(Output const &) {
    return detail::flatten_tuple_t<decltype(detail::fwd_calc_order_t(
        std::tuple<Output const>{}))>{};
}

} // namespace adhoc3

#endif // TAPE_SIZE_FWD_HPP