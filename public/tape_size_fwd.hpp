#ifndef NEW_TAPE_SIZE_FWD_HPP
#define NEW_TAPE_SIZE_FWD_HPP

#include "adhoc.hpp"
#include "constants_type.hpp"
#include "dependency.hpp"
#include "utils.hpp"

#include <tuple>

namespace adhoc {

namespace detail {

template <typename... TypesAlive> struct fwd_calc_order_aux_t3 {
    template <typename... Operations> constexpr static auto call() noexcept;
};

template <> struct fwd_calc_order_aux_t3<> {
    template <typename... Operations> constexpr static auto call() noexcept {
        return std::tuple<Operations...>{};
    }
};

#if __cplusplus >= 202002L

template <constants::detail::AsTemplateArg<double> D, typename... TypesAlive>
struct fwd_calc_order_aux_t3<constants::CD<D>, TypesAlive...> {
    template <typename... Operations> constexpr static auto call() noexcept {
        // we don't add anything to operations because it's a constant
        return fwd_calc_order_aux_t3<TypesAlive...>::template call<
            Operations...>();
    }
};

template <constants::detail::AsTemplateArg<double> D, typename... TypesAlive>
struct fwd_calc_order_aux_t3<constants::CD<D> const, TypesAlive...> {
    template <typename... Operations> constexpr static auto call() noexcept {
        // we don't add anything to operations because it's a constant
        return fwd_calc_order_aux_t3<TypesAlive...>::template call<
            Operations...>();
    }
};

#else

template <std::uint64_t D, typename... TypesAlive>
struct fwd_calc_order_aux_t3<constants::CD<D>, TypesAlive...> {
    template <typename... Operations> constexpr static auto call() noexcept {
        // we don't add anything to operations because it's a constant
        return fwd_calc_order_aux_t3<TypesAlive...>::template call<
            Operations...>();
    }
};

template <std::uint64_t D, typename... TypesAlive>
struct fwd_calc_order_aux_t3<constants::CD<D> const, TypesAlive...> {
    template <typename... Operations> constexpr static auto call() noexcept {
        // we don't add anything to operations because it's a constant
        return fwd_calc_order_aux_t3<TypesAlive...>::template call<
            Operations...>();
    }
};

#endif

template <std::size_t N, typename... TypesAlive>
struct fwd_calc_order_aux_t3<double_t<N>, TypesAlive...> {
    template <typename... Operations> constexpr static auto call() noexcept {
        // we don't add anything to operations because it's an input
        return fwd_calc_order_aux_t3<TypesAlive...>::template call<
            Operations...>();
    }
};

template <std::size_t N, typename... TypesAlive>
struct fwd_calc_order_aux_t3<double_t<N> const, TypesAlive...> {
    template <typename... Operations> constexpr static auto call() noexcept {
        // we don't add anything to operations because it's an input
        return fwd_calc_order_aux_t3<TypesAlive...>::template call<
            Operations...>();
    }
};

template <typename CurrentType, typename NextType, typename... TypesAlive>
struct skip_operation {
    template <typename... Operations> constexpr static auto call() noexcept {
        // current type will come up again because it is included on TypesAlive
        return fwd_calc_order_aux_t3<NextType, TypesAlive...>::template call<
            Operations...>();
    }
};

template <template <class...> class Xvariate, class... Input,
          typename... TypesAlive>
struct fwd_calc_order_aux_t3<Xvariate<Input...>, TypesAlive...> {
    template <typename... Operations> constexpr static auto call() noexcept {
        using this_type = Xvariate<Input...>;
        static_assert(!has_type2<this_type, TypesAlive...>());

        constexpr bool other_types_depend_on_this =
            (depends<TypesAlive, this_type>::call() || ...);

        if constexpr (other_types_depend_on_this) {
            return skip_operation<this_type, TypesAlive...>::template call<
                Operations...>();
        } else {
            return fwd_calc_order_aux_t3<Input..., TypesAlive...>::
                template call<this_type, Operations...>();
        }
    }
};

template <template <class...> class Xvariate, class... Input,
          typename... TypesAlive>
struct fwd_calc_order_aux_t3<Xvariate<Input...> const, TypesAlive...> {
    template <typename... Operations> constexpr static auto call() noexcept {
        using this_type = Xvariate<Input...>;
        static_assert(!has_type2<this_type, TypesAlive...>());

        constexpr bool other_types_depend_on_this =
            (depends<TypesAlive, this_type>::call() || ...);

        if constexpr (other_types_depend_on_this) {
            return skip_operation<this_type, TypesAlive...>::template call<
                Operations...>();
        } else {
            return fwd_calc_order_aux_t3<Input..., TypesAlive...>::
                template call<this_type, Operations...>();
        }
    }
};

} // namespace detail

template <class... Outputs>
constexpr auto fwd_calc_order_t(Outputs const &.../* o */) {
    return detail::fwd_calc_order_aux_t3<Outputs...>::template call();
}

} // namespace adhoc

#endif // NEW_TAPE_SIZE_FWD_HPP