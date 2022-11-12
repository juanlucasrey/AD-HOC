#ifndef EVALUATE_FWD_HPP
#define EVALUATE_FWD_HPP

#include <constants_type3.hpp>
#include <tape_input3.hpp>
#include <tape_intermediate.hpp>
#include <tape_size_fwd.hpp>

namespace adhoc3 {

namespace detail {

template <class Val, class... InputTypes, class... IntermediateTypes>
inline auto get(TapeIntermediate<InputTypes...> const &in,
                TapeIntermediate<IntermediateTypes...> &intermediate) -> double;

template <std::uint64_t D, class... InputTypes, class... IntermediateTypes>
inline auto get(TapeIntermediate<InputTypes...> const &,
                TapeIntermediate<IntermediateTypes...> &, constants::CD<D> &)
    -> double {
    return constants::CD<D>::v();
}

template <std::uint64_t D, class... InputTypes, class... IntermediateTypes>
inline auto get(TapeIntermediate<InputTypes...> const &,
                TapeIntermediate<IntermediateTypes...> &,
                constants::CD<D> const &) -> double {
    return constants::CD<D>::v();
}

template <class Input, class... InputTypes, class... IntermediateTypes>
inline auto get(TapeIntermediate<InputTypes...> const &in,
                TapeIntermediate<IntermediateTypes...> &intermediate, Input &)
    -> double {
    constexpr bool is_input = adhoc2::has_type2<Input, InputTypes...>();
    if constexpr (is_input) {
        return in.get(Input{});
    } else {
        return intermediate.get(Input{});
    }
}

template <class... InputTypes, class... IntermediateTypes>
inline void
evaluate_fwd(TapeIntermediate<InputTypes...> const & /* in */,
             TapeIntermediate<IntermediateTypes...> & /* intermediate */,
             std::tuple<>) {}

template <class... InputTypes, class... IntermediateTypes,
          template <class> class Univariate, class Input,
          class... IntermediateTypesToCalc>
inline void
evaluate_fwd(TapeIntermediate<InputTypes...> const &in,
             TapeIntermediate<IntermediateTypes...> &intermediate,
             std::tuple<Univariate<Input>, IntermediateTypesToCalc...>) {
    intermediate.set(Univariate<Input>{},
                     Univariate<Input>::v(get(in, intermediate, Input{})));
    evaluate_fwd(in, intermediate, std::tuple<IntermediateTypesToCalc...>{});
}

template <class... InputTypes, class... IntermediateTypes,
          template <class, class> class Bivariate, class Input1, class Input2,
          class... IntermediateTypesToCalc>
inline void evaluate_fwd(
    TapeIntermediate<InputTypes...> const &in,
    TapeIntermediate<IntermediateTypes...> &intermediate,
    std::tuple<Bivariate<Input1, Input2>, IntermediateTypesToCalc...>) {
    intermediate.set(
        Bivariate<Input1, Input2>{},
        Bivariate<Input1, Input2>::v(get(in, intermediate, Input1{}),
                                     get(in, intermediate, Input2{})));
    evaluate_fwd(in, intermediate, std::tuple<IntermediateTypesToCalc...>{});
}

template <class... InputTypes, class... IntermediateTypes>
inline void
evaluate_fwd_first(TapeIntermediate<InputTypes...> const &in,
                   TapeIntermediate<IntermediateTypes...> &intermediate) {
    evaluate_fwd(in, intermediate, std::tuple<IntermediateTypes...>{});
}

} // namespace detail

template <class... InputTypes, class... IntermediateTypes>
inline void
evaluate_fwd(detail::TapeIntermediate<InputTypes...> const &in,
             detail::TapeIntermediate<IntermediateTypes...> &intermediate) {
    detail::evaluate_fwd_first(in, intermediate);
}

} // namespace adhoc3

#endif // EVALUATE_FWD_HPP