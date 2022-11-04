#ifndef CONSTANTS_TAPE_INPUT_HPP
#define CONSTANTS_TAPE_INPUT_HPP

#include <dependency.hpp>

#include <array>
#include <type_traits>

namespace adhoc2 {

template <class... Inputs> class TapeInput {
    std::array<double, sizeof...(Inputs)> buff{};

    template <class FirstInput, class... InputsIn>
    constexpr static auto check_all_different() -> bool {
        if constexpr (sizeof...(InputsIn) == 0) {
            return true;
        } else {
            return (!std::is_same_v<FirstInput, InputsIn> && ...) &&
                   (!std::is_same_v<FirstInput, InputsIn const> && ...) &&
                   check_all_different<InputsIn...>();
        }
    }

    template <class FirstInput, class... InputsIn>
    constexpr static auto check_all_independent() -> bool {
        if constexpr (sizeof...(InputsIn) == 0) {
            return (!depends<FirstInput, Inputs>::call() && ...);
        } else {
            return (!depends<FirstInput, Inputs>::call() && ...) &&
                   check_all_independent<InputsIn...>();
        }
    }

    TapeInput() = default;

  public:
    explicit TapeInput(Inputs... /* in */) {
        static_assert(check_all_different<Inputs...>(),
                      "you have redundant inputs on the tape");
        static_assert(check_all_independent<Inputs...>(),
                      "you have dependent inputs on the tape");
    }

    template <class ThisInput> void set(ThisInput, double in) {
        static_assert((std::is_same_v<ThisInput, Inputs> || ...),
                      "you are trying to set a variable that is not registered "
                      "on the tape");

        constexpr auto position_on_buff = idx_type2<ThisInput, Inputs...>();
        this->buff[position_on_buff] = in;
    }

    template <class... InputPrev, class... InputNew>
    friend auto addvar(TapeInput<InputPrev...>, InputNew...)
        -> TapeInput<InputPrev..., InputNew...>;
};

template <class... InputPrev, class... InputNew>
[[nodiscard]] auto addvar(TapeInput<InputPrev...> prevtape, InputNew...)
    -> TapeInput<InputPrev..., InputNew...> {
    auto res = TapeInput<InputPrev..., InputNew...>();

    for (std::size_t i = 0; i < prevtape.buff.size(); ++i) {
        res.buff[i] = prevtape.buff[i];
    }
    return res;
}

} // namespace adhoc2

#endif // CONSTANTS_TAPE_INPUT_HPP