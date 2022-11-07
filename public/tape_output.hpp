#ifndef CONSTANTS_TAPE_OUTPUT_HPP
#define CONSTANTS_TAPE_OUTPUT_HPP

#include <dependency.hpp>

#include <array>
#include <type_traits>

namespace adhoc2 {

template <class... Outputs> class TapeOutput {
    std::array<double, sizeof...(Outputs)> buff{};

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
            return (!depends<FirstInput, Outputs>::call() && ...);
        } else {
            return (!depends<FirstInput, Outputs>::call() && ...) &&
                   check_all_independent<InputsIn...>();
        }
    }

    TapeOutput() = default;

  public:
    explicit TapeOutput(Outputs... /* in */) {
        static_assert(check_all_different<Outputs...>(),
                      "you have redundant inputs on the tape");
        static_assert(check_all_independent<Outputs...>(),
                      "you have dependent inputs on the tape");
    }

    template <class ThisInput> void set(ThisInput, double in) {
        static_assert((std::is_same_v<ThisInput, Outputs> || ...),
                      "you are trying to set a variable that is not registered "
                      "on the tape");

        constexpr auto position_on_buff = idx_type2<ThisInput, Outputs...>();
        this->buff[position_on_buff] = in;
    }

    template <class... InputPrev, class... InputNew>
    friend auto addvar(TapeOutput<InputPrev...>, InputNew...)
        -> TapeOutput<InputPrev..., InputNew...>;
};

template <class... InputPrev, class... InputNew>
[[nodiscard]] auto addvar(TapeOutput<InputPrev...> prevtape, InputNew...)
    -> TapeOutput<InputPrev..., InputNew...> {
    auto res = TapeOutput<InputPrev..., InputNew...>();

    for (std::size_t i = 0; i < prevtape.buff.size(); ++i) {
        res.buff[i] = prevtape.buff[i];
    }
    return res;
}

} // namespace adhoc2

#endif // CONSTANTS_TAPE_OUTPUT_HPP