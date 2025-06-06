#ifndef BRNG_DISTRIBUTION_RNGTEST
#define BRNG_DISTRIBUTION_RNGTEST

template <class UIntType, UIntType Min, UIntType Max> struct RNGtest final {

    using result_type = UIntType;

    explicit RNGtest(result_type value_in) : value(value_in) {}

    template <bool FwdDirection = true>
    inline auto operator()() -> result_type {
        return value;
    }

    static constexpr auto min() -> UIntType { return Min; }
    static constexpr auto max() -> UIntType { return Max; }

  private:
    UIntType value = Min;
};

#endif // BRNG_DISTRIBUTION_RNGTEST
