#ifndef BRNG_DISTRIBUTION_RNGTEST
#define BRNG_DISTRIBUTION_RNGTEST

template <class UIntType, UIntType Min, UIntType Max> struct RNGtest final {

    using value_type = UIntType;
    using result_type = UIntType;

    explicit RNGtest(result_type value_in) : value(value_in) {}

    inline auto operator*() const -> value_type { return value; }
    inline auto operator++() -> RNGtest & { return *this; }
    inline auto operator--() -> RNGtest & { return *this; }

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
