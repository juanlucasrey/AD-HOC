#ifndef BRNG_DISTRIBUTION_UNIFORM_DISTRIBUTION
#define BRNG_DISTRIBUTION_UNIFORM_DISTRIBUTION

#include <limits>

namespace adhoc {

template <class RealType = double, bool LeftClosed = true,
          bool RightClosed = false, RealType a = 0.0, RealType b = 1.0>
class uniform_real_distribution final {
  public:
    using result_type = RealType;

    template <bool FwdDirection = true, typename RNG>
    auto operator()(RNG &g) -> result_type {
        static_assert(RNG::max() > RNG::min());
        constexpr auto full_range_int = RNG::max() - RNG::min();

        constexpr auto full_range =
            static_cast<RealType>(full_range_int) +
            static_cast<RealType>(
                !RightClosed) + // right interval might be closed
            static_cast<RealType>(!LeftClosed); // left interval might be open

        constexpr RealType one_over_range = (b - a) / full_range;

        constexpr auto stub = static_cast<RNG::result_type>(!LeftClosed);

        if constexpr (RNG::min() == stub) {
            if constexpr (a) {
                return (static_cast<RealType>(
                            g.template operator()<FwdDirection>()) *
                        one_over_range) +
                       a;
            } else {
                return static_cast<RealType>(
                           g.template operator()<FwdDirection>()) *
                       one_over_range;
            }

        } else if constexpr (RNG::min() > stub) {
            constexpr auto rem = RNG::min() - stub;

            if constexpr (a) {
                return (static_cast<RealType>(
                            g.template operator()<FwdDirection>() - rem) *
                        one_over_range) +
                       a;
            } else {
                return static_cast<RealType>(
                           g.template operator()<FwdDirection>() - rem) *
                       one_over_range;
            }

        } else if constexpr (std::numeric_limits<
                                 typename RNG::result_type>::max() ==
                             RNG::max()) {
            static_assert(RNG::min() == 0);
            static_assert(stub == 1);

            // we treat this separately because there is a tiny chance it might
            // overflow. we convert to float BEFORE adding to stub
            if constexpr (a) {
                return ((static_cast<RealType>(
                             g.template operator()<FwdDirection>()) +
                         static_cast<RealType>(stub)) *
                        one_over_range) +
                       a;
            } else {
                return (static_cast<RealType>(
                            g.template operator()<FwdDirection>()) +
                        static_cast<RealType>(stub)) *
                       one_over_range;
            }
        } else {
            static_assert(RNG::min() == 0);
            static_assert(stub == 1);

            if constexpr (a) {
                return (static_cast<RealType>(
                            g.template operator()<FwdDirection>() + stub) *
                        one_over_range) +
                       a;
            } else {
                return static_cast<RealType>(
                           g.template operator()<FwdDirection>() + stub) *
                       one_over_range;
            }
        }
    }
};

template <class RNG, class RealType = double, bool LeftClosed = true,
          bool RightClosed = false, RealType a = static_cast<RealType>(0.0),
          RealType b = static_cast<RealType>(1.0)>
class uniform_real {
  public:
    using result_type = RealType;

    explicit uniform_real() : rng_() {}

    template <typename... Args>
    explicit uniform_real(Args const &...args) : rng_(args...) {}

    template <typename... Args>
    explicit uniform_real(Args &...args) : rng_(args...) {}

    template <bool FwdDirection = true>
    inline auto operator()() -> result_type {

        static_assert(RNG::max() > RNG::min());
        constexpr auto full_range_int = RNG::max() - RNG::min();

        constexpr auto full_range =
            static_cast<RealType>(full_range_int) +
            static_cast<RealType>(
                !RightClosed) + // right interval might be closed
            static_cast<RealType>(!LeftClosed); // left interval might be open

        constexpr RealType one_over_range = (b - a) / full_range;

        constexpr auto stub = static_cast<RNG::result_type>(!LeftClosed);

        if constexpr (RNG::min() == stub) {
            if constexpr (a) {
                return (static_cast<RealType>(
                            this->rng_.template operator()<FwdDirection>()) *
                        one_over_range) +
                       a;
            } else {
                return static_cast<RealType>(
                           this->rng_.template operator()<FwdDirection>()) *
                       one_over_range;
            }

        } else if constexpr (RNG::min() > stub) {
            constexpr auto rem = RNG::min() - stub;

            if constexpr (a) {
                return (static_cast<RealType>(
                            this->rng_.template operator()<FwdDirection>() -
                            rem) *
                        one_over_range) +
                       a;
            } else {
                return static_cast<RealType>(
                           this->rng_.template operator()<FwdDirection>() -
                           rem) *
                       one_over_range;
            }

        } else if constexpr (std::numeric_limits<
                                 typename RNG::result_type>::max() ==
                             RNG::max()) {
            static_assert(RNG::min() == 0);
            static_assert(stub == 1);

            // we treat this separately because there is a tiny chance it might
            // overflow. we convert to float BEFORE adding to stub
            if constexpr (a) {
                return ((static_cast<RealType>(
                             this->rng_.template operator()<FwdDirection>()) +
                         static_cast<RealType>(stub)) *
                        one_over_range) +
                       a;
            } else {
                return (static_cast<RealType>(
                            this->rng_.template operator()<FwdDirection>()) +
                        static_cast<RealType>(stub)) *
                       one_over_range;
            }

        } else {
            static_assert(RNG::min() == 0);
            static_assert(stub == 1);

            if constexpr (a) {
                return (static_cast<RealType>(
                            this->rng_.template operator()<FwdDirection>() +
                            stub) *
                        one_over_range) +
                       a;
            } else {
                return static_cast<RealType>(
                           this->rng_.template operator()<FwdDirection>() +
                           stub) *
                       one_over_range;
            }
        }
    };

    void discard(unsigned long long z) { rng_.discard(z); }

    auto operator==(uniform_real const &rhs) const noexcept -> bool {
        return this->rng_ == rhs.rng_;
    }

  private:
    RNG rng_{RNG()};
};

template <class RNG, class RealType = double>
using open = uniform_real<RNG, RealType, false>;

template <class RNG, class RealType = double>
using canonical = uniform_real<RNG, RealType>;

template <class RNG, class RealType = double>
using closed = uniform_real<RNG, RealType, true, true>;

} // namespace adhoc

#endif // BRNG_DISTRIBUTION_UNIFORM_DISTRIBUTION
