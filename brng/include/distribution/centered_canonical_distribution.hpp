#ifndef BRNG_DISTRIBUTION_CENTERED_CANONICAL_DISTRIBUTION
#define BRNG_DISTRIBUTION_CENTERED_CANONICAL_DISTRIBUTION

#include <limits>

namespace adhoc {

template <class RealType = double, bool LeftClosed = true,
          bool RightClosed = false>
class centered_canonical_distribution final {
  public:
    using result_type = RealType;

    template <bool FwdDirection = true, typename RNG>
    auto operator()(RNG &g) -> result_type {
        static_assert(RNG::max() > RNG::min());
        constexpr auto full_range_int = RNG::max() - RNG::min();
        constexpr auto cutoff =
            (full_range_int & 1) + (full_range_int / 2) + RNG::min();

        constexpr auto full_range =
            static_cast<RealType>(full_range_int) +
            static_cast<RealType>(
                !RightClosed) + // right interval might be closed
            static_cast<RealType>(!LeftClosed); // left interval might be open

        constexpr RealType one_over_range = 1.0 / full_range;

        constexpr auto stub_left = static_cast<RNG::result_type>(!LeftClosed);
        constexpr auto stub_right = static_cast<RNG::result_type>(!RightClosed);

        auto const val = g.template operator()<FwdDirection>();
        if (val <= cutoff) {
            if constexpr (RNG::min() == stub_left) {
                return static_cast<RealType>(val) * one_over_range;
            } else if constexpr (RNG::min() > stub_left) {
                constexpr auto rem = RNG::min() - stub_left;
                return static_cast<RealType>(val - rem) * one_over_range;
            } else {
                static_assert(RNG::min() == 0);
                static_assert(stub_left == 1);
                return static_cast<RealType>(val + stub_left) * one_over_range;
            }
        } else {
            if constexpr (stub_right) {
                return -static_cast<RealType>(RNG::max() - val + stub_right) *
                       one_over_range;
            } else {
                return -static_cast<RealType>(RNG::max() - val) *
                       one_over_range;
            }
        }
    }
};

// why do we use (-0.5, -0) union (0, 0.5) instead of (0, 1)? because in IEEE
// 754, subnormal numbers near 0 allow more precision. in this representation,
// we allow more precision in the ", 1)" border, now represented by ", -0)""
template <class RNG, class RealType = double, bool LeftClosed = true,
          bool RightClosed = false>
class centered_canonical {
  public:
    using result_type = RealType;

    explicit centered_canonical() : rng_() {}

    template <typename... Args>
    explicit centered_canonical(Args const &...args) : rng_(args...) {}

    template <typename... Args>
    explicit centered_canonical(Args &...args) : rng_(args...) {}

    template <bool FwdDirection = true>
    inline auto operator()() -> result_type {

        static_assert(RNG::max() > RNG::min());
        constexpr auto full_range_int = RNG::max() - RNG::min();
        constexpr auto cutoff =
            (full_range_int & 1) + (full_range_int / 2) + RNG::min();

        // we use the highest possible float precision. why? because this might
        // be infinity if we use float
        constexpr auto full_range =
            static_cast<double>(full_range_int) +
            static_cast<double>(
                !RightClosed) +               // right interval might be closed
            static_cast<double>(!LeftClosed); // left interval might be open

        // we (potentially) downgrade
        constexpr auto one_over_range = static_cast<RealType>(1.0 / full_range);

        constexpr auto stub_left = static_cast<RNG::result_type>(!LeftClosed);
        constexpr auto stub_right = static_cast<RNG::result_type>(!RightClosed);

        auto const val = this->rng_.template operator()<FwdDirection>();
        if (val <= cutoff) {
            if constexpr (RNG::min() == stub_left) {
                return static_cast<RealType>(val) * one_over_range;
            } else if constexpr (RNG::min() > stub_left) {
                constexpr auto rem = RNG::min() - stub_left;
                return static_cast<RealType>(val - rem) * one_over_range;
            } else {
                static_assert(RNG::min() == 0);
                static_assert(stub_left == 1);
                return static_cast<RealType>(val + stub_left) * one_over_range;
            }
        } else {
            if constexpr (stub_right) {
                return -static_cast<RealType>(RNG::max() - val + stub_right) *
                       one_over_range;
            } else {
                return -static_cast<RealType>(RNG::max() - val) *
                       one_over_range;
            }
        }
    };

    void discard(unsigned long long z) { rng_.discard(z); }

    auto operator==(centered_canonical const &rhs) const noexcept -> bool {
        return this->rng_ == rhs.rng_;
    }

  private:
    RNG rng_{RNG()};
};

template <class RNG, class RealType = double>
using centered_canonical_open = centered_canonical<RNG, RealType, false>;

template <class RNG, class RealType = double>
using centered_canonical_default = centered_canonical<RNG, RealType>;

template <class RNG, class RealType = double>
using centered_canonical_closed = centered_canonical<RNG, RealType, true, true>;

} // namespace adhoc

#endif // BRNG_DISTRIBUTION_CENTERED_CANONICAL_DISTRIBUTION
