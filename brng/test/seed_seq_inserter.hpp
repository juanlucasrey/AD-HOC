#ifndef BRNG_SEED_SEQ_INSERTER
#define BRNG_SEED_SEQ_INSERTER

#include <tools/mask.hpp>

#include <bit>
#include <cstdint>
#include <vector>

namespace adhoc {

class seed_seq_inserter final {
  public:
    using result_type = std::uint32_t;

    template <class UIntType>
    explicit seed_seq_inserter(std::vector<UIntType> const &to_insert) {
        constexpr std::size_t w = std::numeric_limits<UIntType>::digits;
        static_assert(std::has_single_bit(w));

        if constexpr (w == 32) {
            vals = to_insert;
        } else if constexpr (w < 32) {
            vals.resize((to_insert.size() * w + (32 - 1)) / 32);

            std::size_t idx{0};
            std::size_t bits_used{0};
            for (std::size_t i = 0; i < to_insert.size(); ++i) {
                vals[idx] |= static_cast<std::uint32_t>(to_insert[i])
                             << bits_used;
                bits_used += w;
                if (bits_used == 32) {
                    bits_used = 0;
                    ++idx;
                }
            }
        } else {
            constexpr std::size_t ratio = w / 32;
            vals.resize(ratio * to_insert.size());

            std::size_t idx{0};
            for (std::size_t i = 0; i < to_insert.size(); ++i) {
                for (std::size_t j = 0; j < ratio; ++j) {
                    vals[idx++] =
                        static_cast<std::uint32_t>(to_insert[i] >> (32 * j));
                }
            }
        }
    }

    template <class OutputIt> void generate(OutputIt first, OutputIt last) {
        auto input = vals.begin();

        for (; first != last; (void)++first, (void)++input) {
            *first = *input;
        }
    }

  private:
    std::vector<std::uint32_t> vals;
};

} // namespace adhoc

#endif // BRNG_SEED_SEQ_INSERTER
