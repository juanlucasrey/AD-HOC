#ifndef BRNG_SEED_SEQ
#define BRNG_SEED_SEQ

#include <vector>

namespace adhoc {

template <class UIntType> struct seed_seq final {

    using result_type = UIntType;
    template <class OutputIt> void generate(OutputIt first, OutputIt last) {
        auto input = vals.begin();

        for (; first != last; (void)++first, (void)++input)
            *first = *input;
    }

    std::vector<UIntType> vals{};
};

} // namespace adhoc

#endif // BRNG_SEED_SEQ
