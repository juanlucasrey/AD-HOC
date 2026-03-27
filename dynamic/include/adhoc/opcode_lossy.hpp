#ifndef ADHOC_OPCODE_LOSSY_HPP
#define ADHOC_OPCODE_LOSSY_HPP

#include <cstdint>

namespace adhoc {

// Simplified opcodes for lossy tape - only additions and multiplications
enum class LossyOpCode : std::uint8_t {
    ADD_FROM_SEED,
    COPY,          // result = source
    COPY_MINUS,    // result = -source
    ADD,           // result += source
    SUB,           // result -= source
    MINUS_INPLACE, // result = -result
    MUL_INPLACE,   // result *= factor
    MUL_ADD,       // result += factor * source (multiply and accumulate)
    MUL_SET,       // result = factor * source (multiply and set)
};

} // namespace adhoc

#endif // ADHOC_OPCODE_LOSSY_HPP
