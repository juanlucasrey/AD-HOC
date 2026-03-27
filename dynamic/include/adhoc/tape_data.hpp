#ifndef ADHOC_TAPE_DATA_HPP
#define ADHOC_TAPE_DATA_HPP

#include <cstdint>
#include <vector>

namespace adhoc {

enum class OpCode : std::uint8_t {
    ADD,
    SUB,
    MUL,
    ADD_C,
    SUB_C,
    MUL_C,
    NORM, // a.k.a. SQUARE
    ABS,
    INV,
    SQRT,
    POW_C,
    EXP,
    LOG,
    ERF,
    ERFC,
    COS,
    REG_INPUT,
    REG_OUTPUT,
    // IFT
    // IFT_START,
    // IFT_END,
    // REG_INPUT_IFT,
    // REG_OUTPUT_IFT,
};

struct TapeData {
    std::vector<OpCode> ops;
    std::vector<double> vals;
    std::vector<std::size_t> ids;
    size_t next_id{ 0 };
};

} // namespace adhoc

#endif // ADHOC_TAPE_DATA_HPP
