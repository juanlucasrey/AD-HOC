include(../../test_simple/test_simple_cmake.cmake)

set(BRNG
    tools/circular_buffer.cpp
    mersenne_twister_engine.cpp
    linear_congruential_engine.cpp
    subtract_with_carry_engine.cpp
    philox_engine.cpp
    sobol_engine.cpp
)

run_check("${BRNG}")

set(BRNGStatic
    tools/mask.cpp
    tools/modular_arithmetic.cpp
)
compile_check("${BRNGStatic}" false)
