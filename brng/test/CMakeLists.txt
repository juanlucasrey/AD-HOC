include(../../test_simple/test_simple_cmake.cmake)

set(BRNG
    tools/circular_buffer.cpp
    mersenne_twister_engine.cpp
)

run_check("${BRNG}")

set(BRNGStatic
    tools/mask.cpp
)
compile_check("${BRNGStatic}" false)
