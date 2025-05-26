include(../../test_simple/test_simple_cmake.cmake)

set(BRNG
    mersenne_twister_engine.cpp
)

run_check("${BRNG}")
