include(../../test_simple/test_simple_cmake.cmake)

set(BRNG
    tools/circular_buffer.cpp
    tools/bit_find.cpp
    tools/unshift.cpp
    mersenne_twister_engine.cpp
    linear_congruential_engine.cpp
    subtract_with_carry_engine.cpp
    philox_engine.cpp
    sobol_engine.cpp
    lfsr_xor_engine.cpp
)

run_check("${BRNG}")

set(BRNGStatic
    tools/mask.cpp
    tools/modular_arithmetic.cpp
)
compile_check("${BRNGStatic}" false)

set(test_name "dynamic.WELL512a")
add_executable(${test_name} WELL512a.cpp external/WELL512a.c)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include)

set(test_name "dynamic.WELL1024a")
add_executable(${test_name} WELL1024a.cpp external/WELL1024a.c)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include)

set(test_name "dynamic.WELL19937a")
add_executable(${test_name} WELL19937a.cpp external/WELL19937a.c)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include)

set(test_name "dynamic.WELL44497a")
add_executable(${test_name} WELL44497a.cpp external/WELL44497a.c)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include)
