include(../../test_simple/test_simple_cmake.cmake)

file(COPY ${PROJECT_SOURCE_DIR}/brng/test/external/randomgen DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILES_MATCHING PATTERN *.txt)

set(BRNG
    tools/circular_buffer.cpp
    tools/unshift.cpp
    tools/modular_arithmetic.cpp
    distribution/centered_canonical_distribution.cpp
    distribution/normal_distribution.cpp
    linear_congruential_engine.cpp
    subtract_with_carry_engine.cpp
    philox_engine.cpp
    sobol_engine.cpp
    lfsr_xor_engine.cpp
    arc4_engine.cpp
    splitmix_engine.cpp
    chacha_engine.cpp
    gjrand_engine.cpp
    hc128_engine.cpp
)

run_check("${BRNG}")

set(BRNGStatic
    tools/mask.cpp
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

set(test_name "dynamic.MRG31k3p")
add_executable(${test_name} MRG31k3p.cpp external/MRG31k3p.c)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include)

set(test_name "dynamic.MRG32k3a")
add_executable(${test_name} MRG32k3a.cpp external/MRG32k3a.c)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include)

set(test_name "dynamic.MRG32k5a")
add_executable(${test_name} MRG32k5a.cpp external/MRG32k5a.c)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include)

set(test_name "dynamic.MRG63k3a")
add_executable(${test_name} MRG63k3a.cpp external/MRG63k3a.c)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include)

set(test_name "dynamic.PCG_engine")
add_executable(${test_name} PCG_engine.cpp external/pcg_basic.c)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include)

add_subdirectory(external/PractRand-pre0.95)

set(test_name "dynamic.mersenne_twister_engine")
add_executable(${test_name} mersenne_twister_engine.cpp)
target_link_libraries(${test_name} PRIVATE
    PractRand)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include new/PractRand-pre0.95/include new/PractRand-pre0.95/tools)

set(test_name "dynamic.jsf_engine")
add_executable(${test_name} jsf_engine.cpp)
target_link_libraries(${test_name} PRIVATE
    PractRand)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include new/PractRand-pre0.95/include new/PractRand-pre0.95/tools)

set(test_name "dynamic.hc256_engine")
add_executable(${test_name} hc256_engine.cpp)
target_link_libraries(${test_name} PRIVATE
    PractRand)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include new/PractRand-pre0.95/include new/PractRand-pre0.95/tools)

