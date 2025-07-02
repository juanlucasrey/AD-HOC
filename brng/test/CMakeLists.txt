include(../../test_simple/test_simple_cmake.cmake)

file(COPY ${PROJECT_SOURCE_DIR}/brng/test/third_party/randomgen DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILES_MATCHING PATTERN *.txt)
file(COPY ${PROJECT_SOURCE_DIR}/brng/test/third_party/SFMT-src-1.5.1 DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILES_MATCHING PATTERN *.txt)
file(COPY ${PROJECT_SOURCE_DIR}/brng/test/third_party/dSFMT-src-2.2.3 DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILES_MATCHING PATTERN *.txt)
file(COPY ${PROJECT_SOURCE_DIR}/brng/test/third_party/sherif-dear DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILES_MATCHING PATTERN *.csv)

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
    gjrand_engine.cpp
    hc128_engine.cpp
    sfmt_engine.cpp
    dsfmt_engine.cpp
)

run_check("${BRNG}")

set(BRNGStatic
    tools/mask.cpp
)
compile_check("${BRNGStatic}" false)

set(test_name "dynamic.WELL512a")
add_executable(${test_name} WELL512a.cpp third_party/umontreal/WELL512a.c)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include)

set(test_name "dynamic.WELL1024a")
add_executable(${test_name} WELL1024a.cpp third_party/umontreal/WELL1024a.c)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include)

set(test_name "dynamic.WELL19937a")
add_executable(${test_name} WELL19937a.cpp third_party/umontreal/WELL19937a.c)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include)

set(test_name "dynamic.WELL44497a")
add_executable(${test_name} WELL44497a.cpp third_party/umontreal/WELL44497a.c)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include)

set(test_name "dynamic.MRG31k3p")
add_executable(${test_name} MRG31k3p.cpp third_party/umontreal/MRG31k3p.c)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include)

set(test_name "dynamic.MRG32k3a")
add_executable(${test_name} MRG32k3a.cpp third_party/umontreal/MRG32k3a.c)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include)

set(test_name "dynamic.MRG32k5a")
add_executable(${test_name} MRG32k5a.cpp third_party/umontreal/MRG32k5a.c)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include)

set(test_name "dynamic.MRG63k3a")
add_executable(${test_name} MRG63k3a.cpp third_party/umontreal/MRG63k3a.c)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include)

set(test_name "dynamic.PCG_engine")
add_executable(${test_name} PCG_engine.cpp third_party/pcg-c-basic/pcg_basic.c)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include)

add_subdirectory(third_party/PractRand-pre0.95)

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

set(test_name "dynamic.chacha_engine")
add_executable(${test_name} chacha_engine.cpp)
target_link_libraries(${test_name} PRIVATE
    PractRand)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include new/PractRand-pre0.95/include new/PractRand-pre0.95/tools)

set(test_name "dynamic.romu_engine")
add_executable(${test_name} romu_engine.cpp third_party/romu/romu.cpp)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include)

set(test_name "dynamic.sfc_engine")
add_executable(${test_name} sfc_engine.cpp)
target_link_libraries(${test_name} PRIVATE
    PractRand)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include new/PractRand-pre0.95/include new/PractRand-pre0.95/tools)

set(test_name "dynamic.efiix_engine")
add_executable(${test_name} efiix_engine.cpp)
target_link_libraries(${test_name} PRIVATE
    PractRand)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include new/PractRand-pre0.95/include new/PractRand-pre0.95/tools)

set(test_name "dynamic.xsm_engine")
add_executable(${test_name} xsm_engine.cpp)
target_link_libraries(${test_name} PRIVATE
    PractRand)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include new/PractRand-pre0.95/include new/PractRand-pre0.95/tools)

set(test_name "dynamic.isaac_engine")
add_executable(${test_name} isaac_engine.cpp third_party/burtleburtle/isaac.c)
target_link_libraries(${test_name} PRIVATE
    PractRand)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include new/PractRand-pre0.95/include new/PractRand-pre0.95/tools)
