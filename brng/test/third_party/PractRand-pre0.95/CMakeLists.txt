
set(PractRandFiles
    src/math.cpp
    src/non_uniform.cpp
    src/platform_specifics.cpp
    src/rand.cpp
    src/sha2.cpp
    src/test_batteries.cpp
    src/tests.cpp
    src/RNGs/arbee.cpp
    src/RNGs/chacha.cpp
    src/RNGs/efiix.cpp
    src/RNGs/hc256.cpp
    src/RNGs/isaac.cpp
    src/RNGs/jsf.cpp
    src/RNGs/mt19937.cpp
    src/RNGs/rarns.cpp
    src/RNGs/salsa.cpp
    src/RNGs/sfc.cpp
    src/RNGs/sha2_based_pool.cpp
    src/RNGs/trivium.cpp
    src/RNGs/xsm.cpp
    src/RNGs/other/fibonacci.cpp
    src/RNGs/other/indirection.cpp
    src/RNGs/other/mult.cpp
    src/RNGs/other/rng_sets.cpp
    src/RNGs/other/simple.cpp
    src/RNGs/other/transform.cpp
)

# disable warnings for these external files
set_source_files_properties(
  ${PractRandFiles}
  PROPERTIES
  COMPILE_FLAGS "-w"
)

add_library(
    PractRand
    "${PractRandFiles}"
)

target_include_directories(PractRand SYSTEM PUBLIC include)
