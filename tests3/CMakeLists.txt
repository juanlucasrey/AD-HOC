# set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD 20)
# add_definitions(-DCODELOGGER)

# find_package(Eigen3 REQUIRED)

if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    # without this parameter loops will not be vectorised
    # set_source_files_properties(random/tools/scrambler.cpp PROPERTIES COMPILE_FLAGS "$<$<CONFIG:Release,RelWithDebInfo>:--param vect-max-version-for-alias-checks=100>")
    # set_source_files_properties(random/tools/scrambler.cpp PROPERTIES COMPILE_FLAGS "$<$<CONFIG:Release,RelWithDebInfo>:-fopt-info-vec-all --param vect-max-version-for-alias-checks=100>")
    # set_source_files_properties(least_squares.t.cpp PROPERTIES COMPILE_FLAGS "-Wno-useless-cast")
endif()

add_executable(
  adhoc_test3
  dependency.t.cpp
  ders.t.cpp
  raw.t.cpp
  differential_operator.t.cpp
  differential_operator_2.t.cpp
  tape.t.cpp
  tape2.t.cpp
  calc_tree.t.cpp
  partition/binomial_coefficient.t.cpp
)

# fo coverage
# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-instr-generate -fcoverage-mapping")

# SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wextra -Wconversion -Wdouble-promotion -Wundef -Wdelete-non-virtual-dtor -Wnon-virtual-dtor -Wparentheses -Wcast-align -Wwrite-strings -Wredundant-decls -Wunused-parameter")

if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unsafe-loop-optimizations -Wuseless-cast -Wlogical-op")
endif()


target_include_directories(adhoc_test3 SYSTEM PUBLIC /usr/local/include /usr/local/include/eigen3 ../public3)
target_link_directories(adhoc_test3 PUBLIC /usr/local/lib)

find_package(OpenCL REQUIRED)
target_link_libraries(
  adhoc_test3
  gtest
  gtest_main
  "-framework OpenCL"
)

include(GoogleTest)
gtest_discover_tests(adhoc_test3)
