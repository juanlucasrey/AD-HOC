# set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD 20)
# add_definitions(-DCODELOGGER)

add_executable(
  adhoc_test
  leafs.t.cpp
  constants.t.cpp
  tuple_utils.t.cpp
  static_tape.t.cpp
  evaluate_fwd.t.cpp
  evaluate_bwd.t.cpp
  tape_size_fwd.t.cpp
  tape_size_bwd.t.cpp
  mcsim.t.cpp
  BS.t.cpp
  yc_test.cpp
  vector.cpp
  # opencl only for clang at the moment
  # opencl.t.cpp
)

# fo coverage
# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-instr-generate -fcoverage-mapping")

# SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wextra -Wconversion -Wdouble-promotion -Wundef -Wdelete-non-virtual-dtor -Wnon-virtual-dtor -Wparentheses -Wcast-align -Wwrite-strings -Wredundant-decls -Wunused-parameter")

if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unsafe-loop-optimizations -Wuseless-cast -Wlogical-op -Wno-psabi")
endif()


target_include_directories(adhoc_test SYSTEM PUBLIC /usr/local/include  ../public)
target_link_directories(adhoc_test PUBLIC /usr/local/lib)

find_package(OpenCL REQUIRED)
target_link_libraries(
  adhoc_test
  gtest
  gtest_main
  "-framework OpenCL"
)

include(GoogleTest)
gtest_discover_tests(adhoc_test)
