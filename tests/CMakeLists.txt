set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD 20)
# add_definitions(-DCODELOGGER)

add_executable(
  adhoc_test
  new/leafs.t.cpp
  new/static_tape.t.cpp
  new/tape_size_fwd.t.cpp
  BS.t.cpp
  adhoc.t.cpp
  constants.t.cpp
  static_tape.t.cpp
  tape_size_bwd.t.cpp
  tape_size_fwd.t.cpp
  evaluate_fwd.t.cpp
  mcsim.t.cpp
  # opencl only for clang at the moment
  # opencl.t.cpp
)

if (NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # googletest should be installed via brew, but it only works with clang compiled code
    # To ensure it works with other compilers (mostly GCC), I need to download and statically link.
    include(FetchContent)
    FetchContent_Declare(
      googletest
      GIT_REPOSITORY https://github.com/google/googletest.git
      GIT_TAG main
    )

    FetchContent_MakeAvailable(googletest)
endif()

# SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wextra -Wconversion -Wdouble-promotion -Wundef -Wdelete-non-virtual-dtor -Wnon-virtual-dtor -Wparentheses -Wcast-align -Wwrite-strings -Wredundant-decls -Wunused-parameter")

if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unsafe-loop-optimizations -Wuseless-cast -Wlogical-op")
endif()


target_include_directories(adhoc_test SYSTEM PUBLIC /usr/local/include  ../public)
target_link_directories(adhoc_test SYSTEM PUBLIC /usr/local/lib)

find_package(OpenCL REQUIRED)
target_link_libraries(
  adhoc_test
  gtest
  gtest_main
  "-framework OpenCL"
)

include(GoogleTest)
gtest_discover_tests(adhoc_test)
