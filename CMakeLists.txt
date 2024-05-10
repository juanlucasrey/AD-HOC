cmake_minimum_required(VERSION 3.21)
set(CMAKE_CXX_STANDARD 17)

project(adhoc)

enable_testing()

# googletest should be installed via brew, but it only works with clang compiled code
# To ensure it works with other compilers (mostly GCC), I need to download and statically link.
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG main
)

FetchContent_MakeAvailable(googletest)

add_subdirectory(tests)
add_subdirectory(tests3)
add_subdirectory(tests4)
