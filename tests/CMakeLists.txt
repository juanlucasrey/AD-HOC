set(CMAKE_CXX_STANDARD 17)

add_executable(
  adhoc_test
  test.cpp
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

target_include_directories(adhoc_test SYSTEM PUBLIC /usr/local/include  ../public)
target_link_directories(adhoc_test SYSTEM PUBLIC /usr/local/lib)

target_link_libraries(
  adhoc_test
  gtest
  gtest_main
)

include(GoogleTest)
gtest_discover_tests(adhoc_test)
