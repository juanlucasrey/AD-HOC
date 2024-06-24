# set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD 20)

add_executable(
  adhoc_test4
  adhoc.t.cpp
  differential_operator.t.cpp
  univariate_vs_fd.t.cpp
  univariate_vs_python.t.cpp
  saturation_arithmetic.t.cpp
  name.t.cpp
)

# fo coverage
# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-instr-generate -fcoverage-mapping")

# SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wextra -Wconversion -Wdouble-promotion -Wundef -Wdelete-non-virtual-dtor -Wnon-virtual-dtor -Wparentheses -Wcast-align -Wwrite-strings -Wredundant-decls -Wunused-parameter")

# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconstexpr-steps=1271242 -fconstexpr-depth=3000")

if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unsafe-loop-optimizations -Wuseless-cast -Wlogical-op")
endif()


target_include_directories(adhoc_test4 SYSTEM PUBLIC /usr/local/include /usr/local/include/eigen3 ../public4)
target_link_directories(adhoc_test4 PUBLIC /usr/local/lib)

find_package(OpenCL REQUIRED)
target_link_libraries(
  adhoc_test4
  gtest
  gtest_main
  "-framework OpenCL"
)

file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

include(GoogleTest)
gtest_discover_tests(adhoc_test4)
