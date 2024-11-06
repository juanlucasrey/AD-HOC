# set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD 20)
# set(CMAKE_CXX_STANDARD 23)

add_executable(
  adhoc_test4
  finite_differences.cpp
  utils/univariate.t.cpp
  utils/index_sequence.t.cpp
  monomial_included.t.cpp
  differential_operator.t.cpp
  univariate_vs_fd.t.cpp
  univariate_vs_python.t.cpp
  order.t.cpp
  calc_tree.t.cpp
  sort.t.cpp
  utils/tuple.t.cpp
  backpropagator.t.cpp
  black_scholes1.t.cpp
  black_scholes2.t.cpp
  # black_scholes3.t.cpp
  backpropagator_tools.t.cpp
)

# fo coverage
# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-instr-generate -fcoverage-mapping")

# SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wextra -Wconversion -Wdouble-promotion -Wundef -Wdelete-non-virtual-dtor -Wnon-virtual-dtor -Wparentheses -Wcast-align -Wwrite-strings -Wredundant-decls -Wunused-parameter")

# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconstexpr-steps=1271242 -fconstexpr-depth=3000")

if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unsafe-loop-optimizations -Wuseless-cast -Wlogical-op")
endif()

# if(CMAKE_CXX_COMPILER_ID MATCHES Clang)
#   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconstexpr-steps=1271242 -fconstexpr-depth=2300")
# endif()

target_include_directories(adhoc_test4 SYSTEM PUBLIC /usr/local/include /usr/local/include/eigen3 ../include)
# target_include_directories(adhoc_test4 SYSTEM PUBLIC /usr/local/include /usr/local/include/eigen3 ../public4 ../../boost_1_85_0)
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
