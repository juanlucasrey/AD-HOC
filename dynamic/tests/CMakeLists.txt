# Header-only tests
set(test_name "dynamic.functions.header_only")
add_executable(${test_name} functions.cpp)
target_compile_definitions(${test_name} PRIVATE ADHOC_HEADER_ONLY)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include)

set(test_name "dynamic.checkpoint.header_only")
add_executable(${test_name} checkpoint.cpp)
target_compile_definitions(${test_name} PRIVATE ADHOC_HEADER_ONLY)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include)

set(test_name "dynamic.credit_curve.header_only")
add_executable(${test_name} credit_curve.cpp)
target_compile_definitions(${test_name} PRIVATE ADHOC_HEADER_ONLY)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include)

# Create adhoc library from adhoc.cpp
add_library(adhoc_lib STATIC ../include/adhoc.cpp)
target_include_directories(adhoc_lib PUBLIC ../include)

# Library-linked tests
set(test_name "dynamic.functions.linked")
add_executable(${test_name} functions.cpp)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include)
target_link_libraries(${test_name} PRIVATE adhoc_lib)

set(test_name "dynamic.checkpoint.linked")
add_executable(${test_name} checkpoint.cpp)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include)
target_link_libraries(${test_name} PRIVATE adhoc_lib)

set(test_name "dynamic.credit_curve.linked")
add_executable(${test_name} credit_curve.cpp)
add_test(NAME ${test_name} COMMAND ${test_name})
target_include_directories(${test_name} SYSTEM PUBLIC ../include)
target_link_libraries(${test_name} PRIVATE adhoc_lib)
