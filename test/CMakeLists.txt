function(compile_check filelist fail)
    foreach(testfile IN LISTS filelist)
        file(RELATIVE_PATH target "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/${testfile}")
        string(REPLACE .cpp "" target ${target})
        string(REPLACE / "." target ${target})
        set(test_name "static.${target}")
        add_executable(${test_name} "${testfile}")
        set_target_properties(${test_name} PROPERTIES EXCLUDE_FROM_ALL true EXCLUDE_FROM_DEFAULT_BUILD true)
        add_test(NAME ${test_name} COMMAND ${CMAKE_COMMAND} --build "${CMAKE_BINARY_DIR}" --target ${test_name})
        if (fail)
            set_tests_properties(${test_name} PROPERTIES WILL_FAIL true)
        endif()
        target_include_directories(${test_name} SYSTEM PUBLIC ../include)
    endforeach()
endfunction()

function(run_check filelist)
    foreach(testfile IN LISTS filelist)
        file(RELATIVE_PATH target "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/${testfile}")
        string(REPLACE .cpp "" target ${target})
        string(REPLACE / "." target ${target})
        set(test_name "dynamic.${target}")
        add_executable(${test_name} "${testfile}")
        add_test(NAME ${test_name} COMMAND ${test_name})
        target_include_directories(${test_name} SYSTEM PUBLIC ../include)
    endforeach()
endfunction()

# for coverage
# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-instr-generate -fcoverage-mapping")
# SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")

if(NOT CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wextra -Wconversion -Wdouble-promotion -Wundef -Wdelete-non-virtual-dtor -Wnon-virtual-dtor -Wparentheses -Wcast-align -Wwrite-strings -Wredundant-decls -Wunused-parameter")
    if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unsafe-loop-optimizations -Wuseless-cast -Wlogical-op")
    endif()
endif()

# if(CMAKE_CXX_COMPILER_ID MATCHES Clang)
#   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconstexpr-steps=1271242 -fconstexpr-depth=2300")
# endif()

# list of static checks that compile
set(StaticWorking
    name.cpp
    dependency.cpp
    combinatorics/multinomial_coefficient_index_sequence.cpp
    combinatorics/multinomial_coefficient.cpp
    utils/index_sequence.cpp
    utils/tuple.cpp
    sort.cpp
    differential_operator.cpp
    half_square.cpp
    backpropagator.cpp
    backpropagator_tools.cpp
    order.cpp
    monomial_included.cpp
)

# list of static checks that compile (only working with GCC at the moment)
set(StaticWorkingGCC
    name_gcc.cpp
)

# list of static checks that do not compile (and shouldn't)
set(StaticFailing
    calc_tree_fail_not_on_tree.cpp
    calc_tree_fail_not_input.cpp
)

# list of normal, dynamic tests
set(Dynamic
    black_scholes1.cpp
    black_scholes2.cpp
    univariate_vs_python.cpp
    utils/univariate.cpp
    univariate_vs_fd.cpp
    calc_tree.cpp
    copy_ordered.cpp
)

# list of normal, dynamic tests (only working with GCC at the moment)
set(DynamicGCC
    black_scholes3.cpp
    univariate_vs_python_gcc.cpp
    univariate_vs_fd_gcc.cpp
)

compile_check("${StaticWorking}" false)
compile_check("${StaticFailing}" true)
if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    compile_check("${StaticWorkingGCC}" false)
endif()

run_check("${Dynamic}")
if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    run_check("${DynamicGCC}")
endif()

file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
