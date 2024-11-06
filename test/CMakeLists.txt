function(compile_check filelist fail)
    foreach(testfile IN LISTS filelist)
        file(RELATIVE_PATH target "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/${testfile}")
        string(REPLACE .cpp "" target ${target})
        string(REPLACE / "." target ${target})
        set(test_name "static.${target}")
        # message(STATUS "testfile=${testfile}")
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

# list of static checks that compile
set(StaticWorking
    differential_operator.t.cpp
    differential_operator2.t.cpp
    name.cpp
    combinatorics/multinomial_coefficient_index_sequence.cpp
    combinatorics/multinomial_coefficient.cpp
)

# list of static checks that do not compile (and shouldn't)
set(StaticFailing
    differential_operator.t.fail.cpp)

# list of normal, dynamic tests
set(Dynamic
    black_scholes1.cpp)

compile_check("${StaticWorking}" false)
compile_check("${StaticFailing}" true)
run_check("${Dynamic}")
