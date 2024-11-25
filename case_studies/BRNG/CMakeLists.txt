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

set(BRNG
    mersenne_twister_engine.cpp
)

run_check("${BRNG}")
