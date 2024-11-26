cmake_minimum_required(VERSION 3.21)
set(CMAKE_CXX_STANDARD 20)

project(adhoc)

# only affects MacOS. MacOSX15.sdk is broken for GCC
if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    set(CMAKE_OSX_SYSROOT /Library/Developer/CommandLineTools/SDKs/MacOSX14.5.sdk)
endif()

enable_testing()
add_subdirectory(test)
add_subdirectory(case_studies)
