project(radial_diff)

cmake_minimum_required (VERSION 3.5)

add_definitions(-std=c++11)

option (BUILD_TESTS "Enable to build tests" ON)
option (BUILD_EXAMPLES "Enable to build examples" ON)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

set(sources main.cpp)

find_package ( Boost )

include_directories ( SYSTEM ${Boost_INCLUDE_DIRS} )

include_directories(third-party)
include_directories(include)

add_executable(RadialDiff ${sources})

#if(BUILD_TESTS)
    add_subdirectory(third-party/googletest-release-1.8.1)
    add_subdirectory(tests)
#endif()

    add_subdirectory(examples)
