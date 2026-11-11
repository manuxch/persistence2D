cmake_minimum_required(VERSION 3.20)
project(persistence2D VERSION 0.2.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD          20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS        OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# GUDHI (header-only) and the Boost headers it needs. Installed GUDHI exports GUDHIConfig.cmake; a source tree
# can be given with -DGUDHI_DIR=<dir containing GUDHIConfig.cmake> or -DGUDHI_INCLUDE_DIRS=<include dir>.
if(NOT GUDHI_INCLUDE_DIRS)
    find_package(GUDHI 3.10 REQUIRED)
endif()
find_package(Boost 1.71 REQUIRED CONFIG)
find_package(Threads REQUIRED)

add_library(p2d_core STATIC
    src/grid.cpp
    src/distance.cpp
    src/persistence.cpp
    src/continuum.cpp
    src/height_exact.cpp
    src/descriptors.cpp
    src/params.cpp
    src/output.cpp
)
target_include_directories(p2d_core PUBLIC include)
target_include_directories(p2d_core SYSTEM PRIVATE ${GUDHI_INCLUDE_DIRS})
target_link_libraries(p2d_core PUBLIC Threads::Threads PRIVATE Boost::headers)
target_compile_definitions(p2d_core PUBLIC P2D_VERSION="${PROJECT_VERSION}")

add_executable(persistence2D src/main.cpp)
target_link_libraries(persistence2D PRIVATE p2d_core)

add_executable(test_persistence2D tests/test_persistence2D.cpp)
target_link_libraries(test_persistence2D PRIVATE p2d_core)

# -march=native tunes for the build CPU; turn it off to build on a cluster login node for other nodes.
option(P2D_NATIVE "Compile with -march=native in Release" ON)
# A fully static executable runs on machines without GUDHI, Boost or a recent libstdc++ (e.g. a cluster).
option(P2D_STATIC "Link persistence2D statically" OFF)

foreach(tgt p2d_core persistence2D test_persistence2D)
    target_compile_options(${tgt} PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang>:
            -Wall -Wextra -Wpedantic
            $<$<CONFIG:Release>:-O3>
            $<$<AND:$<CONFIG:Release>,$<BOOL:${P2D_NATIVE}>>:-march=native>
            $<$<CONFIG:Debug>:-g -fsanitize=address,undefined>
        >)
    target_link_options(${tgt} PRIVATE
        $<$<AND:$<CXX_COMPILER_ID:GNU,Clang>,$<CONFIG:Debug>>:-fsanitize=address,undefined>)
endforeach()
if(P2D_STATIC)
    target_link_options(persistence2D PRIVATE -static)
endif()

enable_testing()
add_test(NAME persistence2D_unit COMMAND test_persistence2D)
find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND)
    add_test(NAME persistence2D_cli
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_cli.py
                     $<TARGET_FILE:persistence2D> ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
endif()

message(STATUS "persistence2D : ${PROJECT_VERSION}")
message(STATUS "GUDHI         : ${GUDHI_VERSION} (${GUDHI_INCLUDE_DIRS})")
message(STATUS "Build type    : ${CMAKE_BUILD_TYPE}")
