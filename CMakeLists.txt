cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(PIR)

add_definitions(-DSOURCE_DIR="${CMAKE_SOURCE_DIR}")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-p -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# STRINGS avoids having the 2 newline characters at the end of the string.
# Alternatively it's possible to use file(READ ...) and then
# string(REGEX REPLACE "\n$" "" HELIB_VERSION "${HELIB_VERSION}")
file(STRINGS "./HElibPublicKeySwitch/VERSION" HELIB_VERSION)

find_package(helib "${HELIB_VERSION}" EXACT REQUIRED)

if(helib_FOUND)
    message(STATUS "HElib found at: ${helib_DIR}")
else()
    message(FATAL_ERROR "HElib not found!")
endif()


add_subdirectory(src)
add_subdirectory(src/benchmarks)
