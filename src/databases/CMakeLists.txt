# Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.10)

# Project name
project(Databases)

# Set the include directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Add the source files
set(SOURCES
    aes_plaintext_database.cpp
    plaintext_database.cpp
    FHE_SIMD_database.cpp
    FHE_SIMD_counting_query.cpp
    FHE_SIMD_ctxt_packing.cpp
    FHE_SIMD_DML.cpp
    FHE_SIMD_MAF_query.cpp
    FHE_SIMD_PRS_query.cpp
    FHE_SIMD_range_query.cpp
    FHE_SIMD_setup.cpp
    FHE_SIMD_similarity_query.cpp
    FHE_disk_database.cpp
    tools.cpp
    tools.hpp
    comparator.cpp
    comparator.hpp
    ../globals.hpp
)

# Create the library
add_library(Databases ${SOURCES})
target_link_libraries(Databases helib)
