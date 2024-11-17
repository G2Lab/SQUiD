find_package(benchmark REQUIRED)
find_package(OpenSSL REQUIRED)

# Add include directories for OpenSSL
include_directories(${OPENSSL_INCLUDE_DIR})

add_subdirectory(databases)

add_library(SQUiD
    globals.hpp
)

target_link_libraries(SQUiD helib)
target_link_libraries(SQUiD benchmark::benchmark)
target_link_libraries(SQUiD Databases)
target_link_libraries(SQUiD ${OPENSSL_LIBRARIES})

#Add JSON package
include(FetchContent)
FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz)
FetchContent_MakeAvailable(json)
target_link_libraries(SQUiD nlohmann_json::nlohmann_json)

#Add CURL library
find_package(CURL REQUIRED)
target_link_libraries(SQUiD CURL::libcurl)


option(BUILD_SQUID_API "Build the API for SQUiD" OFF)
if(BUILD_SQUID_API)
    add_subdirectory(API)
endif()


add_executable(main main.cpp)
add_executable(squid SQUiD.cpp)
add_executable(data_owner_encryptor data_owner_encryptor.cpp)

target_link_libraries(main SQUiD)
target_link_libraries(squid SQUiD)
target_link_libraries(data_owner_encryptor SQUiD)

enable_testing()
add_subdirectory(tests)