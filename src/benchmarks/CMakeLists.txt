add_executable(bench main_benchmark.cpp)
add_executable(misc misc_benchmark.cpp)

target_link_libraries(bench SQUiD)
target_link_libraries(misc SQUiD)

set_target_properties(bench PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(misc PROPERTIES EXCLUDE_FROM_ALL TRUE)