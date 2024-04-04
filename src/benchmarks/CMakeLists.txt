
add_executable(setup setup_benchmark.cpp)
add_executable(bench main_benchmark.cpp)
add_executable(misc misc_benchmark.cpp)

target_link_libraries(setup GenomicPIR)
target_link_libraries(bench GenomicPIR)
target_link_libraries(misc GenomicPIR)

set_target_properties(setup PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(bench PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(misc PROPERTIES EXCLUDE_FROM_ALL TRUE)