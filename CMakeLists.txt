cmake_minimum_required(VERSION 3.25)
project(cpp_language_feature)

set(CMAKE_CXX_STANDARD 14)

add_executable(if_test if_test.c)
add_executable(gep_test gep_test.c)
add_executable(cast_test cast_test.c)
add_executable(const_test const_test.c)
add_executable(function_test function_test.c)
