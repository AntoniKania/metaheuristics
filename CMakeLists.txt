cmake_minimum_required(VERSION 3.26)
project(metaheuristics)

set(CMAKE_CXX_STANDARD 17)
find_package(TBB REQUIRED)

add_executable(metaheuristics main.cpp)
target_link_libraries(metaheuristics PRIVATE TBB::tbb)