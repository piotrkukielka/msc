cmake_minimum_required(VERSION 3.15)
project(msc_0505)

set(CMAKE_CXX_STANDARD 17)

add_executable(msc_0505 main.cpp)
target_sources(msc_0505 PUBLIC src/solvers.cpp)
target_include_directories(msc_0505 PUBLIC headers)
