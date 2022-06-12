cmake_minimum_required(VERSION 3.16)
project(Atif CXX)

set(CMAKE_C_STANDARD 1)
set(CMAKE_CXX_STANDARD 17)

FILE(GLOB MyCppSources ../source/*.cpp)

add_executable(AtifExe ${MyCppSources})

target_include_directories(AtifExe PUBLIC ../source)
