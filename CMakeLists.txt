# minimum cmake requirement
cmake_minimum_required(VERSION 3.16)

# project name
project(Acorn)

# set output path, C++ standard, and compile commands export
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)

# set debug flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg -Og -Wall -Wextra -flarge-source-files -fprofile-arcs -ftest-coverage")

# set release flags
set(CMAKE_CXX_FLAGS_RELEASE "-mavx -s -O3 -Wall -Wextra -flarge-source-files")

# include directories
include_directories(include lib)

# find the necessary packages
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)

add_subdirectory(program/integral)
