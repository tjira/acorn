# minimum cmake requirement
cmake_minimum_required(VERSION 3.16)

# project name
project(Acorn)

# include the fetch plugin
include(FetchContent)

# set the c++ standard
set(CMAKE_CXX_STANDARD 20)

# set policies and compile commands
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# set output path
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

# set the base compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg -Og -Wall -Wextra --coverage")
set(CMAKE_CXX_FLAGS_RELEASE "-s -O3 -flto=auto -march=native")

# include directories
include_directories(example/diagram external/include external/include/torch/csrc/api/include include lib)

# find the necessary packages
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)

# find installed libraries
find_library(LIBTORCH_CPU NAMES torch_cpu PATHS external/lib)
find_library(LIBTORCH NAMES torch PATHS external/lib)
find_library(LIBFFTW NAMES fftw3 PATHS external/lib)
find_library(LIBINT NAMES int2 PATHS external/lib)
find_library(LIBC10 NAMES c10 PATHS external/lib)

# add the base library
add_library(acorn_base STATIC src/expression.cpp src/linalg.cpp src/timer.cpp)

# link the base library
target_link_libraries(acorn_base Eigen3::Eigen)

# add subprograms
add_subdirectory(program/expression)
add_subdirectory(program/transform)
add_subdirectory(program/integral)
add_subdirectory(program/cdyn)
add_subdirectory(program/qdyn)
add_subdirectory(program/ci)
add_subdirectory(program/hf)
add_subdirectory(program/mp)

# enable testing
enable_testing()

# add general tests
add_test(NAME hf                   COMMAND ${CMAKE_MAKE_PROGRAM} hf                   WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME mp2                  COMMAND ${CMAKE_MAKE_PROGRAM} mp2                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME mp3                  COMMAND ${CMAKE_MAKE_PROGRAM} mp3                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME fci                  COMMAND ${CMAKE_MAKE_PROGRAM} fci                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME qdyn_1d_HO_imaginary COMMAND ${CMAKE_MAKE_PROGRAM} qdyn_1d_HO_imaginary WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME qdyn_1d_HO_real      COMMAND ${CMAKE_MAKE_PROGRAM} qdyn_1d_HO_real      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME qdyn_2d_HO_imaginary COMMAND ${CMAKE_MAKE_PROGRAM} qdyn_2d_HO_imaginary WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME qdyn_2d_HO_real      COMMAND ${CMAKE_MAKE_PROGRAM} qdyn_2d_HO_real      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME qdyn_3d_HO_imaginary COMMAND ${CMAKE_MAKE_PROGRAM} qdyn_3d_HO_imaginary WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME qdyn_3d_HO_real      COMMAND ${CMAKE_MAKE_PROGRAM} qdyn_3d_HO_real      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)

# add expected results of general tests
set_property(TEST hf                   PROPERTY PASS_REGULAR_EXPRESSION "FINAL SINGLE POINT ENERGY: -74.965901")
set_property(TEST mp2                  PROPERTY PASS_REGULAR_EXPRESSION "FINAL SINGLE POINT ENERGY: -75.004854")
set_property(TEST mp3                  PROPERTY PASS_REGULAR_EXPRESSION "FINAL SINGLE POINT ENERGY: -75.015455")
set_property(TEST fci                  PROPERTY PASS_REGULAR_EXPRESSION "FINAL SINGLE POINT ENERGY: -75.020410")
set_property(TEST qdyn_1d_HO_imaginary PROPERTY PASS_REGULAR_EXPRESSION "FINAL WFN 02 ENERGY: 2.500001"        )
set_property(TEST qdyn_1d_HO_real      PROPERTY PASS_REGULAR_EXPRESSION "FINAL WFN 00 ENERGY: 1.124931"        )
set_property(TEST qdyn_2d_HO_imaginary PROPERTY PASS_REGULAR_EXPRESSION "FINAL WFN 02 ENERGY: 2.000001"        )
set_property(TEST qdyn_2d_HO_real      PROPERTY PASS_REGULAR_EXPRESSION "FINAL WFN 00 ENERGY: 2.249863"        )
set_property(TEST qdyn_3d_HO_imaginary PROPERTY PASS_REGULAR_EXPRESSION "FINAL WFN 02 ENERGY: 2.500001"        )
set_property(TEST qdyn_3d_HO_real      PROPERTY PASS_REGULAR_EXPRESSION "FINAL WFN 00 ENERGY: 3.374795"        )

# add random tests
add_test(NAME qdyn_1d_tully_1         COMMAND ${CMAKE_MAKE_PROGRAM} qdyn_1d_tully_1         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME cdyn_lz-dia_1d_tully_1  COMMAND ${CMAKE_MAKE_PROGRAM} cdyn_lz-dia_1d_tully_1  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME cdyn_lz-adia_1d_tully_1 COMMAND ${CMAKE_MAKE_PROGRAM} cdyn_lz-adia_1d_tully_1 WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME qdyn_1d_tully_2         COMMAND ${CMAKE_MAKE_PROGRAM} qdyn_1d_tully_2         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME cdyn_lz-dia_1d_tully_2  COMMAND ${CMAKE_MAKE_PROGRAM} cdyn_lz-dia_1d_tully_2  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME cdyn_lz-adia_1d_tully_2 COMMAND ${CMAKE_MAKE_PROGRAM} cdyn_lz-adia_1d_tully_2 WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME qdyn_1d_ds_1            COMMAND ${CMAKE_MAKE_PROGRAM} qdyn_1d_ds_1            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME cdyn_lz-dia_1d_ds_1     COMMAND ${CMAKE_MAKE_PROGRAM} cdyn_lz-dia_1d_ds_1     WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME cdyn_lz-adia_1d_ds_1    COMMAND ${CMAKE_MAKE_PROGRAM} cdyn_lz-adia_1d_ds_1    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME qdyn_1d_ts_1            COMMAND ${CMAKE_MAKE_PROGRAM} qdyn_1d_ts_1            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME cdyn_lz-dia_1d_ts_1     COMMAND ${CMAKE_MAKE_PROGRAM} cdyn_lz-dia_1d_ts_1     WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME cdyn_lz-adia_1d_ts_1    COMMAND ${CMAKE_MAKE_PROGRAM} cdyn_lz-adia_1d_ts_1    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)

# add expected results of random tests
set_property(TEST qdyn_1d_tully_1         PROPERTY PASS_REGULAR_EXPRESSION "ADIABATIC STATE 0 POP: 0.410505")
set_property(TEST cdyn_lz-dia_1d_tully_1  PROPERTY PASS_REGULAR_EXPRESSION "ADIABATIC STATE 0 POP: 0.509000")
set_property(TEST cdyn_lz-adia_1d_tully_1 PROPERTY PASS_REGULAR_EXPRESSION "ADIABATIC STATE 0 POP: 0.488000")
set_property(TEST qdyn_1d_tully_2         PROPERTY PASS_REGULAR_EXPRESSION "ADIABATIC STATE 0 POP: 0.024834")
set_property(TEST cdyn_lz-dia_1d_tully_2  PROPERTY PASS_REGULAR_EXPRESSION "ADIABATIC STATE 0 POP: 0.175000")
set_property(TEST cdyn_lz-adia_1d_tully_2 PROPERTY PASS_REGULAR_EXPRESSION "ADIABATIC STATE 0 POP: 0.142000")
set_property(TEST qdyn_1d_ds_1            PROPERTY PASS_REGULAR_EXPRESSION "ADIABATIC STATE 0 POP: 0.960702")
set_property(TEST cdyn_lz-dia_1d_ds_1     PROPERTY PASS_REGULAR_EXPRESSION "ADIABATIC STATE 0 POP: 0.916000")
set_property(TEST cdyn_lz-adia_1d_ds_1    PROPERTY PASS_REGULAR_EXPRESSION "ADIABATIC STATE 0 POP: 0.942000")
set_property(TEST qdyn_1d_ts_1            PROPERTY PASS_REGULAR_EXPRESSION "ADIABATIC STATE 0 POP: 0.720077")
set_property(TEST cdyn_lz-dia_1d_ts_1     PROPERTY PASS_REGULAR_EXPRESSION "ADIABATIC STATE 0 POP: 0.708000")
set_property(TEST cdyn_lz-adia_1d_ts_1    PROPERTY PASS_REGULAR_EXPRESSION "ADIABATIC STATE 0 POP: 0.720000")
