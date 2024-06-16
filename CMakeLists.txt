# minimum cmake requirement
cmake_minimum_required(VERSION 3.16)

# project name
project(Acorn)

# set policies
set(CMAKE_POLICY_DEFAULT_CMP0153 OLD)
set(CMAKE_POLICY_DEFAULT_CMP0148 OLD)
set(CMAKE_POLICY_DEFAULT_CMP0146 OLD)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

# set the standard and compile commands export
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)

# disable deprecation warning caused by fftw
set(CMAKE_WARN_DEPRECATED OFF CACHE BOOL "" FORCE)

# set the output paths
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG   ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY         ${PROJECT_SOURCE_DIR}/bin)

# set the build variables
set(BUILD_SHARED_LIBS OFF)
set(BUILD_TESTS OFF)

# disable eigen qt tests
set(EIGEN_TEST_NOQT ON)

# set the base compiler flags
set(CMAKE_CXX_FLAGS_RELEASE "-march=native -s -O3 -funroll-loops")
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg -Og -Wall -Wextra --coverage")

if (STATIC)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -static")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -static")
endif()

# set platform specific flags
if (WIN32)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wa,-mbig-obj")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wa,-mbig-obj")
else()
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto=auto")
endif()

if (WIN32)
    # include the fetch plugin
    include(FetchContent)

    # declare the libint library on github releases
    FetchContent_Declare(libint SYSTEM DOWNLOAD_EXTRACT_TIMESTAMP TRUE URL https://github.com/evaleev/libint/releases/download/v2.8.2/libint-2.8.2-mpqc4.tgz)

    # declare the eigen release for linear algebra
    FetchContent_Declare(eigen SYSTEM DOWNLOAD_EXTRACT_TIMESTAMP TRUE URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz)

    # declare the fftw library on the fftw webpage
    FetchContent_Declare(fftw SYSTEM DOWNLOAD_EXTRACT_TIMESTAMP TRUE URL https://www.fftw.org/fftw-3.3.10.tar.gz)

    # download these libraries
    FetchContent_MakeAvailable(eigen fftw libint)

    # read the problematic libint file
    file(READ ${libint_SOURCE_DIR}/include/libint2/boys.h CONTENT)

    # fix the bug in the file
    string(REPLACE "#ifdef _MSC_VER" "#ifdef WIN32" CONTENT "${CONTENT}")

    # write back to the file
    file(WRITE ${libint_SOURCE_DIR}/include/libint2/boys.h "${CONTENT}")

    # include additional directories
    include_directories(${CMAKE_BINARY_DIR}/include ${fftw_SOURCE_DIR}/api)
endif()

# include directories
include_directories(include lib)

# find the necessary packages
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)

# add the base library
add_library(acorn_base STATIC src/expression.cpp src/fourier.cpp src/linalg.cpp src/system.cpp src/timer.cpp src/wavefunction.cpp)

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
set_property(TEST qdyn_1d_HO_imaginary PROPERTY PASS_REGULAR_EXPRESSION "FINAL WFN 2 ENERGY: 2.500001"         )
set_property(TEST qdyn_1d_HO_real      PROPERTY PASS_REGULAR_EXPRESSION "FINAL WFN 0 ENERGY: 1.124931"         )
set_property(TEST qdyn_2d_HO_imaginary PROPERTY PASS_REGULAR_EXPRESSION "FINAL WFN 2 ENERGY: 2.000001"         )
set_property(TEST qdyn_2d_HO_real      PROPERTY PASS_REGULAR_EXPRESSION "FINAL WFN 0 ENERGY: 2.249863"         )
set_property(TEST qdyn_3d_HO_imaginary PROPERTY PASS_REGULAR_EXPRESSION "FINAL WFN 2 ENERGY: 2.500001"         )
set_property(TEST qdyn_3d_HO_real      PROPERTY PASS_REGULAR_EXPRESSION "FINAL WFN 0 ENERGY: 3.374795"         )

# add random tests
add_test(NAME random_qdyn_1d_ds_1              COMMAND ${CMAKE_MAKE_PROGRAM} random_qdyn_1d_ds_1              WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME random_qdyn_1d_ts_1              COMMAND ${CMAKE_MAKE_PROGRAM} random_qdyn_1d_ts_1              WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME random_qdyn_1d_tully_1           COMMAND ${CMAKE_MAKE_PROGRAM} random_qdyn_1d_tully_1           WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME random_cdyn_1d_diabatic_ds_1     COMMAND ${CMAKE_MAKE_PROGRAM} random_cdyn_1d_diabatic_ds_1     WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME random_cdyn_1d_diabatic_ts_1     COMMAND ${CMAKE_MAKE_PROGRAM} random_cdyn_1d_diabatic_ts_1     WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME random_cdyn_1d_diabatic_tully_1  COMMAND ${CMAKE_MAKE_PROGRAM} random_cdyn_1d_diabatic_tully_1  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME random_cdyn_1d_adiabatic_ds_1    COMMAND ${CMAKE_MAKE_PROGRAM} random_cdyn_1d_adiabatic_ds_1    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME random_cdyn_1d_adiabatic_ts_1    COMMAND ${CMAKE_MAKE_PROGRAM} random_cdyn_1d_adiabatic_ts_1    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)
add_test(NAME random_cdyn_1d_adiabatic_tully_1 COMMAND ${CMAKE_MAKE_PROGRAM} random_cdyn_1d_adiabatic_tully_1 WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example)

# add expected results of random tests
set_property(TEST random_qdyn_1d_ds_1              PROPERTY PASS_REGULAR_EXPRESSION "ADIABATIC STATE 0 POP: 0.052824")
set_property(TEST random_qdyn_1d_ts_1              PROPERTY PASS_REGULAR_EXPRESSION "ADIABATIC STATE 0 POP: 0.720077")
set_property(TEST random_qdyn_1d_tully_1           PROPERTY PASS_REGULAR_EXPRESSION "ADIABATIC STATE 0 POP: 0.676839")
set_property(TEST random_cdyn_1d_adiabatic_ds_1    PROPERTY PASS_REGULAR_EXPRESSION "ADIABATIC STATE 0 POP: 0.072000")
set_property(TEST random_cdyn_1d_adiabatic_ts_1    PROPERTY PASS_REGULAR_EXPRESSION "ADIABATIC STATE 0 POP: 0.713000")
set_property(TEST random_cdyn_1d_adiabatic_tully_1 PROPERTY PASS_REGULAR_EXPRESSION "ADIABATIC STATE 0 POP: 0.504000")
set_property(TEST random_cdyn_1d_diabatic_ds_1     PROPERTY PASS_REGULAR_EXPRESSION "DIABATIC STATE 0 POP: 0.917000" )
set_property(TEST random_cdyn_1d_diabatic_ts_1     PROPERTY PASS_REGULAR_EXPRESSION "DIABATIC STATE 0 POP: 0.177000" )
set_property(TEST random_cdyn_1d_diabatic_tully_1  PROPERTY PASS_REGULAR_EXPRESSION "DIABATIC STATE 0 POP: 0.445000" )
