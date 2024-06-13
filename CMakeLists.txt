# minimum cmake requirement
cmake_minimum_required(VERSION 3.16)

# project name
project(Acorn)

# set output path, C++ standard, and compile commands export
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)

# set the compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg -Og -Wall -Wextra -fprofile-arcs -ftest-coverage")
set(CMAKE_CXX_FLAGS_RELEASE "-march=native -s -O3 -funroll-loops -flto=auto")

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
add_test(NAME hf                   COMMAND make -C ${PROJECT_SOURCE_DIR}/example hf                  )
add_test(NAME mp2                  COMMAND make -C ${PROJECT_SOURCE_DIR}/example mp2                 )
add_test(NAME mp3                  COMMAND make -C ${PROJECT_SOURCE_DIR}/example mp3                 )
add_test(NAME fci                  COMMAND make -C ${PROJECT_SOURCE_DIR}/example fci                 )
add_test(NAME aqdyn_1d_imaginary   COMMAND make -C ${PROJECT_SOURCE_DIR}/example aqdyn_1d_imaginary  )
add_test(NAME aqdyn_1d_real        COMMAND make -C ${PROJECT_SOURCE_DIR}/example aqdyn_1d_real       )
add_test(NAME aqdyn_2d_imaginary   COMMAND make -C ${PROJECT_SOURCE_DIR}/example aqdyn_2d_imaginary  )
add_test(NAME aqdyn_2d_real        COMMAND make -C ${PROJECT_SOURCE_DIR}/example aqdyn_2d_real       )
add_test(NAME aqdyn_3d_imaginary   COMMAND make -C ${PROJECT_SOURCE_DIR}/example aqdyn_3d_imaginary  )
add_test(NAME aqdyn_3d_real        COMMAND make -C ${PROJECT_SOURCE_DIR}/example aqdyn_3d_real       )

# add expected results of general tests
set_property(TEST hf                   PROPERTY PASS_REGULAR_EXPRESSION "FINAL SINGLE POINT ENERGY: -74.965901")
set_property(TEST mp2                  PROPERTY PASS_REGULAR_EXPRESSION "FINAL SINGLE POINT ENERGY: -75.004854")
set_property(TEST mp3                  PROPERTY PASS_REGULAR_EXPRESSION "FINAL SINGLE POINT ENERGY: -75.015455")
set_property(TEST fci                  PROPERTY PASS_REGULAR_EXPRESSION "FINAL SINGLE POINT ENERGY: -75.020410")
set_property(TEST aqdyn_1d_imaginary   PROPERTY PASS_REGULAR_EXPRESSION "FINAL WFN 2 ENERGY: 2.500001"         )
set_property(TEST aqdyn_1d_real        PROPERTY PASS_REGULAR_EXPRESSION "FINAL WFN 0 ENERGY: 0.625207"         )
set_property(TEST aqdyn_2d_imaginary   PROPERTY PASS_REGULAR_EXPRESSION "FINAL WFN 2 ENERGY: 2.000001"         )
set_property(TEST aqdyn_2d_real        PROPERTY PASS_REGULAR_EXPRESSION "FINAL WFN 0 ENERGY: 1.250415"         )
set_property(TEST aqdyn_3d_imaginary   PROPERTY PASS_REGULAR_EXPRESSION "FINAL WFN 2 ENERGY: 2.500001"         )
set_property(TEST aqdyn_3d_real        PROPERTY PASS_REGULAR_EXPRESSION "FINAL WFN 0 ENERGY: 1.875623"         )

# add random tests
add_test(NAME random_qdyn_1d_ds_1              COMMAND make -C ${PROJECT_SOURCE_DIR}/example random_qdyn_1d_ds_1             )
add_test(NAME random_qdyn_1d_ts_1              COMMAND make -C ${PROJECT_SOURCE_DIR}/example random_qdyn_1d_ts_1             )
add_test(NAME random_qdyn_1d_tully_1           COMMAND make -C ${PROJECT_SOURCE_DIR}/example random_qdyn_1d_tully_1          )
add_test(NAME random_cdyn_1d_diabatic_ds_1     COMMAND make -C ${PROJECT_SOURCE_DIR}/example random_cdyn_1d_diabatic_ds_1    )
add_test(NAME random_cdyn_1d_diabatic_ts_1     COMMAND make -C ${PROJECT_SOURCE_DIR}/example random_cdyn_1d_diabatic_ts_1    )
add_test(NAME random_cdyn_1d_diabatic_tully_1  COMMAND make -C ${PROJECT_SOURCE_DIR}/example random_cdyn_1d_diabatic_tully_1 )
add_test(NAME random_cdyn_1d_adiabatic_ds_1    COMMAND make -C ${PROJECT_SOURCE_DIR}/example random_cdyn_1d_adiabatic_ds_1   )
add_test(NAME random_cdyn_1d_adiabatic_ts_1    COMMAND make -C ${PROJECT_SOURCE_DIR}/example random_cdyn_1d_adiabatic_ts_1   )
add_test(NAME random_cdyn_1d_adiabatic_tully_1 COMMAND make -C ${PROJECT_SOURCE_DIR}/example random_cdyn_1d_adiabatic_tully_1)

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
