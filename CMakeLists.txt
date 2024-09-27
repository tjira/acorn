# minimum cmake requirement
cmake_minimum_required(VERSION 3.16)

# project name
project(Acorn)

# set the c++ standard
set(CMAKE_CXX_STANDARD 20)

# set policies and compile commands
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# set output path
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

# set the base compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg -Og -Wall -Wextra -Wl,-rpath,'$ORIGIN' --coverage")
set(CMAKE_CXX_FLAGS_RELEASE "-s -O3 -flto=auto -march=native -Wl,-rpath,'$ORIGIN'")

# include the system directories
include_directories(SYSTEM external/include external/include/eigen3 external/include/torch/csrc/api/include)

# include the project directories
include_directories(include example/diagram)

# find the necessary packages
find_package(OpenMP REQUIRED)

# include the eigen directory
include_directories(${EIGEN3_INCLUDE_DIR})

# find installed libraries
find_library(LIBTORCH_CPU NAMES torch_cpu PATHS external/lib NO_DEFAULT_PATH)
find_library(LIBTORCH NAMES torch PATHS external/lib NO_DEFAULT_PATH)
find_library(LIBFFTW NAMES fftw PATHS external/lib NO_DEFAULT_PATH)
find_library(LIBINT NAMES int PATHS external/lib NO_DEFAULT_PATH)
find_library(LIBC10 NAMES c10 PATHS external/lib NO_DEFAULT_PATH)

# add the executable
add_executable(acorn
    src/classicaldynamics.cpp
    src/configurationinteraction.cpp
    src/coupledcluster.cpp
    src/export.cpp
    src/expression.cpp
    src/fouriertransform.cpp
    src/fewestswitches.cpp
    src/hartreefock.cpp
    src/input.cpp
    src/integral.cpp
    src/landauzener.cpp
    src/main.cpp
    src/mollerplesset.cpp
    src/quantumdynamics.cpp
    src/system.cpp
    src/timer.cpp
    src/transform.cpp
    src/utility.cpp
    src/wavefunction.cpp
)

# add compile definitions
target_compile_definitions(acorn PRIVATE DATADIR="${PROJECT_SOURCE_DIR}")

# link the libraries
target_link_libraries(acorn ${LIBC10} ${LIBFFTW} ${LIBINT} ${LIBTORCH_CPU} ${LIBTORCH})

# link OpenMP
if(OpenMP_CXX_FOUND)
    target_link_libraries(acorn OpenMP::OpenMP_CXX)
endif()

# enable testing
enable_testing()

# add general tests
add_test(NAME hf                          COMMAND ${PROJECT_SOURCE_DIR}/bin/acorn -i ${PROJECT_SOURCE_DIR}/example/input/hf.json                          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
add_test(NAME mp2                         COMMAND ${PROJECT_SOURCE_DIR}/bin/acorn -i ${PROJECT_SOURCE_DIR}/example/input/mp2.json                         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
add_test(NAME mp3                         COMMAND ${PROJECT_SOURCE_DIR}/bin/acorn -i ${PROJECT_SOURCE_DIR}/example/input/mp3.json                         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
add_test(NAME mp4                         COMMAND ${PROJECT_SOURCE_DIR}/bin/acorn -i ${PROJECT_SOURCE_DIR}/example/input/mp4.json                         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
add_test(NAME fci                         COMMAND ${PROJECT_SOURCE_DIR}/bin/acorn -i ${PROJECT_SOURCE_DIR}/example/input/fci.json                         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
add_test(NAME ccd                         COMMAND ${PROJECT_SOURCE_DIR}/bin/acorn -i ${PROJECT_SOURCE_DIR}/example/input/ccd.json                         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
add_test(NAME ccd-t                       COMMAND ${PROJECT_SOURCE_DIR}/bin/acorn -i ${PROJECT_SOURCE_DIR}/example/input/ccd-t.json                       WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
add_test(NAME ccsd                        COMMAND ${PROJECT_SOURCE_DIR}/bin/acorn -i ${PROJECT_SOURCE_DIR}/example/input/ccsd.json                        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
add_test(NAME ccsd-t                      COMMAND ${PROJECT_SOURCE_DIR}/bin/acorn -i ${PROJECT_SOURCE_DIR}/example/input/ccsd-t.json                      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
add_test(NAME qdyn_1d_imaginary_one-state COMMAND ${PROJECT_SOURCE_DIR}/bin/acorn -i ${PROJECT_SOURCE_DIR}/example/input/qdyn_1d_imaginary_one-state.json WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
add_test(NAME qdyn_2d_imaginary_one-state COMMAND ${PROJECT_SOURCE_DIR}/bin/acorn -i ${PROJECT_SOURCE_DIR}/example/input/qdyn_2d_imaginary_one-state.json WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
add_test(NAME qdyn_1d_real_one-state      COMMAND ${PROJECT_SOURCE_DIR}/bin/acorn -i ${PROJECT_SOURCE_DIR}/example/input/qdyn_1d_real_one-state.json      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
add_test(NAME qdyn_2d_real_one-state      COMMAND ${PROJECT_SOURCE_DIR}/bin/acorn -i ${PROJECT_SOURCE_DIR}/example/input/qdyn_2d_real_one-state.json      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
add_test(NAME qdyn_1d_real_two-state      COMMAND ${PROJECT_SOURCE_DIR}/bin/acorn -i ${PROJECT_SOURCE_DIR}/example/input/qdyn_1d_real_two-state.json      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
add_test(NAME cdyn_1d_two-state_fs-adia   COMMAND ${PROJECT_SOURCE_DIR}/bin/acorn -i ${PROJECT_SOURCE_DIR}/example/input/cdyn_1d_two-state_fs-adia.json   WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
add_test(NAME cdyn_1d_two-state_lz-adia   COMMAND ${PROJECT_SOURCE_DIR}/bin/acorn -i ${PROJECT_SOURCE_DIR}/example/input/cdyn_1d_two-state_lz-adia.json   WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
add_test(NAME cdyn_1d_two-state_lz-dia    COMMAND ${PROJECT_SOURCE_DIR}/bin/acorn -i ${PROJECT_SOURCE_DIR}/example/input/cdyn_1d_two-state_lz-dia.json    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

# add expected results of general tests
set_property(TEST hf                          PROPERTY PASS_REGULAR_EXPRESSION "FINAL HF ENERGY: -74.965901"                            )
set_property(TEST mp2                         PROPERTY PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -75.004854"                           )
set_property(TEST mp3                         PROPERTY PASS_REGULAR_EXPRESSION "FINAL MP3 ENERGY: -75.015455"                           )
set_property(TEST mp4                         PROPERTY PASS_REGULAR_EXPRESSION "FINAL MP4 ENERGY: -75.018703"                           )
set_property(TEST fci                         PROPERTY PASS_REGULAR_EXPRESSION "FINAL FCI ENERGY: -75.020410"                           )
set_property(TEST ccd                         PROPERTY PASS_REGULAR_EXPRESSION "FINAL CCD ENERGY: -75.020006"                           )
set_property(TEST ccd-t                       PROPERTY PASS_REGULAR_EXPRESSION "FINAL CCD\\(T\\) ENERGY: -75.020084"                    )
set_property(TEST ccsd                        PROPERTY PASS_REGULAR_EXPRESSION "FINAL CCSD ENERGY: -75.020284"                          )
set_property(TEST ccsd-t                      PROPERTY PASS_REGULAR_EXPRESSION "FINAL CCSD\\(T\\) ENERGY: -75.020351"                   )
set_property(TEST qdyn_1d_imaginary_one-state PROPERTY PASS_REGULAR_EXPRESSION "FINAL ENERGY OF WAVEFUNCTION 03 \\(ITP\\):     2.500001")
set_property(TEST qdyn_2d_imaginary_one-state PROPERTY PASS_REGULAR_EXPRESSION "FINAL ENERGY OF WAVEFUNCTION 03 \\(ITP\\):     2.000001")
set_property(TEST qdyn_1d_real_one-state      PROPERTY PASS_REGULAR_EXPRESSION "FINAL ENERGY OF WAVEFUNCTION 01 \\(RTP\\):     1.124931")
set_property(TEST qdyn_2d_real_one-state      PROPERTY PASS_REGULAR_EXPRESSION "FINAL ENERGY OF WAVEFUNCTION 01 \\(RTP\\):     2.249863")
set_property(TEST qdyn_1d_real_two-state      PROPERTY PASS_REGULAR_EXPRESSION "FINAL DIABATIC  STATE 01 POPULATION \\(RTP\\): 0.039297")
set_property(TEST cdyn_1d_two-state_fs-adia   PROPERTY PASS_REGULAR_EXPRESSION "FINAL ADIABATIC POPULATION OF STATE 01: 0.965000"       )
set_property(TEST cdyn_1d_two-state_lz-adia   PROPERTY PASS_REGULAR_EXPRESSION "FINAL ADIABATIC POPULATION OF STATE 01: 0.934000"       )
set_property(TEST cdyn_1d_two-state_lz-dia    PROPERTY PASS_REGULAR_EXPRESSION "FINAL DIABATIC POPULATION OF STATE 01: 0.069000"        )
