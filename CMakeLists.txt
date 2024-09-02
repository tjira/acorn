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
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg -Og -Wall -Wextra -Wl,-rpath,'$ORIGIN' --coverage -flarge-source-files")
set(CMAKE_CXX_FLAGS_RELEASE "-s -O3 -flto=auto -march=native -Wl,-rpath,'$ORIGIN' -flarge-source-files")

# include the system directories
include_directories(SYSTEM external/include external/include/torch/csrc/api/include)

# include the project directories
include_directories(include example/diagram)

# find the necessary packages
find_package(Eigen3 REQUIRED)
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
    src/configurationinteraction.cpp
    src/coupledcluster.cpp
    src/export.cpp
    src/expression.cpp
    src/fouriertransform.cpp
    src/hartreefock.cpp
    src/input.cpp
    src/integral.cpp
    src/main.cpp
    src/mollerplesset.cpp
    src/quantumdynamics.cpp
    src/system.cpp
    src/timer.cpp
    src/transform.cpp
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
