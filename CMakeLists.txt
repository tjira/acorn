# minimum cmake requirement
cmake_minimum_required(VERSION 3.18)

# project name
project(Acorn)

# include plugins
include(FetchContent)

# set cmake variables
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)

# set compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg -Og -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-mavx -s -O3")

# declare libraries
FetchContent_Declare(argparse SYSTEM GIT_REPOSITORY https://github.com/p-ranav/argparse.git GIT_TAG ac4c578f6020d5164f66ccb26b4727dea657a12b)
FetchContent_Declare(exprtk GIT_REPOSITORY https://github.com/ArashPartow/exprtk.git GIT_TAG f46bffcd6966d38a09023fb37ba9335214c9b959)
FetchContent_Declare(json GIT_REPOSITORY https://github.com/nlohmann/json.git GIT_TAG a259ecc51e1951e12f757ce17db958e9881e9c6c)

# fetch the libraries
FetchContent_MakeAvailable(argparse exprtk json)

# include directories
include_directories(include ${exprtk_SOURCE_DIR} ${argparse_SOURCE_DIR}/include ${json_SOURCE_DIR}/include)

# find OpenMP
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)
find_package(FFTW3 REQUIRED)

# add acorn library
add_library(acorn
    src/bagel.cpp
    src/determinant.cpp
    src/eigen.cpp
    src/expression.cpp
    src/integral.cpp
    src/method.cpp
    src/modelsolver.cpp
    src/modelsystem.cpp
    src/numpy.cpp
    src/orca.cpp
    src/population.cpp
    src/printer.cpp
    src/restrictedconfigurationinteraction.cpp
    src/restrictedhartreefock.cpp
    src/restrictedmollerplesset.cpp
    src/system.cpp
    src/timer.cpp
    src/transform.cpp
    src/unrestrictedhartreefock.cpp
)
target_link_libraries(acorn int2.a Eigen3::Eigen)

# link OpenMP if release specified
if (CMAKE_BUILD_TYPE STREQUAL "Release")
    target_link_libraries(acorn OpenMP::OpenMP_CXX)
endif()

# add acorn executable
add_executable(acorn_exe src/main.cpp)
target_link_libraries(acorn_exe acorn fftw3)

# set the folder where the basis set files are located
target_compile_definitions(acorn_exe PRIVATE DATADIR="${PROJECT_SOURCE_DIR}")

# set the executable name
set_target_properties(acorn_exe PROPERTIES OUTPUT_NAME "acorn")

# add tools
add_executable(acorn_centerxyz tool/centerxyz.cpp)

# TESTING ==========================================================================================

# enable testing
enable_testing()

if (BUILD_TESTING)
    # add test executables
    add_executable(test_rhf_energy test/rhf_energy.cpp)
    add_executable(test_rmp2_energy test/rmp2_energy.cpp)
    add_executable(test_rfci_energy test/rfci_energy.cpp)
    add_executable(test_uhf_energy test/uhf_energy.cpp)

    # link the test executables to the library
    target_link_libraries(test_rhf_energy acorn)
    target_link_libraries(test_rmp2_energy acorn)
    target_link_libraries(test_rfci_energy acorn)
    target_link_libraries(test_uhf_energy acorn)

    # change the output folder for the test executables
    set_target_properties(test_rhf_energy PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/test)
    set_target_properties(test_rmp2_energy PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/test)
    set_target_properties(test_rfci_energy PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/test)
    set_target_properties(test_uhf_energy PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/test)
endif()

# include test file
include(CMakeTests.txt)
