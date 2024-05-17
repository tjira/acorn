# minimum cmake requirement
cmake_minimum_required(VERSION 3.18)

# project name
project(Acorn)

# set policies
set(CMAKE_POLICY_DEFAULT_CMP0148 OLD)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

# set the cmake variables
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_WARN_DEPRECATED OFF CACHE BOOL "" FORCE)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)

# set the build variables
set(BUILD_SHARED_LIBS OFF)
set(BUILD_TESTS OFF)

# disable eigen qt tests
set(EIGEN_TEST_NOQT ON)

# set glfw cmake flags
set(GLFW_BUILD_EXAMPLES OFF)
set(GLFW_BUILD_TESTS OFF)
set(GLFW_BUILD_DOCS OFF)

# set compile flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg -Og -Wall -Wextra -flarge-source-files -fprofile-arcs -ftest-coverage")
set(CMAKE_CXX_FLAGS_RELEASE "-mavx -s -O3 -Wall -Wextra -flarge-source-files")
if (WIN32)
    string(APPEND CMAKE_CXX_FLAGS_RELEASE " -mwindows -static -Wa,-mbig-obj")
    string(APPEND CMAKE_CXX_FLAGS_DEBUG " -mwindows -static -Wa,-mbig-obj")
endif()

# includes and testing
include(FetchContent)
enable_testing()

# declare libraries
FetchContent_Declare(libint SYSTEM DOWNLOAD_EXTRACT_TIMESTAMP TRUE URL https://github.com/evaleev/libint/releases/download/v2.8.2/libint-2.8.2-mpqc4.tgz)
FetchContent_Declare(argparse SYSTEM GIT_REPOSITORY https://github.com/p-ranav/argparse.git GIT_TAG ac4c578f6020d5164f66ccb26b4727dea657a12b)
FetchContent_Declare(exprtk SYSTEM GIT_REPOSITORY https://github.com/ArashPartow/exprtk.git GIT_TAG f46bffcd6966d38a09023fb37ba9335214c9b959)
FetchContent_Declare(eigen SYSTEM GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git GIT_TAG 6f9ad7da6122fdb4197c0b43dfec09ec3525305e)
FetchContent_Declare(json SYSTEM GIT_REPOSITORY https://github.com/nlohmann/json.git GIT_TAG a259ecc51e1951e12f757ce17db958e9881e9c6c)
FetchContent_Declare(glad SYSTEM GIT_REPOSITORY https://github.com/Dav1dde/glad.git GIT_TAG 2348b07c1ab4504d60398713781d8a57880234fa)
FetchContent_Declare(glfw SYSTEM GIT_REPOSITORY https://github.com/glfw/glfw.git GIT_TAG 3eaf1255b29fdf5c2895856c7be7d7185ef2b241)
FetchContent_Declare(glm SYSTEM GIT_REPOSITORY https://github.com/g-truc/glm.git GIT_TAG 47585fde0c49fa77a2bf2fb1d2ead06999fd4b6e)
FetchContent_Declare(fftw SYSTEM DOWNLOAD_EXTRACT_TIMESTAMP TRUE URL https://www.fftw.org/fftw-3.3.10.tar.gz)

# fetch the libraries
FetchContent_MakeAvailable(argparse exprtk fftw json)
if (STANDALONE)
    FetchContent_MakeAvailable(eigen libint)
endif()

# include directories
include_directories(include ${CMAKE_BINARY_DIR}/include ${exprtk_SOURCE_DIR} ${fftw_SOURCE_DIR}/api ${argparse_SOURCE_DIR}/include ${json_SOURCE_DIR}/include)

# find the necessary packages
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)

# add acorn library
add_library(acorn
    src/bagel.cpp
    src/determinant.cpp
    src/eigen.cpp
    src/expression.cpp
    src/integral.cpp
    src/matrixofmatrices.cpp
    src/method.cpp
    src/modelsolver.cpp
    src/modelsystem.cpp
    src/numpy.cpp
    src/orca.cpp
    src/population.cpp
    src/printer.cpp
    src/propagator.cpp
    src/restrictedconfigurationinteraction.cpp
    src/restrictedhartreefock.cpp
    src/restrictedmollerplesset.cpp
    src/system.cpp
    src/timer.cpp
    src/transform.cpp
    src/unrestrictedhartreefock.cpp
)

# link libraries to the acorn library
target_link_libraries(acorn fftw3 libint2$<$<NOT:$<BOOL:${STANDALONE}>>:.a> Eigen3::Eigen $<$<CONFIG:Release>:OpenMP::OpenMP_CXX>)

# add acorn executable
add_executable(acorn_exe src/main.cpp)
target_link_libraries(acorn_exe acorn)

# set the folder where the basis set files are located
target_compile_definitions(acorn_exe PRIVATE DATADIR="${PROJECT_SOURCE_DIR}")
target_compile_definitions(acorn PRIVATE DATADIR="${PROJECT_SOURCE_DIR}")

# set the executable name
set_target_properties(acorn_exe PROPERTIES OUTPUT_NAME "acorn")

if (BUILD_TOOLS)
    # fetch the libraries for tools
    FetchContent_MakeAvailable(glad glfw glm)

    # generate glad library
    add_subdirectory(${glad_SOURCE_DIR}/cmake)
    glad_add_library(glad API gl:core=4.2)

    # add tools
    add_executable(acorn_centerxyz tool/centerxyz.cpp)
    add_executable(acorn_makemov tool/makemov.cpp)
    add_executable(acorn_view tool/view.cpp)

    # link the tools
    target_link_libraries(acorn_view glad glfw glm)
endif()

# list of test names
set(TESTS
    test/rhf.cpp
    test/rhf_mulliken.cpp
    test/rhf_grad.cpp
    test/rmp2.cpp
    test/rmp3.cpp
    test/rfci.cpp
    test/rcisd.cpp
    test/rcid.cpp
    test/rcis.cpp
    test/uhf.cpp
    test/qad1d_imag.cpp
    test/qad1d_real.cpp
    test/qad2d_imag.cpp
    test/qad2d_real.cpp
    test/qnd1d_real.cpp
    test/nad1d.cpp
)

# create the test sourcelist
create_test_sourcelist(TESTLIST acorn_test.cpp ${TESTS})

if (BUILD_TESTING)
    # add the test executable and link the libraries
    add_executable(acorn_test acorn_test.cpp ${TESTS})
    target_link_libraries(acorn_test acorn)

    # generate the tests
    foreach (TEST ${TESTS})
        get_filename_component(TESTNAME ${TEST} NAME_WE)
        add_test(NAME test_${TESTNAME} COMMAND acorn_test "test/${TESTNAME}")
    endforeach ()
endif()
