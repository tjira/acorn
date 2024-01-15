# minimum cmake requirement
cmake_minimum_required(VERSION 3.18)

# project name
project(Acorn)

# include plugins
include(FetchContent)

# set cmake variables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)

# set compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg -Og -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-mavx -s -O3")

# declare libraries
FetchContent_Declare(argparse SYSTEM GIT_REPOSITORY https://github.com/p-ranav/argparse.git GIT_TAG ac4c578f6020d5164f66ccb26b4727dea657a12b)
FetchContent_Declare(json GIT_REPOSITORY https://github.com/nlohmann/json.git GIT_TAG a259ecc51e1951e12f757ce17db958e9881e9c6c)

# fetch the libraries
FetchContent_MakeAvailable(argparse json)

# include directories
include_directories(include ${argparse_SOURCE_DIR}/include ${json_SOURCE_DIR}/include)

# find OpenMP
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)

# add acorn executable
add_executable(acorn
    src/integral.cpp
    src/main.cpp
    src/method.cpp
    src/restrictedhartreefock.cpp
    src/restrictedmollerplesset.cpp
    src/system.cpp
    src/timer.cpp
    src/transform.cpp
)
target_compile_definitions(acorn PRIVATE DATADIR="${PROJECT_SOURCE_DIR}")
target_link_libraries(acorn int2.a Eigen3::Eigen OpenMP::OpenMP_CXX)
