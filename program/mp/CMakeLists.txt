# add the executable
add_executable(acorn_mp src/main.cpp ../../src/system.cpp ../../src/timer.cpp)

# include program include directories
target_include_directories(acorn_mp PRIVATE include)

# set the compile definitions
target_compile_definitions(acorn_mp PUBLIC PTABLE="${PTABLE}")

# link libraries
target_link_libraries(acorn_mp Eigen3::Eigen)
