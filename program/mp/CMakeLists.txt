# add the executable
add_executable(acorn_mp src/main.cpp ../../src/linalg.cpp ../../src/system.cpp ../../src/timer.cpp)

# include program include directories
target_include_directories(acorn_mp PRIVATE include)

# link libraries
target_link_libraries(acorn_mp Eigen3::Eigen)
