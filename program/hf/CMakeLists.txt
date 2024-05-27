# add the executable
add_executable(acorn_hf src/main.cpp ../../src/system.cpp ../../src/timer.cpp)

# include program include directories
target_include_directories(acorn_hf PRIVATE include)

# link libraries
target_link_libraries(acorn_hf libint2.a Eigen3::Eigen)
