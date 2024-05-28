# add the executable
add_executable(acorn_hf src/main.cpp ../../src/system.cpp ../../src/timer.cpp)

# include program include directories
target_include_directories(acorn_hf PRIVATE include)

# set the compile definitions
target_compile_definitions(acorn_hf PUBLIC PTABLE="${PTABLE}")

# link libraries
target_link_libraries(acorn_hf Eigen3::Eigen)
