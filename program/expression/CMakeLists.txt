# add the executable
add_executable(acorn_expression src/main.cpp ../../src/expression.cpp ../../src/timer.cpp)

# include program include directories
target_include_directories(acorn_expression PRIVATE include)

# link libraries
target_link_libraries(acorn_expression Eigen3::Eigen)
