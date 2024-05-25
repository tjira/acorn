# add the executable
add_executable(acorn_expression src/expression.cpp src/main.cpp)

# include program include directories
target_include_directories(acorn_expression PRIVATE include)

# link the libraries
target_link_libraries(acorn_expression acorn_base)
