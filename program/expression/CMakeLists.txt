# add the executable
add_executable(acorn_expression src/main.cpp)

# include program include directories
target_include_directories(acorn_expression PRIVATE include)
