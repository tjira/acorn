# add the executable
add_executable(acorn_ci src/main.cpp)

# include program include directories
target_include_directories(acorn_ci PRIVATE include)

# link libraries
target_link_libraries(acorn_ci acorn_base)
