# add the executable
add_executable(acorn_mp src/main.cpp)

# include program include directories
target_include_directories(acorn_mp PRIVATE include)

# link libraries
target_link_libraries(acorn_mp acorn_base)
