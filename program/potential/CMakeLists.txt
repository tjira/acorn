# add the executable
add_executable(acorn_potential src/main.cpp)

# include program include directories
target_include_directories(acorn_potential PRIVATE include)

# link the libraries
target_link_libraries(acorn_potential acorn_base)
