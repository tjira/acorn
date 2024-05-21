# add the executable
add_executable(acorn_integral src/integral.cpp src/main.cpp)

# include program include directories
target_include_directories(acorn_integral PRIVATE include)

# link the libraries
target_link_libraries(acorn_integral acorn_base libint2.a)
