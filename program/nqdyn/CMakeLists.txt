# add the executable
add_executable(acorn_nqdyn src/main.cpp)

# include program include directories
target_include_directories(acorn_nqdyn PRIVATE include)

# link the libraries
target_link_libraries(acorn_nqdyn acorn_base)
