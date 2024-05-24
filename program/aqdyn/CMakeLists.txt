# add the executable
add_executable(acorn_aqdyn src/main.cpp)

# include program include directories
target_include_directories(acorn_aqdyn PRIVATE include)

# link the libraries
target_link_libraries(acorn_aqdyn acorn_base)
