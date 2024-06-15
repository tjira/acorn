# add the executable
add_executable(acorn_nqdyn src/main.cpp)

# include program include directories
target_include_directories(acorn_nqdyn PRIVATE include)

# link libraries
target_link_libraries(acorn_nqdyn acorn_base fftw3$<$<NOT:$<BOOL:WIN32>>:.a>)
