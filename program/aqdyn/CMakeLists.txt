# add the executable
add_executable(acorn_aqdyn src/main.cpp)

# include program include directories
target_include_directories(acorn_aqdyn PRIVATE include)

# link libraries
target_link_libraries(acorn_aqdyn acorn_base fftw3$<$<NOT:$<BOOL:WIN32>>:.a>)
