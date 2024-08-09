# add the executable
add_executable(acorn_qdyn src/main.cpp src/fourier.cpp src/wavefunction.cpp)

# include program include directories
target_include_directories(acorn_qdyn PRIVATE include)

# link libraries
target_link_libraries(acorn_qdyn ${LIBFFTW})
