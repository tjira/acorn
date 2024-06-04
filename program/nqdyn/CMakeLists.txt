# add the executable
add_executable(acorn_nqdyn src/main.cpp ../../src/fourier.cpp ../../src/linalg.cpp ../../src/timer.cpp ../../src/wavefunction.cpp)

# include program include directories
target_include_directories(acorn_nqdyn PRIVATE include)

# link libraries
target_link_libraries(acorn_nqdyn fftw3.a Eigen3::Eigen)
