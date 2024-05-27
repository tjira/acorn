# add the executable
add_executable(acorn_aqdyn src/main.cpp ../../src/fourier.cpp ../../src/timer.cpp ../../src/wavefunction.cpp)

# include program include directories
target_include_directories(acorn_aqdyn PRIVATE include)

# link libraries
target_link_libraries(acorn_aqdyn fftw3.a Eigen3::Eigen)
