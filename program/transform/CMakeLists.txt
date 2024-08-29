# add the executable
add_executable(acorn_transform src/main.cpp src/transform.cpp)

# include program include directories
target_include_directories(acorn_transform PRIVATE include)

# link the libraries
target_link_libraries(acorn_transform ${LIBC10} ${LIBTORCH} ${LIBTORCH_CPU})

# link OpenMP
if(OpenMP_CXX_FOUND)
    target_link_libraries(acorn_transform OpenMP::OpenMP_CXX)
endif()
