# add the executable
add_executable(acorn_cdyn src/main.cpp src/cdyn.cpp src/lz.cpp)

# include program include directories
target_include_directories(acorn_cdyn PRIVATE include)

# link libraries
target_link_libraries(acorn_cdyn ${LIBC10} ${LIBTORCH} ${LIBTORCH_CPU})

# link OpenMP
if(OpenMP_CXX_FOUND)
    target_link_libraries(acorn_cdyn OpenMP::OpenMP_CXX)
endif()
