# add the executable
add_executable(acorn_ci src/main.cpp src/ci.cpp)

# include program include directories
target_include_directories(acorn_ci PRIVATE include)

# link libraries
target_link_libraries(acorn_ci ${LIBC10} ${LIBTORCH} ${LIBTORCH_CPU})

# link OpenMP
if(OpenMP_CXX_FOUND)
    target_link_libraries(acorn_ci OpenMP::OpenMP_CXX)
endif()
