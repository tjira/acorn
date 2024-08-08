# add the executable
add_executable(acorn_mp src/main.cpp src/mp.cpp)

# include program include directories
target_include_directories(acorn_mp PRIVATE include)

# link libraries
target_link_libraries(acorn_mp acorn_base ${LIBC10} ${LIBTORCH} ${LIBTORCH_CPU})

# link OpenMP
if(OpenMP_CXX_FOUND AND CMAKE_BUILD_TYPE STREQUAL "Release")
    target_link_libraries(acorn_mp OpenMP::OpenMP_CXX)
endif()
