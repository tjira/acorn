# add the executable
add_executable(acorn_hf src/main.cpp src/hf.cpp src/system.cpp)

# include program include directories
target_include_directories(acorn_hf PRIVATE include)

# link libraries
target_link_libraries(acorn_hf ${LIBC10} ${LIBTORCH} ${LIBTORCH_CPU})

# link OpenMP
if(OpenMP_CXX_FOUND)
    target_link_libraries(acorn_hf OpenMP::OpenMP_CXX)
endif()
