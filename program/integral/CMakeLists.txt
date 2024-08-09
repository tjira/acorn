# add the executable
add_executable(acorn_integral src/integral.cpp src/main.cpp)

# include program include directories
target_include_directories(acorn_integral PRIVATE include)

# add compile definitions
target_compile_definitions(acorn_integral PRIVATE DATADIR="${PROJECT_SOURCE_DIR}")

# link the libraries
target_link_libraries(acorn_integral ${LIBC10} ${LIBINT} ${LIBTORCH} ${LIBTORCH_CPU})

# link OpenMP
if(OpenMP_CXX_FOUND)
    target_link_libraries(acorn_integral OpenMP::OpenMP_CXX)
endif()
