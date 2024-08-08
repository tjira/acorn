# add the executable
add_executable(acorn_cc src/main.cpp src/cc.cpp)

# include program include directories
target_include_directories(acorn_cc PRIVATE include)

# link libraries
target_link_libraries(acorn_cc acorn_base ${LIBC10} ${LIBTORCH} ${LIBTORCH_CPU})
