# add the executable
add_executable(acorn_integral src/integral.cpp src/main.cpp)

# include program include directories
target_include_directories(acorn_integral PRIVATE include)

# add compile definitions
target_compile_definitions(acorn_integral PRIVATE DATADIR="${PROJECT_SOURCE_DIR}")

# link the libraries
target_link_libraries(acorn_integral acorn_base libint2$<$<NOT:$<BOOL:WIN32>>:.a>)
