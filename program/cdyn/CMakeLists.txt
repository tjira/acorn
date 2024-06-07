# add the executable
add_executable(acorn_cdyn src/main.cpp src/landauzener.cpp)

# include program include directories
target_include_directories(acorn_cdyn PRIVATE include)

# link libraries
target_link_libraries(acorn_cdyn acorn_base)
