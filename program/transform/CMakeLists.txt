# add the executable
add_executable(acorn_transform src/main.cpp src/transform.cpp)

# include program include directories
target_include_directories(acorn_transform PRIVATE include)
