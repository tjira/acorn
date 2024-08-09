# add the executable
add_executable(acorn_hf src/main.cpp src/system.cpp)

# include program include directories
target_include_directories(acorn_hf PRIVATE include)

# link OpenMP
if(OpenMP_CXX_FOUND)
    target_link_libraries(acorn_hf OpenMP::OpenMP_CXX)
endif()
