cmake_minimum_required(VERSION 3.13)
project(EhrenfestModel)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

add_subdirectory(extern/computational_physics_lib)
add_subdirectory(src)
