cmake_minimum_required(VERSION 3.22.2)

project(plots)

#find optional packages

find_package(Matplot++ CONFIG REQUIRED)


add_executable(plots main.cpp)

#link libraries
target_link_libraries(plots Matplot++::matplot)


#set to debug mode
set(CMAKE_C_FLAGS_DEBUG "-g -DDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG")

