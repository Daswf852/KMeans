cmake_minimum_required(VERSION 3.16)

include_directories(
    /usr/include
    /usr/include/opencv4
    lib/
)

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -ggdb")

project(palette)
find_package(OpenCV REQUIRED)
project(bench)
find_package(OpenCV REQUIRED)

add_executable(palette main.cpp)
target_link_libraries(palette ${OpenCV_LIBS} pthread)

add_executable(bench benchmark.cpp)
target_link_libraries(bench pthread)
