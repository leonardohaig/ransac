cmake_minimum_required(VERSION 3.10)
project(ransac)

set(CMAKE_CXX_STANDARD 11)

find_package(OpenCV REQUIRED)

aux_source_directory(./ SOURCE_FILE)

#add_executable(ransac main.cpp)
add_executable(ransac ${SOURCE_FILE})

target_link_libraries(ransac ${OpenCV_LIBS})


