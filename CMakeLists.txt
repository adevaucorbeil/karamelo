cmake_minimum_required(VERSION 3.10)

project(karamelo)

include_directories(${PROJECT_SOURCE_DIR})


file(GLOB MyCSources *.cpp)
add_executable(karamelo ${MyCSources})

set(CMAKE_CXX_FLAGS "-DWITHOUT_NUMPY -O3 -march=native")
target_link_libraries(${PROJECT_NAME} -lpython2.7 -lgzstream -lz)
