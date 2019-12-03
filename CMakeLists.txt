cmake_minimum_required(VERSION 3.12)

set (CMAKE_CXX_STANDARD 11)



project(karamelo)

include_directories(${PROJECT_SOURCE_DIR})


file(GLOB MyCSources *.cpp)
add_executable(karamelo ${MyCSources})

find_package(Boost 1.67.0 COMPONENTS system filesystem)
find_package(OpenMP)
include_directories(${Boost_INCLUDE_DIR})

set(CMAKE_CXX_FLAGS "-DWITHOUT_NUMPY -O3 -march=native")
if(OpenMP_CXX_FOUND)
  target_link_libraries(${PROJECT_NAME} PUBLIC -lpython2.7 -lgzstream -lz PUBLIC OpenMP::OpenMP_CXX ${Boost_LIBRARIES})
endif()


