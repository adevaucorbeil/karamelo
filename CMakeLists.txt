cmake_minimum_required(VERSION 3.10)


project(karamelo)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
include_directories(${PROJECT_SOURCE_DIR})


file(GLOB MyCSources *.cpp)
add_executable(karamelo ${MyCSources})

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

target_compile_options(${PROJECT_NAME} PRIVATE ${MPI_CXX_COMPILE_FLAGS} "-DWITHOUT_NUMPY -march=native")
target_link_libraries(${PROJECT_NAME} ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS} -lpython2.7 -lgzstream -lz)
