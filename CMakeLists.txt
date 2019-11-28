cmake_minimum_required(VERSION 3.10)


project(karamelo)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
include_directories(${PROJECT_SOURCE_DIR})


file(GLOB MyCSources *.cpp)
add_executable(karamelo ${MyCSources})

target_compile_options(${PROJECT_NAME} PRIVATE ${MPI_CXX_COMPILE_FLAGS} "-DWITHOUT_NUMPY -O3 -march=native")
target_link_libraries(${PROJECT_NAME} ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS} -lpython2.7 -lgzstream -lz)
