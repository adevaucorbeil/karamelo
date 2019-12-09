cmake_minimum_required(VERSION 3.10)

set (CMAKE_CXX_STANDARD 14)


project(karamelo)

find_package(MPI REQUIRED)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${MPI_INCLUDE_PATH})
include_directories(${PROJECT_SOURCE_DIR})


file(GLOB MyCSources *.cpp)
add_executable(karamelo ${MyCSources})

find_package(Boost 1.67.0 COMPONENTS system filesystem)
find_package(OpenMP)
include_directories(${Boost_INCLUDE_DIR})

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

target_compile_options(${PROJECT_NAME} PRIVATE ${MPI_CXX_COMPILE_FLAGS} "-DWITHOUT_NUMPY -march=native")

if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS} -lpython2.7 -lgzstream -lz PUBLIC OpenMP::OpenMP_CXX ${Boost_LIBRARIES} Eigen3::Eigen)
endif()
