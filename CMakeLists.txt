cmake_minimum_required(VERSION 3.5)

set (CMAKE_CXX_STANDARD 17)

project(karamelo)

#find_library(LIBNVTOOLSEXT nvToolsExt PATHS ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

add_definitions(-D_POSIX_C_SOURCE=200809L)

add_subdirectory(doc EXCLUDE_FROM_ALL)

add_subdirectory(src)
add_subdirectory(graphics)
add_subdirectory(kokkos_demo)
add_subdirectory(third-party)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

target_compile_options(${PROJECT_NAME} PRIVATE ${MPI_CXX_COMPILE_FLAGS})
if(WIN32)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")
  set(CMAKE_CXX_FLAGS_PROFILING "${CMAKE_CXX_FLAGS_PROFILING} -O2 -g")
else()
  target_compile_options(${PROJECT_NAME} PRIVATE "-march=native")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
  set(CMAKE_CXX_FLAGS_PROFILING "${CMAKE_CXX_FLAGS_PROFILING} -O3 -g")
endif()
