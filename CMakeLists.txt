cmake_minimum_required(VERSION 3.5)

set (CMAKE_CXX_STANDARD 14)

project(karamelo)

find_library(LIBNVTOOLSEXT nvToolsExt PATHS ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set(CMAKE_CUDA_COMPILER "/apps/nvidia-hpc-sdk/21.3/Linux_x86_64/21.3/compilers/bin/nvcc")
set(CMAKE_CUDA_COMPILER_VERSION "10.2.89")
set(CMAKE_CUDA_COMPILER_TOOLKIT_ROOT "/apps/nvidia-hpc-sdk/21.3/Linux_x86_64/21.3/cuda/")
set(CMAKE_CUDA_COMPILER_LIBRARY_ROOT "/apps/nvidia-hpc-sdk/21.3/Linux_x86_64/21.3/cuda/")

set(CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES "/apps/nvidia-hpc-sdk/21.3/Linux_x86_64/21.3/cuda/10.2/include")

add_definitions(-D_POSIX_C_SOURCE=200809L)

set(EIGEN_BUILD_DOC OFF CACHE INTERNAL "")
set(GLFW_BUILD_DOCS OFF CACHE INTERNAL "")

add_subdirectory(src)
add_subdirectory(doc EXCLUDE_FROM_ALL)
add_subdirectory(graphics)
add_subdirectory(kokkos_demo)
add_subdirectory(third-party/glad)
add_subdirectory(third-party/glfw)
add_subdirectory(third-party/gzstream)
add_subdirectory(third-party/kokkos)

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
