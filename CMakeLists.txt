cmake_minimum_required(VERSION 3.5)

set (CMAKE_CXX_STANDARD 17)

find_package(Git)
if(Git_FOUND)
  message("Git found: ${GIT_EXECUTABLE}")
endif()

project(karamelo)

find_library(LIBNVTOOLSEXT nvToolsExt PATHS ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
include_directories(${PROJECT_SOURCE_DIR})

set(CMAKE_CUDA_COMPILER "/apps/nvidia-hpc-sdk/21.3/Linux_x86_64/21.3/compilers/bin/nvcc")
set(CMAKE_CUDA_COMPILER_VERSION "10.2.89")
set(CMAKE_CUDA_COMPILER_TOOLKIT_ROOT "/apps/nvidia-hpc-sdk/21.3/Linux_x86_64/21.3/cuda/")
set(CMAKE_CUDA_COMPILER_LIBRARY_ROOT "/apps/nvidia-hpc-sdk/21.3/Linux_x86_64/21.3/cuda/")

set(CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES "/apps/nvidia-hpc-sdk/21.3/Linux_x86_64/21.3/cuda/10.2/include")

# the commit's SHA1, and whether the building workspace was dirty or not
execute_process(COMMAND
  "${GIT_EXECUTABLE}" describe --match=NeVeRmAtCh --always --abbrev=40 --dirty
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_SHA1
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

# the date of the commit
execute_process(COMMAND
  "${GIT_EXECUTABLE}" log -1 --format=%ad --date=local
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_DATE
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

# the subject of the commit
execute_process(COMMAND
  "${GIT_EXECUTABLE}" log -1 --format=%s
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_COMMIT_SUBJECT
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

# generate version.cpp
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/version.cpp.in" "${CMAKE_CURRENT_SOURCE_DIR}/src/version.cpp" @ONLY)

add_definitions(-D_POSIX_C_SOURCE=200809L)

file(GLOB MyCSources src/*.cpp)
add_executable(karamelo ${MyCSources})

set(EIGEN_BUILD_DOC OFF CACHE INTERNAL "")

add_subdirectory(doc EXCLUDE_FROM_ALL)
add_subdirectory(graphics)
add_subdirectory(kokkos_demo)
add_subdirectory(open_acc_demo)
add_subdirectory(third-party/eigen)
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

target_link_libraries(${PROJECT_NAME} PUBLIC ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS} gzstream eigen)
