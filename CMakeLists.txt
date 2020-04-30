cmake_minimum_required(VERSION 3.10)

set (CMAKE_CXX_STANDARD 14)
find_package(Git)
if(Git_FOUND)
  message("Git found: ${GIT_EXECUTABLE}")
endif()

project(karamelo)

find_package(MPI REQUIRED)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${MPI_INCLUDE_PATH})
include_directories(${PROJECT_SOURCE_DIR})


file(GLOB MyCSources src/*.cpp)
add_executable(karamelo ${MyCSources})

find_package(Boost 1.67.0 COMPONENTS system filesystem)
find_package(OpenMP)
include_directories(${Boost_INCLUDE_DIR})

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")



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
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/version.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/src/version.cpp" @ONLY)

list(APPEND SOURCES "${CMAKE_CURRENT_BINARY_DIR}/src/version.cpp" version.h)


target_compile_options(${PROJECT_NAME} PRIVATE ${MPI_CXX_COMPILE_FLAGS} "-DWITHOUT_NUMPY -march=native")

if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS} -lpython2.7 -lgzstream -lz PUBLIC OpenMP::OpenMP_CXX ${Boost_LIBRARIES} Eigen3::Eigen)
endif()

