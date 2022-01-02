file(GLOB_RECURSE SOURCE_FILES *.cpp)
add_executable(karamelo ${SOURCE_FILES})

file(GLOB_RECURSE HEADER_FILES *.h)
foreach(HEADER_FILE ${HEADER_FILES})
  get_filename_component(HEADER_DIRECTORY ${HEADER_FILE} PATH)
  list(APPEND HEADER_DIRECTORIES ${HEADER_DIRECTORY})
endforeach()
list(REMOVE_DUPLICATES HEADER_DIRECTORIES)
target_include_directories(karamelo PUBLIC ${HEADER_DIRECTORIES})

target_link_libraries(karamelo PUBLIC MPI::MPI_CXX Kokkos::kokkos gzstream)

add_subdirectory(style)
add_subdirectory(version)
