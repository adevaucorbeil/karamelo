add_executable(kokkos_demo kokkos_demo.cpp)

target_include_directories(kokkos_demo PUBLIC ./
                                              ./linalg)

target_link_libraries(kokkos_demo PUBLIC Kokkos::kokkos graphics)

add_executable(rdf rdf.cpp)
target_include_directories(rdf PUBLIC ./)
target_link_libraries(rdf PUBLIC Kokkos::kokkos)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lnvToolsExt")