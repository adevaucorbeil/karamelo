add_executable(kokkos_demo kokkos_demo.cpp)

target_include_directories(kokkos_demo PUBLIC ./
                                              ./linalg)

target_link_libraries(kokkos_demo PUBLIC Kokkos::kokkos graphics)