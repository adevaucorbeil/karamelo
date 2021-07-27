add_library(gzstream SHARED gzstream.C)

target_include_directories(gzstream PUBLIC ./)

target_link_libraries(gzstream PUBLIC -lz)