add_library(gzstream STATIC gzstream.C)

target_include_directories(gzstream PUBLIC ./)

if(WIN32)
  add_library(zlib STATIC IMPORTED)
  set_target_properties(zlib PROPERTIES
    IMPORTED_LOCATION "C:/Program Files/Zlib/lib/zlib/zlibstatic.lib")
  target_include_directories(zlib INTERFACE "C:/Program Files/Zlib/include/zlib")

  target_link_libraries(gzstream PUBLIC zlib)
else()
  target_link_libraries(gzstream PUBLIC -lz)
endif()