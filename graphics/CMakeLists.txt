add_library(graphics STATIC graphics.cpp)

target_include_directories(graphics PUBLIC ./)

target_link_libraries(graphics PUBLIC glad glfw)