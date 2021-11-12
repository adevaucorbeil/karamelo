add_executable(open_acc_demo open_acc_demo.cpp)

target_include_directories(open_acc_demo PUBLIC ./
                                                ./linalg)

find_package(OpenACC REQUIRED)
target_link_libraries(open_acc_demo PUBLIC OpenACC::OpenACC_CXX graphics)