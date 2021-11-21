add_executable(open_acc_demo open_acc_demo.cpp)

target_include_directories(open_acc_demo PUBLIC ./
                                                ./linalg)

find_package(OpenACC REQUIRED)
target_link_libraries(open_acc_demo PUBLIC OpenACC::OpenACC_CXX graphics)
target_compile_options(open_acc_demo PUBLIC -DUSE_COUNTING_ITERATOR
                                            -acc
                                            -ta=tesla:managed,lineinfo
                                            -Minfo=accel
                                            -lnvToolsExt)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lnvToolsExt -DUSE_COUNTING_ITERATOR -acc -ta=tesla:managed,lineinfo -Minfo=accel -lnvToolsExt -L/apps/nvidia-hpc-sdk/21.3/Linux_x86_64/21.3/cuda/11.2/lib64")