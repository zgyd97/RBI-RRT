# cmake_minimum_required(VERSION 3.10)
# project(RBI-RRT)

# set(CMAKE_CXX_STANDARD 14)

# include_directories(
#     include
#     /usr/include/eigen3
#     /usr/include/ompl
# )

# # 指定源文件
# set(SOURCE_FILES
#     src/RBI-RRT.cpp
#     src/kdtree.c
# )

# # 创建可执行文件
# add_executable(RBI-RRT ${SOURCE_FILES})

# find_package(Eigen3 REQUIRED)

# target_link_libraries(RBI-RRT
#     ompl
#     Eigen3::Eigen
# )

# =================dirrt-connect===============
cmake_minimum_required(VERSION 3.10)
project(DIRRT_Connect)

set(CMAKE_CXX_STANDARD 14)

include_directories(
    include
    /usr/include/eigen3
    /usr/include/ompl
)

# 指定源文件
set(SOURCE_FILES
    src/DIRRT_Connect.cpp
    src/kdtree.c
)

# 创建可执行文件
add_executable(DIRRT_Connect ${SOURCE_FILES})

find_package(Eigen3 REQUIRED)

target_link_libraries(DIRRT_Connect
    ompl
    Eigen3::Eigen
)



