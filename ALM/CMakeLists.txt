cmake_minimum_required(VERSION 3.16)
project(devQP)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED)
find_package(MUMPS REQUIRED)
#find_package(MATLAB COMPONENTS MX_LIBRARY)
find_package(MATIO REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})
#include_directories(${MUMPS_INCLUDE_DIR})
include_directories(${MATLAB_INCLUDE_DIRS})
include_directories(${MATIO_INCLUDE_DIRS})

message("${CMAKE_MODULE_PATH}")
message("${EIGEN3_INCLUDE_DIR}")

#set(EIG_MEXES almbound)

#FOREACH(t ${EIG_MEXES})
#    matlab_add_mex(NAME ${t} SRC ${t}.cpp	LINK_TO ${MATLAB_LIBRARIES} ${MUMPS_LIBRARIES})
#ENDFOREACH(t ${EIG_MEXES})

#add_executable(devQP function.h main.cpp)
#include_directories(${PROJECT_SOURCE_DIR})
#add_executable(mextest mextest.cpp)
#add_executable(simple simple.cc)
add_executable(main main.cpp)
target_link_libraries (main ${MATIO_LIBRARIES} ${MUMPS_LIBRARIES})
