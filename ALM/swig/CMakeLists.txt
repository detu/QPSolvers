cmake_minimum_required(VERSION 3.16)
project(devQP)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_POSITION_INDEPENDENT_CODE ON)


## Defult build type: Debug ##
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()
## Debug definition: ##
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_definitions("-DDEBUG")
endif()

find_package(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

find_package(Eigen3 REQUIRED)
find_package(MATLAB COMPONENTS MX_LIBRARY)
find_package(MKL REQUIRED)

set(MUMPS_HOME /home/detu/ThirdParty-Mumps/Debug)
set(MUMPS_INCLUDE_DIR ${MUMPS_HOME}/include/coin-or/mumps)

include_directories(${EIGEN3_INCLUDE_DIR})
#include_directories(${MUMPS_INCLUDE_DIRS})
include_directories(${MATLAB_INCLUDE_DIRS})
include_directories(${MKL_INCLUDE_DIRS})

message("${CMAKE_MODULE_PATH}")
message("${EIGEN3_INCLUDE_DIR}")
message("${fmt_INCLUDE_DIR}")
message("${MKL_INCLUDE_DIRS}")

#add_library(eigenMumps SHARED)

#link_directories(
#        ${MUMPS_HOME}/lib
#)
#set(EIG_MEXES almbound)
set(EIG_MEXES qpAL)

FOREACH(t ${EIG_MEXES})
  #matlab_add_mex(NAME ${t} SRC ${t}.cpp	LINK_TO ${MATLAB_LIBRARIES} libcoinmumps.so.2.1.0 ${fmt_LIBRARY} ${MKL_LIBRARY} gfortran)
  # matlab_add_mex(NAME ${t} SRC ${t}.cpp LINK_TO ${MATLAB_LIBRARIES} ${MUMPS_LIBRARIES} ${fmt_LIBRARY} ${MKL_LIBRARY})
ENDFOREACH(t ${EIG_MEXES})


add_executable(main main.cpp)
target_link_libraries(main ${MATIO_LIBRARIES} ${fmt_LIBRARY} ${MKL_LIBRARY} ${MUMPS_LIBRARIES} gfortran)

