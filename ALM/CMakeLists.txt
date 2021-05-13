cmake_minimum_required(VERSION 3.16)
project(devQP)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
set(CMAKE_CXX_STANDARD 17)

## Defult build type: Debug ##
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()
## Debug definition: ##
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_definitions("-DDEBUG")
endif()

find_package(Eigen3 REQUIRED)
#find_package(MUMPS REQUIRED)
#find_package(MATLAB COMPONENTS MX_LIBRARY)
find_package(MATIO REQUIRED)

# -- IPOPT Solver
set(IPOPT_HOME /home/detu/CoinIpopt)

include_directories(${EIGEN3_INCLUDE_DIR})
#include_directories(${MUMPS_INCLUDE_DIR})
include_directories(${MATLAB_INCLUDE_DIRS})
include_directories(${MATIO_INCLUDE_DIRS})

# CoinIpopt
include_directories(
    ${IPOPT_HOME}/${CMAKE_BUILD_TYPE}/include/coin
    ${IPOPT_HOME}/${CMAKE_BUILD_TYPE}/include/coin/ThirdParty
    ${IPOPT_HOME}/${CMAKE_BUILD_TYPE}/include/extra
)
link_directories(
    ${IPOPT_HOME}/${CMAKE_BUILD_TYPE}/lib
)

message("${CMAKE_MODULE_PATH}")
message("${EIGEN3_INCLUDE_DIR}")

#set(MUMPS_INCLUDE_DIR /usr/local/include/coin-or/mumps)
#set(EIG_MEXES almbound)

#FOREACH(t ${EIG_MEXES})
#    matlab_add_mex(NAME ${t} SRC ${t}.cpp	LINK_TO ${MATLAB_LIBRARIES} ${MUMPS_LIBRARIES})
#ENDFOREACH(t ${EIG_MEXES})

#add_executable(devQP function.h main.cpp)
#include_directories(${PROJECT_SOURCE_DIR})
#add_executable(mextest mextest.cpp)
#add_executable(simple simple.cc)
add_executable(main main.cpp)
target_link_libraries(main ${MATIO_LIBRARIES})
target_link_libraries(main libcoinmumps.a)
target_link_libraries(main libcoinlapack.a)
target_link_libraries(main libcoinblas.a)
target_link_libraries(main libcoinmetis.a)
target_link_libraries(main gfortran)
