cmake_minimum_required(VERSION 3.16)
project(devQP)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_BUILD_TYPE Debug)
#SET(CMAKE_CXX_FLAGS  "-fPIC")

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
find_package(MATIO)
find_package(fmt)
find_package(MKL REQUIRED)
#find_package(MPI REQUIRED)
#find_package(LAPACK REQUIRED)

#set(BLA_VENDOR Intel10_64lp)
# -- IPOPT Solver
#set(IPOPT_HOME /home/detu/CoinIpopt)
set(MUMPS_HOME /home/detu/ThirdParty-Mumps/Debug)
set(MUMPS_INCLUDE_DIR ${MUMPS_HOME}/include/coin-or/mumps)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${MUMPS_INCLUDE_DIR})
#include_directories(${MATLAB_INCLUDE_DIRS})
include_directories(${MATIO_INCLUDE_DIRS})
include_directories(${fmt_INCLUDE_DIR})
include_directories(${MKL_INCLUDE_DIR})
#include_directories(${MPI_INCLUDE_PATH})

# CoinIpopt
#include_directories(
#    ${IPOPT_HOME}/${CMAKE_BUILD_TYPE}/include/coin
#    ${IPOPT_HOME}/${CMAKE_BUILD_TYPE}/include/coin/ThirdParty
#    ${IPOPT_HOME}/${CMAKE_BUILD_TYPE}/include/extra
#)
#link_directories(
#    ${IPOPT_HOME}/${CMAKE_BUILD_TYPE}/lib
#)

message("${CMAKE_MODULE_PATH}")
message("${EIGEN3_INCLUDE_DIR}")
message("${fmt_INCLUDE_DIR}")
message("${MUMPS_INCLUDE_DIR}")
#message("${MPI_INCLUDE_PATH}")

link_directories(
#        ${MUMPS_HOME}/libseq
        ${MUMPS_HOME}/lib
)
#set(MUMPS_LIBRARIES libmpiseq.a libdmumps.a libmumps_common.a)
#set(MUMPS_LIBRARIES libcmumps.a libdmumps.a libmumps_common.a libpord.a libsmumps.a libzmumps.a)
#set(MUMPS_LIBRARIES libmpiseq.a libdmumps.a libmumps_common.a)
#set(EIG_MEXES almbound)
#set(EIG_MEXES qpAL)

#FOREACH(t ${EIG_MEXES})
#    matlab_add_mex(NAME ${t} SRC ${t}.cpp	LINK_TO ${MATLAB_LIBRARIES} ${MUMPS_LIBRARIES} ${fmt_LIBRARY} ${MKL_LIBRARY})
    #matlab_add_mex(NAME ${t} SRC ${t}.cpp	LINK_TO ${MATLAB_LIBRARIES} ${fmt_LIBRARY} libcoinmumps.a libcoinlapack.a libcoinblas.a libcoinmetis.a gfortran)
#ENDFOREACH(t ${EIG_MEXES})

#add_executable(devQP function.h main.cpp)
#include_directories(${PROJECT_SOURCE_DIR})
#add_executable(mextest mextest.cpp)
#add_executable(simple simple.cc)
add_executable(main main.cpp)
#target_link_libraries(main ${MATIO_LIBRARIES} ${fmt_LIBRARY} ${MKL_LIBRARIES} ${MUMPS_LIBRARIES} ${MPI_LIBRARIES} ${LAPACK_LIBRARIES})
#target_link_libraries(main ${MATIO_LIBRARIES} ${fmt_LIBRARY} ${MKL_LIBRARIES} ${MUMPS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(main ${MATIO_LIBRARIES} ${fmt_LIBRARY} ${MKL_LIBRARIES} libcoinmumps.so.2.1.0)
#target_link_libraries(main libmpiseq.a)
#target_link_libraries(main libdmumps.a)
#target_link_libraries(main libcoinmumps.a)
#target_link_libraries(main libcoinlapack.a)
#target_link_libraries(main libcoinblas.a)
#target_link_libraries(main libcoinmetis.a)
#target_link_libraries(main gfortran)
