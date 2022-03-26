cmake_minimum_required(VERSION 3.16)
project(devQP)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#set(CMAKE_CXX_FLAGS_DEBUG "-fno-stack-protector")

find_package(Eigen3 REQUIRED)
find_package(MATLAB COMPONENTS MX_LIBRARY)
find_package(MATIO)
find_package(fmt)
find_package(MKL REQUIRED)
find_package(MUMPS REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${MUMPS_INCLUDE_DIRS})
include_directories(${MATLAB_INCLUDE_DIRS})
include_directories(${MATIO_INCLUDE_DIRS})
include_directories(${fmt_INCLUDE_DIR})
include_directories(${MKL_INCLUDE_DIRS})

message("${MUMPS_LIBRARIES}")
message("${MKL_LIBRARIES}")

set(EIG_MEXES qpAL)

FOREACH(t ${EIG_MEXES})
  matlab_add_mex(NAME ${t} SRC ${t}.cpp LINK_TO ${MATLAB_LIBRARIES} ${MUMPS_LIBRARIES} ${fmt_LIBRARY} ${MKL_LIBRARIES})
ENDFOREACH(t ${EIG_MEXES})

add_executable(main main.cpp)
target_link_libraries(main ${MATIO_LIBRARIES} ${fmt_LIBRARY} ${MKL_LIBRARIES} ${MUMPS_LIBRARIES})

