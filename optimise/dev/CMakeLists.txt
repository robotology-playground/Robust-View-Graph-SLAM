cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(ceres_test)
set(CMAKE_VERBOSE_MAKEFILE 1)
#set(LIBS_TYPE STATIC)

#set(CERES_DIR /home/tariq/Dev/ceres-solver/build)
find_package(Ceres REQUIRED)

set(EIGEN_DIR /home/tariq/Dev/eigen/build)

set(GLOG_DIR /home/tariq/Dev/glog/build)

set(GFLAGS_DIR /home/tariq/Dev/gflags/build)

link_directories(${GLOG_DIR}/lib ${GFLAGS_DIR}/lib)

include_directories(${CERES_INCLUDE_DIRS} ${GLOG_DIR}/include ${GFLAGS_DIR}/include ${EIGEN_DIR}/include/eigen3)

#include(CheckCXXCompilerFlag)
#check_cxx_compiler_flag("-std=c++11" COMPILER_SUPPORTS_CXX11)

add_executable(test_ceres test.cpp)
target_link_libraries(test_ceres ${CERES_LIBRARIES} ${GLOG_DIR}/lib/libglog.so)
