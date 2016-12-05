cmake_minimum_required(VERSION 3.2.2)

SET(PROJECTNAME vgSLAM_build_graph)
PROJECT(${PROJECTNAME})

# Make CMake aware of the cmake folder for local FindXXX scripts,
# append rather than set in case the user has passed their own
# additional paths via -D.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../cmake")

message(${CMAKE_MODULE_PATH})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Set postfixes for generated libraries based on buildtype.
set(CMAKE_RELEASE_POSTFIX "")
set(CMAKE_DEBUG_POSTFIX "-debug")

# Flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")

# Options
option(YARP "Enable use of YARP." ON)
option(OpenCV "Enable use of OpenCV." ON)
option(ICUB "Enable use of iCub." ON)
# Enable the use of Eigen as a sparse linear algebra library for
# solving the nonlinear least squares problems. 
# Simplicial Cholesky factorization in Eigen is licensed under the LGPL.
# This does not affect the covariance estimation algorithm, as it
# depends on the sparse QR factorization algorithm, which is licensed
# under the MPL.
OPTION(EIGENSPARSE "Enable Eigen as a sparse linear algebra library" ON)

# Unset a variable, cache variable, or environment variable.
unset(PWG_COMPILE_OPTIONS)
unset(PWG_LIBRARIES)

# YARP.
find_package( YARP REQUIRED )
if (YARP_FOUND)
	message("-- Found YARP library: ${YARP_LIBRARIES}")
	include_directories(${YARP_INCLUDE_DIRS})
else (YARP_FOUND)
	message("-- Did not find YARP library.")
endif (YARP_FOUND)

# iCub
find_package( ICUB REQUIRED )
if (ICUB_FOUND)
	message("-- Found iCub library: ${ICUB_LIBRARIES}")
include_directories(${ICUB_INCLUDE_DIRS})
else (ICUB_FOUND)
	message("-- Did not find YARP library.")
endif (ICUB_FOUND)

# OpenCV.
find_package(OpenCV REQUIRED)
if (OpenCV_FOUND)
	message("-- Found OPENCV library: ${OpenCV_LIBS}")
else (OpenCV_FOUND)
    	message("-- Did not find OPENCV library.")
endif (OpenCV_FOUND)


set(SOURCE_DIR ./src)
set(INCLUDE_DIR ./include)
include_directories(${SOURCE_DIR} ${INCLUDE_DIR} ../src ../include)
file(GLOB Executable_SOURCES ${SOURCE_DIR}/*.cpp ../src/featureselector.cpp)
file(GLOB Executable_HEADERS ${INCLUDE_DIR}/*.h  ../include/featureselector.h)
add_executable(${PROJECTNAME} ${Executable_SOURCES} ${Executable_HEADERS})

target_link_libraries(${PROJECTNAME} ${YARP_LIBRARIES} ${ICUB_LIBRARIES}  ${OpenCV_LIBS})

if(${CMAKE_VERSION} VERSION_LESS 3.1)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-std=c++11" COMPILER_SUPPORTS_CXX11)
  check_cxx_compiler_flag("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
  if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
  else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
  endif()
else()
  target_compile_features(${PROJECTNAME} PRIVATE cxx_range_for)
endif()
