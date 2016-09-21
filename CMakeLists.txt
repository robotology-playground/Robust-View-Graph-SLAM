# Copyright: 2016 iCub Facility, Istituto Italiano di Tecnologia
# Author: Tariq Abuhashim
# CopyPolicy: Released under the terms of the GNU GPL v2.0.
# 

cmake_minimum_required(VERSION 3.2.2)

SET(PROJECTNAME vgSLAM)
PROJECT(${PROJECTNAME})

# Make CMake aware of the cmake folder for local FindXXX scripts,
# append rather than set in case the user has passed their own
# additional paths via -D.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

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
option(ICUB "Enable use of ICUB." ON)
option(OpenCV "Enable use of OpenCV." ON)
option(LAPACK "Enable use of LAPACK." ON)
option(SUITESPARSE "Enable SuiteSparse." ON)
option(CXSPARSE "Enable CXSparse." ON)
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
if (YARP)
	find_package( YARP REQUIRED )
	if (YARP_FOUND)
		message("-- Found YARP library: ${YARP_LIBRARIES}")
    	include_directories(${YARP_INCLUDE_DIR})
	else (YARP_FOUND)
    	message("-- Did not find YARP library.")
	endif (YARP_FOUND)
endif (YARP)

# ICUB.
if (ICUB)
	find_package( ICUB REQUIRED )
	if (ICUB_FOUND)
		message("-- Found ICUB library: ${ICUB_LIBRARIES}")
		include_directories(${ICUB_INCLUDE_DIRS})
	else (ICUB_FOUND)
		message("-- Did not find ICUB library.")
	endif (ICUB_FOUND)
endif (ICUB)

# OpenCV.
if (OpenCV)
	find_package(OpenCV REQUIRED)
	if (OpenCV_FOUND)
		message("-- Found OPENCV library: ${OpenCV_INCLUDE_DIRS}")
    	include_directories(${OpenCV_INCLUDE_DIRS})
	else (OpenCV_FOUND)
    	message("-- Did not find OPENCV library.")
	endif (OpenCV_FOUND)
endif (OpenCV)


# LAPACK (& BLAS).
if (LAPACK)
	find_package(LAPACK REQUIRED )
	if (LAPACK_FOUND)
		message("-- Found LAPACK library: ${LAPACK_LIBRARIES}")
	else (LAPACK_FOUND)
		message("-- Did not find LAPACK library, disabling LAPACK support.")
	endif (LAPACK_FOUND)
	find_package(BLAS REQUIRED )
	if (BLAS_FOUND)
		message("-- Found BLAS library: ${BLAS_LIBRARIES}")
	else (BLAS_FOUND)
		message("-- Did not find BLAS library, disabling LAPACK support.")
	endif (BLAS_FOUND)
	if (NOT (LAPACK_FOUND AND BLAS_FOUND))
	set(LAPACK OFF)
	list(APPEND PWG_COMPILE_OPTIONS PWG_NO_LAPACK)
	endif (NOT (LAPACK_FOUND AND BLAS_FOUND))
else (LAPACK)
	message("-- Building without LAPACK.")
	list(APPEND PWG_COMPILE_OPTIONS PWG_NO_LAPACK)
endif (LAPACK)

# Eigen.
#find_path(Eigen_INCLUDE_DIR NAMES signature_of_eigen3_matrix_library PATHS $ENV{EIGEN_ROOT} /usr/include/eigen3 /usr/local/include/eigen3)
#if(NOT EXISTS ${Eigen_INCLUDE_DIR})
#    message( "could not find Eigen - should be the path containing the Eigen 3 directory and /Eigen/Eigen, /Eigen/Core, etc..." )
#else(NOT EXISTS ${Eigen_INCLUDE_DIR})
#    message("-- Found Eigen: ${Eigen_INCLUDE_DIR}")
#    include_directories(${Eigen_INCLUDE_DIR})
#endif(NOT EXISTS ${Eigen_INCLUDE_DIR})

# Eigen.
find_package(Eigen REQUIRED)
if (EIGEN_FOUND)
	if (EIGEN_VERSION VERSION_LESS 3.1.0)
		message(FATAL_ERROR "-- PwgOptimiser requires Eigen version >= 3.1.0 in order "
		"that Eigen/SparseCore be available, detected version of Eigen is: "
      "${EIGEN_VERSION}")
	endif (EIGEN_VERSION VERSION_LESS 3.1.0)

	message("-- Found Eigen version ${EIGEN_VERSION}: ${EIGEN_INCLUDE_DIRS}")
	# Ensure that only MPL2 licensed code is part of the default build.
	message("")
	message("   ===============================================================")
	if (EIGENSPARSE)
		list(APPEND PWG_COMPILE_OPTIONS PWG_USE_EIGEN_SPARSE)
		message("   Enabling the use of Eigen as a sparse linear algebra library ")
		message("   for solving the nonlinear least squares problems. Enabling ")
		message("   this option results in an LGPL licensed version of ")
		message("   solver as the Simplicial Cholesky factorization in Eigen")
		message("   is licensed under the LGPL. ")
		if (EIGEN_VERSION VERSION_LESS 3.2.2)
			message("   WARNING:")
			message("")
			message("   Your version of Eigen is older than version 3.2.2.")
			message("   The performance of SPARSE_NORMAL_CHOLESKY and SPARSE_SCHUR")
			message("   linear solvers will suffer. ")
		endif (EIGEN_VERSION VERSION_LESS 3.2.2)
	else (EIGENSPARSE)
		message("   Disabling the use of Eigen as a sparse linear algebra library.")
		message("   This does not affect the covariance estimation algorithm ")
		message("   which can still use the EIGEN_SPARSE_QR algorithm.")
		add_definitions(-DEIGEN_MPL2_ONLY)
	endif (EIGENSPARSE)
	message("   ===============================================================")
	message("")
endif (EIGEN_FOUND)
if (EIGEN_FOUND)
	include_directories(${EIGEN_INCLUDE_DIRS})
endif (EIGEN_FOUND)

# SuiteSparse.
if (SUITESPARSE AND NOT LAPACK)
	# If user has disabled LAPACK, but left SUITESPARSE ON, turn it OFF,
	# LAPACK controls whether Ceres will be linked, directly or indirectly
	# via SuiteSparse to LAPACK.
	message("-- Disabling SuiteSparse as use of LAPACK has been disabled, "
	"turn ON LAPACK to enable (optional) building with SuiteSparse.")
	set(SUITESPARSE OFF)
endif (SUITESPARSE AND NOT LAPACK)
if (SUITESPARSE)
	# By default, if SuiteSparse and all dependencies are found, Ceres is
	# built with SuiteSparse support.

	# Check for SuiteSparse and dependencies.
	find_package(SuiteSparse)
	if (SUITESPARSE_FOUND)
	
		# On Ubuntu the system install of SuiteSparse (v3.4.0) up to at least
		# Ubuntu 13.10 cannot be used to link shared libraries.
		if (BUILD_SHARED_LIBS AND SUITESPARSE_IS_BROKEN_SHARED_LINKING_UBUNTU_SYSTEM_VERSION)
			message(FATAL_ERROR "You are attempting to build Ceres as a shared "
        	"library on Ubuntu using a system package install of SuiteSparse "
        	"3.4.0. This package is broken and does not support the "
        	"construction of shared libraries (you can still build Ceres as "
        	"a static library).  If you wish to build a shared version of Ceres "
        	"you should uninstall the system install of SuiteSparse "
        	"(libsuitesparse-dev) and perform a source install of SuiteSparse "
        	"(we recommend that you use the latest version), "
        	"see http://ceres-solver.org/building.html for more information.")
		endif (BUILD_SHARED_LIBS AND SUITESPARSE_IS_BROKEN_SHARED_LINKING_UBUNTU_SYSTEM_VERSION)
		
		# By default, if all of SuiteSparse's dependencies are found, PWG is
		# built with SuiteSparse support.
		message("-- Found SuiteSparse ${SUITESPARSE_VERSION}, "
 		"building with SuiteSparse.")
 		
	else (SUITESPARSE_FOUND)
	
		# Disable use of SuiteSparse if it cannot be found and continue.
		message("-- Did not find all SuiteSparse dependencies, disabling "
		"SuiteSparse support.")
		set(SUITESPARSE OFF)
		list(APPEND PWG_COMPILE_OPTIONS PWG_NO_SUITESPARSE)
		
	endif (SUITESPARSE_FOUND)
else (SUITESPARSE)
	message("-- Building without SuiteSparse.")
	list(APPEND PWG_COMPILE_OPTIONS PWG_NO_SUITESPARSE)
endif (SUITESPARSE)
if (SUITESPARSE)
	include_directories(${SUITESPARSE_INCLUDE_DIRS})
endif (SUITESPARSE)

# CXSparse.
if (CXSPARSE)
	# Don't search with REQUIRED as we can continue without CXSparse.
	find_package(CXSparse)
	if (CXSPARSE_FOUND)
		# By default, if CXSparse and all dependencies are found, Ceres is
		# built with CXSparse support.
		message("-- Found CXSparse version: ${CXSPARSE_VERSION}, "
      "building with CXSparse.")
	else (CXSPARSE_FOUND)
		# Disable use of CXSparse if it cannot be found and continue.
		message("-- Did not find CXSparse, Building without CXSparse.")
		set(CXSPARSE OFF)
		list(APPEND PWG_COMPILE_OPTIONS PWG_NO_CXSPARSE)
	endif (CXSPARSE_FOUND)
else (CXSPARSE)
	message("-- Building without CXSparse.")
	list(APPEND PWG_COMPILE_OPTIONS PWG_NO_CXSPARSE)
	# Mark as advanced (remove from default GUI view) the CXSparse search
	# variables in case user enabled CXSPARSE, FindCXSparse did not find it, so
	# made search variables visible in GUI for user to set, but then user disables
	# CXSPARSE instead of setting them.
	mark_as_advanced(FORCE CXSPARSE_INCLUDE_DIR CXSPARSE_LIBRARY)
endif (CXSPARSE)
if (CXSPARSE)
	include_directories(${CXSPARSE_INCLUDE_DIRS})
endif (CXSPARSE)

# Ensure that the user understands they have disabled all sparse libraries.
if (NOT SUITESPARSE AND NOT CXSPARSE AND NOT EIGENSPARSE)
	message("   ===============================================================")
	message("   Compiling without any sparse library: SuiteSparse, CXSparse ")
	message("   & Eigen (Sparse) are all disabled or unavailable.  No sparse ")
	message("   linear solvers (SPARSE_NORMAL_CHOLESKY & SPARSE_SCHUR)")
	message("   will be available when Ceres is used.")
	message("   ===============================================================")
endif(NOT SUITESPARSE AND NOT CXSPARSE AND NOT EIGENSPARSE)

# Include libraries
#cholmod
find_library(cholmod_LIB NAMES cholmod)
if (EXISTS ${cholmod_LIB})
	message("-- cholmod libs: ${cholmod_LIB}")
	#add_library(${cholmod_LIB})
endif(EXISTS ${cholmod_LIB})
# lapack
# suitesparseconfig
find_library(SSC_LIB NAMES suitesparseconfig)
if (EXISTS ${SSC_LIB})
	message("-- SuiteSparse_config libs: ${SSC_LIB}")
	#add_library(${SSC_LIB})
endif(EXISTS ${SSC_LIB})
# colamd
find_library(COLAMD_LIB NAMES colamd)
if (EXISTS ${COLAMD_LIB})
	message("-- Colamd libs: ${COLAMD_LIB}")
	#add_library(${COLAMD_LIB})
endif(EXISTS ${COLAMD_LIB})
# amd
find_library(AMD_LIB NAMES amd)
if (EXISTS ${AMD_LIB})
	message("-- AMD libs: ${AMD_LIB}")
endif(EXISTS ${AMD_LIB})

#-lglfw -O3 -Wall -lCg -lCgGL -lpthread -lglut -lGL -lGLU -lX11 -lXxf86vm -lXrandr -lXi -lopencv_core -lopencv_flann -lopencv_highgui
#-lgfortran -lceres -lcholmod -lcolamd -lamd -lcamd -lcxsparse -lsuitesparseconfig -lgomp -lglog -lpthread

# MATLAB
#find_package( matlab )
#find_path(MATLAB_INCLUDE_DIRS NAMES extern PATHS /usr/local/MATLAB/R2016a/extern/include)
find_path(MATLAB_INCLUDE_DIRS NAMES extern PATHS ./include/matlab)
if(NOT EXISTS ${MATLAB_INCLUDE_DIRS})
	message( "could not find MATLAB." )
else(NOT EXISTS ${MATLAB_INCLUDE_DIRS})
	message("-- Found MATLAB: ${MATLAB_INCLUDE_DIRS}")
	include_directories(${MATLAB_INCLUDE_DIRS})
endif(NOT EXISTS ${MATLAB_INCLUDE_DIRS})
#set(MATLAB_EXTERN /usr/local/MATLAB/R2016a/extern/include)
set(MATLAB_EXTERN ./include/matlab)
message("-- Using MATLAB_EXTERN: ${MATLAB_EXTERN}")
include_directories(${MATLAB_EXTERN})

# vlfeat.
set(VLFEAT_INCLUDE_DIR /home/ngenesio/Downloads/vlfeat-0.9.20)
message("-- Using VLFeat: ${VLFEAT_INCLUDE_DIR}")
include_directories(${VLFEAT_INCLUDE_DIR})
find_library(VLFEAT_LIB NAMES vl PATHS /home/ngenesio/Downloads/vlfeat-0.9.20/bin/glnxa64)
#set(VLFEAT_LIB /home/tabuhashim/Dev/vlfeat-0.9.20/bin/glnxa64)
if (EXISTS ${VLFEAT_LIB})
	message("-- VLFEAT libs: ${VLFEAT_LIB}")
endif(EXISTS ${VLFEAT_LIB})

# GSL.
#find_library(GSL_LIB NAMES gsl)
#if (EXISTS ${GSL_LIB})
#	message("-- GSL libs: ${GSL_LIB}")
#endif(EXISTS ${GSL_LIB})

# PwgOptimiser sources
#set(Executable_SOURCES ./src/main.cpp ./src/PwgOptimiser.cpp ./src/RecoverMoments.cpp)
#set(Executable_HEADERS ./include/PwgOptimiser.h ./include/RecoverMoments.h)
set(SOURCE_DIR ./src)
set(INCLUDE_DIR ./include)
include_directories(${SOURCE_DIR} ${INCLUDE_DIR})
file(GLOB Executable_SOURCES ${SOURCE_DIR}/*.cpp)
file(GLOB Executable_HEADERS ${INCLUDE_DIR}/*.h)
add_executable(${PROJECTNAME} ${Executable_SOURCES})

target_link_libraries(${PROJECTNAME} ${YARP_LIBRARIES} ${ICUB_LIBRARIES} ${OpenCV_LIBS} ${cholmod_LIB} ${COLAMD_LIB} ${AMD_LIB} ${SSC_LIB} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${VLFEAT_LIB})# ${GSL_LIB})

if(${CMAKE_VERSION} VERSION_LESS 3.1)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-std=c++11" COMPILER_SUPPORTS_CXX11)
  check_cxx_compiler_flag("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
  if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ")
  elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
  else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
  endif()
else()
  target_compile_features(${PROJECTNAME} PRIVATE cxx_range_for)
endif()
