
# Copyright (C) 2016 iCub Facility - Istituto Italiano di Tecnologia
# Authors: Tariq Abuhashim, Nicolo' Genesio
# E-mail: t.abuhashim@gmail.com, nicogene@hotmail.it
# Date: Nov 2016
# Acknowledgement: This research has received funding from the European Union’s 
# Seventh Framework Programme for research, technological development and demonstration 
# under grant agreement No. 611909(KoroiBot).
# License: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT

cmake_minimum_required(VERSION 3.2.2)

set(PROJECT_OPTIMISE optimise)
PROJECT(${PROJECT_OPTIMISE})

# Make CMake aware of the cmake folder for local FindXXX scripts,
# append rather than set in case the user has passed their own
# additional paths via -D.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

message("----------${PROJECT_OPTIMISE}----------")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

########### Set postfixes for generated libraries based on buildtype.
set(CMAKE_RELEASE_POSTFIX "-release")
set(CMAKE_DEBUG_POSTFIX "-debug")

########### Flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")

########### Options
OPTION( BUILD_SHARED_LIBS 	"Enable shared libraries." 							OFF )
OPTION( YARP 				"Enable use of YARP." 								ON )
OPTION( ICUB 				"Enable use of ICUB." 								ON )
OPTION( OpenCV 				"Enable use of OpenCV." 							ON )
OPTION( LAPACK 				"Enable use of LAPACK." 							ON )
OPTION( EIGENSPARSE			"Enable Eigen as a sparse linear algebra library"	ON )
OPTION( SUITESPARSE 		"Enable SuiteSparse." 								ON )
OPTION( CXSparse 			"Enable CXSparse." 									ON )

########### YARP.
if(YARP)
	find_package( YARP REQUIRED )
	if(YARP_FOUND)
		message("-- Found YARP library: ${YARP_LIBRARIES}")
        include_directories(${YARP_INCLUDE_DIRS})
		SET(PROJECT_LIBS ${PROJECT_LIBS} ${YARP_LIBRARIES})
	else()
    	message("-- Did not find YARP library.")
	endif()
endif()

########### ICUB.
if(ICUB)
	find_package( ICUB REQUIRED )
	if(ICUB_FOUND)
		message("-- Found ICUB library: ${ICUB_LIBRARIES}")
		include_directories(${ICUB_INCLUDE_DIRS})
		SET(PROJECT_LIBS ${PROJECT_LIBS} ${ICUB_LIBRARIES})
	else()
		message("-- Did not find ICUB library.")
	endif()
endif()

########### OpenCV.
if(OpenCV)
	find_package( OpenCV REQUIRED )
	if(OpenCV_FOUND)
		message("-- Found OPENCV library: ${OpenCV_LIBS}")
		SET(PROJECT_LIBS ${PROJECT_LIBS} ${OpenCV_LIBS})
	else()
        message("-- Did not find OPENCV library.")
	endif()
endif()

########### LAPACK (& BLAS).
if(LAPACK)
	find_package(LAPACK REQUIRED )
	if(LAPACK_FOUND)
		message("-- Found LAPACK library: ${LAPACK_LIBRARIES}")
		SET(PROJECT_LIBS ${PROJECT_LIBS} ${LAPACK_LIBRARIES})
	else()
		message("-- Did not find LAPACK library, disabling LAPACK support.")
	endif()
	find_package(BLAS REQUIRED )
	if(BLAS_FOUND)
		message("-- Found BLAS library: ${BLAS_LIBRARIES}")
		SET(PROJECT_LIBS ${PROJECT_LIBS} ${BLAS_LIBRARIES})
	else()
		message("-- Did not find BLAS library, disabling LAPACK support.")
	endif()
	if(NOT (LAPACK_FOUND AND BLAS_FOUND))
	set(LAPACK OFF)
	#list(APPEND PWG_COMPILE_OPTIONS PWG_NO_LAPACK)
	endif()
else()
	message("-- Building without LAPACK.")
	#list(APPEND PWG_COMPILE_OPTIONS PWG_NO_LAPACK)
endif()

########### Eigen.
find_package(Eigen REQUIRED)
if(EIGEN_FOUND)
	if(EIGEN_VERSION VERSION_LESS 3.1.0)
		message(FATAL_ERROR "-- PwgOptimiser requires Eigen version >= 3.1.0 in order "
		"that Eigen/SparseCore be available, detected version of Eigen is: "
      "${EIGEN_VERSION}")
	endif()

	message("-- Found Eigen version ${EIGEN_VERSION}: ${EIGEN_INCLUDE_DIRS}")
	# Ensure that only MPL2 licensed code is part of the default build.
	
	if(EIGENSPARSE)
		list(APPEND PWG_COMPILE_OPTIONS PWG_USE_EIGEN_SPARSE)
		#message("")
		#message("   ===============================================================")
		#message("   Enabling the use of Eigen as a sparse linear algebra library ")
		#message("   for solving the nonlinear least squares problems.")
		#message("   Enabling this option results in an LGPL licensed version of ")
		#message("   solver as the Simplicial Cholesky factorization in Eigen")
		#message("   is licensed under the LGPL. ")
		#message("   ===============================================================")
		#message("")
		if(EIGEN_VERSION VERSION_LESS 3.2.2)
			message("")
			message("   ===============================================================")
			message("   WARNING:")
			message("")
			message("   Your version of Eigen is older than version 3.2.2.")
			message("   The performance of SPARSE_NORMAL_CHOLESKY and SPARSE_SCHUR")
			message("   linear solvers will suffer. ")
			message("   ===============================================================")
			message("")
		endif()
	else()
		message("")
		message("   ===============================================================")
		message("   Disabling the use of Eigen as a sparse linear algebra library.")
		#message("   This does not affect the covariance estimation algorithm ")
		#message("   which can still use the EIGEN_SPARSE_QR algorithm.")
		message("   ===============================================================")
		message("")
		add_definitions(-DEIGEN_MPL2_ONLY)
	endif()
	
endif()
if(EIGEN_FOUND)
	include_directories(${EIGEN_INCLUDE_DIRS})
endif()

########### SuiteSparse.
if(SUITESPARSE AND NOT LAPACK)
	# If user has disabled LAPACK, but left SUITESPARSE ON, turn it OFF,
	# LAPACK controls whether Ceres will be linked, directly or indirectly
	# via SuiteSparse to LAPACK.
	message("-- Disabling SuiteSparse as use of LAPACK has been disabled, "
	"turn ON LAPACK to enable (optional) building with SuiteSparse.")
	set(SUITESPARSE OFF)
endif()
if(SUITESPARSE)
	# By default, if SuiteSparse and all dependencies are found, Ceres is
	# built with SuiteSparse support.

	# Check for SuiteSparse and dependencies.
	find_package(SuiteSparse)
	if(SUITESPARSE_FOUND)
	
		# On Ubuntu the system install of SuiteSparse (v3.4.0) up to at least
		# Ubuntu 13.10 cannot be used to link shared libraries.
		if(BUILD_SHARED_LIBS AND SUITESPARSE_IS_BROKEN_SHARED_LINKING_UBUNTU_SYSTEM_VERSION)
			message(FATAL_ERROR "You are attempting to build Ceres as a shared "
        	"library on Ubuntu using a system package install of SuiteSparse "
        	"3.4.0. This package is broken and does not support the "
        	"construction of shared libraries (you can still build Ceres as "
        	"a static library).  If you wish to build a shared version of Ceres "
        	"you should uninstall the system install of SuiteSparse "
        	"(libsuitesparse-dev) and perform a source install of SuiteSparse "
        	"(we recommend that you use the latest version), "
        	"see http://ceres-solver.org/building.html for more information.")
		endif()
		
		# By default, if all of SuiteSparse's dependencies are found, PWG is
		# built with SuiteSparse support.
		message("-- Found SuiteSparse ${SUITESPARSE_VERSION}, "
 		"building with SuiteSparse.")
		#We can delete these checks
		if(AMD_FOUND)
			message("-- Found AMD: ${AMD_LIBRARY}")
		else()
			message("AMD not found")
		endif()
		if(CAMD_FOUND)
			message("-- Found CAMD: ${CAMD_LIBRARY}")
		else()
			message("CAMD not found:")
		endif()
		if(CCOLAMD_FOUND)
			message("-- Found CCOLAMD: ${CCOLAMD_LIBRARY}")
		else()
			message("CCOLAMD not found")
		endif()
		if(COLAMD_FOUND)
			message("-- Found COLAMD: ${COLAMD_LIBRARY}")
		else()
			message("COLAMD not found")
		endif()
		if(SUITESPARSE_CONFIG_FOUND)
			message("-- Found SUITESPARSE_CONFIG: ${SUITESPARSE_CONFIG_LIBRARY}")
		else()
			message("SUITESPARSE_CONFIG not found")
		endif()
		if(METIS_FOUND)
			message("-- Found METIS: ${METIS_LIBRARY}")
		else()
			message("METIS not found")
		endif()	
	else()
	
		# Disable use of SuiteSparse if it cannot be found and continue.
		message("-- Did not find all SuiteSparse dependencies, disabling "
		"SuiteSparse support.")
		set(SUITESPARSE OFF)
		#list(APPEND PWG_COMPILE_OPTIONS PWG_NO_SUITESPARSE)
		
	endif()
else()
	message("-- Building without SuiteSparse.")
	#list(APPEND PWG_COMPILE_OPTIONS PWG_NO_SUITESPARSE)
endif()
if(SUITESPARSE)
	include_directories(${SUITESPARSE_INCLUDE_DIRS})
endif()

########### CXSparse.
if(CXSPARSE)
	# Don't search with REQUIRED as we can continue without CXSparse.
	find_package(CXSparse)
	if(CXSPARSE_FOUND)
		# By default, if CXSparse and all dependencies are found, Ceres is
		# built with CXSparse support.
		message("-- Found CXSparse version: ${CXSPARSE_VERSION}, building with CXSparse.")
	else()
		# Disable use of CXSparse if it cannot be found and continue.
		message("-- Did not find CXSparse, Building without CXSparse.")
		set(CXSPARSE OFF)
		#list(APPEND PWG_COMPILE_OPTIONS PWG_NO_CXSPARSE)
	endif()
else()
	message("-- Building without CXSparse.")
	#list(APPEND PWG_COMPILE_OPTIONS PWG_NO_CXSPARSE)
	# Mark as advanced (remove from default GUI view) the CXSparse search
	# variables in case user enabled CXSPARSE, FindCXSparse did not find it, so
	# made search variables visible in GUI for user to set, but then user disables
	# CXSPARSE instead of setting them.
	mark_as_advanced(FORCE CXSPARSE_INCLUDE_DIR CXSPARSE_LIBRARY)
endif()
if(CXSPARSE)
	include_directories(${CXSPARSE_INCLUDE_DIRS})
endif()
# Ensure that the user understands they have disabled all sparse libraries.
if(NOT SUITESPARSE AND NOT CXSPARSE AND NOT EIGENSPARSE)
	message("   ===============================================================")
	message("   Compiling without any sparse library: SuiteSparse, CXSparse ")
	message("   & Eigen (Sparse) are all disabled or unavailable.  No sparse ")
	message("   linear solvers (SPARSE_NORMAL_CHOLESKY & SPARSE_SCHUR)")
	message("   will be available when Ceres is used.")
	message("   ===============================================================")
endif()

###########cholmod
find_package(Cholmod)
if(CHOLMOD_FOUND)
	message("-- Found Cholmod: ${CHOLMOD_LIBRARIES}")
	SET(PROJECT_LIBS ${PROJECT_LIBS} ${PARMETIS_LIBRARIES})
else()
	message("Cholmod not found")
endif()
## metis 4.0 to be removed
########### metis-4.0
#set(METIS4_LIB ${suitesparse}/metis-4.0)
find_library(METIS4_LIB
			 PATHS ../SuiteSparse/metis-4.0
			 NAMES libmetis.a)
if(EXISTS ${METIS4_LIB})
	message("-- METIS-4.0 libs: ${METIS4_LIB}")
endif()

########### parmetis
find_package(PARMETIS)
if(PARMETIS_FOUND)
	message("-- Found PARMETIS: ${PARMETIS_LIBRARIES}")
	SET(PROJECT_LIBS ${PROJECT_LIBS} ${PARMETIS_LIBRARIES})
else()
  message("PARMETIS not found")
endif()

########### GSL.
find_package(GSL)
if(GSL_FOUND)
	message("-- Found GSL: ${GSL_LIBRARIES}")
  SET(PROJECT_LIBS ${PROJECT_LIBS} ${GSL_LIBRARIES})
else()
  message("GSL not found")
endif()

#-lglfw -O3 -Wall -lCg -lCgGL -lpthread -lglut -lGL -lGLU -lX11 -lXxf86vm -lXrandr -lXi -lopencv_core -lopencv_flann -lopencv_highgui
#-lgfortran -lceres -lcholmod -lcolamd -lamd -lcamd -lcxsparse -lsuitesparseconfig -lgomp -lglog -lpthread

########### MATLAB
set(MATLAB_INCLUDE_DIRS ../matlab/include)
message("-- Using MATLAB: ${MATLAB_INCLUDE_DIRS}")
include_directories(${MATLAB_INCLUDE_DIRS})
#find_package( matlab )
#find_path(MATLAB_INCLUDE_DIR NAMES extern PATHS /usr/local/MATLAB/R2016a/extern/include)
#find_library(MATLAB_LIB NAMES matlab)
#if(EXISTS ${MATLAB_LIB})
#	message("-- Found MATLAB: ${MATLAB_INCLUDE_DIR}")
#	include_directories(${MATLAB_INCLUDE_DIR})
#endif(EXISTS ${MATLAB_LIB})
#set(MATLAB_EXTERN /usr/local/MATLAB/R2016a/extern/include)
#set(MATLAB_EXTERN ./include/matlab)
#message("-- Using MATLAB_EXTERN: ${MATLAB_EXTERN}")
#include_directories(${MATLAB_EXTERN})

########### vlfeat.
find_package(VLFeat)
if(VLFEAT_FOUND)
	message("-- Using VLFeat: ${VLFEAT_INCLUDE_DIRS}")
	include_directories(${VLFEAT_INCLUDE_DIRS})
	SET(PROJECT_LIBS ${PROJECT_LIBS} ${VLFEAT_LIBRARIES})
else()
	message("-- VLFeat not found, set the VLFEAT_ROOT to the root folder of VLFeat ")
endif()

########### PwgOptimiser sources
#set(Executable_SOURCES ./src/main.cpp ./src/PwgOptimiser.cpp ./src/RecoverMoments.cpp)
#set(Executable_HEADERS ./include/PwgOptimiser.h ./include/RecoverMoments.h)
set(SOURCE_DIR ./src)
set(INCLUDE_DIR ./include)
include_directories(${SOURCE_DIR} ${INCLUDE_DIR})
file(GLOB Executable_SOURCES ${SOURCE_DIR}/*.cpp)
file(GLOB Executable_HEADERS ${INCLUDE_DIR}/*.h)
add_executable(${PROJECT_OPTIMISE} ${Executable_SOURCES})

target_link_libraries(${PROJECT_OPTIMISE} ${YARP_LIBRARIES} ${ICUB_LIBRARIES} ${OpenCV_LIBS}
											${CHOLMOD_LIBRARIES} ${COLAMD_LIBRARY} ${CCOLAMD_LIBRARY} ${AMD_LIBRARY}
										  ${CAMD_LIBRARY} ${SUITESPARSE_CONFIG_LIBRARY} ${LAPACK_LIBRARIES}
											${VLFEAT_LIBRARIES} ${${METIS_LIBRARY}} ${METIS4_LIB} ${PARMETIS_LIBRARIES}
										  ${GSL_LIB} ${BLAS_LIBRARIES})
#target_link_libraries(${PROJECT_OPTIMISE} ${PROJECT_LIBS})

###########
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
  target_compile_features(${PROJECT_OPTIMISE} PRIVATE cxx_range_for)
endif()
