
# Copyright (C) 2016 iCub Facility - Istituto Italiano di Tecnologia
# Authors: Tariq Abuhashim, Nicolo' Genesio
# E-mail: t.abuhashim@gmail.com, nicogene@hotmail.it
# Date: Nov 2016
# Acknowledgement: This research has received funding from the European Union’s 
# Seventh Framework Programme for research, technological development and demonstration 
# under grant agreement No. 611909(KoroiBot).
# License: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT

cmake_minimum_required(VERSION 3.2.2)

# ---[ Using cmake scripts and modules
# Make CMake aware of the cmake folder for local FindXXX scripts,
# append rather than set in case the user has passed their own
# additional paths via -D.
#SET(VLFEAT_ROOT /home/tariq/Dev/matlab)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

OPTION( BUILD_CAMERA_GRAPH	"Build camera graph"			OFF )
OPTION( BUILD_OPTIMISER		"Build optimisation"			OFF )
OPTION( BUILD_EXAMPLES		"Build examples"				OFF )
OPTION( BUILD_SHARED_LIBS	"Build as shared libraries."	OFF)
OPTION( OpenCV				"Enable use of OpenCV."			ON )
OPTION( VLFeat				"Enable use of VLFEAT"			ON )
OPTION( EIGENSPARSE			"Enable use of Eigen."			ON )

# ---[ OpenCV.
if(OpenCV)
	find_package(OpenCV REQUIRED)
	if (OpenCV_FOUND)
		message("-- Found OPENCV library: ${OpenCV_LIBS}")
	else (OpenCV_FOUND)
		message("-- Did not find OPENCV library.")
	endif (OpenCV_FOUND)
endif(OpenCV)

# ---[ vlfeat.
if(VLFeat)
	find_package(VLFeat REQUIRED)
	if(VLFEAT_FOUND)
		message("-- Using VLFeat: ${VLFEAT_INCLUDE_DIRS}")
		include_directories(${VLFEAT_INCLUDE_DIRS})
		SET(PROJECT_LIBS ${PROJECT_LIBS} ${VLFEAT_LIBRARIES})
	else(VLFEAT_FOUND)
		message("-- VLFeat not found, set the VLFEAT_ROOT to the root folder of VLFeat ")
	endif(VLFEAT_FOUND)
endif(VLFeat)

# ---[ Eigen.
if(EIGENSPARSE)
	find_package(Eigen REQUIRED)
	if(EIGEN_FOUND)
		if(EIGEN_VERSION VERSION_LESS 3.1.0)
			message(FATAL_ERROR "-- PwgOptimiser requires Eigen version >= 3.1.0 in order "
								"that Eigen/SparseCore be available, detected version of Eigen is: "
      							"${EIGEN_VERSION}")
		endif()
		#message("-- Found Eigen version ${EIGEN_VERSION}: ${EIGEN_INCLUDE_DIRS}")
		# Ensure that only MPL2 licensed code is part of the default build.
		list(APPEND PWG_COMPILE_OPTIONS PWG_USE_EIGEN_SPARSE)
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
		include_directories(${EIGEN_INCLUDE_DIRS})
	else()
		message("")
		message("   ===============================================================")
		message("   Disabling the use of Eigen as a sparse linear algebra library.")
		message("   ===============================================================")
		message("")
		add_definitions(-DEIGEN_MPL2_ONLY)
	endif()
endif()

# ---[ build camera graph.
if (BUILD_CAMERA_GRAPH)
	add_subdirectory(graph)
endif (BUILD_CAMERA_GRAPH)

# ---[ optimiser.
if (BUILD_OPTIMISER)
	add_subdirectory(optimise)
endif (BUILD_OPTIMISER)

# ---[ examples.
if (BUILD_EXAMPLES)
	add_subdirectory(examples)
endif (BUILD_EXAMPLES)
