
# Copyright (C) 2016 iCub Facility - Istituto Italiano di Tecnologia
# Authors: Tariq Abuhashim, Nicolo' Genesio
# E-mail: t.abuhashim@gmail.com, nicogene@hotmail.it
# Date: Nov 2016
# Acknowledgement: This research has received funding from the European Union’s 
# Seventh Framework Programme for research, technological development and demonstration 
# under grant agreement No. 611909(KoroiBot).
# License: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT

cmake_minimum_required(VERSION 3.2.2)

#SET(PROJECTNAME vgSLAM)
#PROJECT(${PROJECTNAME})

########################################################################
################# ONLY change inside this region #######################
########################################################################
SET(suitesparse /home/tariq/Dev/suitesparse)
SET(vlfeat /home/tariq/Dev/vlfeat-0.9.20)
SET(matlab /home/tariq/Dev/matlab)
########################################################################
########################################################################
########################################################################

OPTION( BUILD_CAMERA_GRAPH	"Build camera graph"	OFF )
OPTION( BUILD_OPTIMISER		"Enable optimisation"	OFF )

if (BUILD_OPTIMISER)
	add_subdirectory(optimise)
endif (BUILD_OPTIMISER)

if (BUILD_CAMERA_GRAPH)
	add_subdirectory(graph)
endif (BUILD_CAMERA_GRAPH)
