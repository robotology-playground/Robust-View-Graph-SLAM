#.rst:
# FindPARMETIS
# -------
#
# Find the PARMETIS.
#
# Once done this will define the following variables::
#
#   PARMETIS_LIBRARIES       - PARMETIS libraries 
#   PARMETIS_FOUND           - if false, you cannot build anything that requires PARMETIS

include(FindPackageHandleStandardArgs)
include(SelectLibraryConfigurations)

find_library(PARMETIS_LIBRARY_RELEASE
             NAMES parmetis
             DOC "PARMETIS library file (release version)")
find_library(PARMETIS_LIBRARY_DEBUG
             NAMES parmetisd 
             DOC "PARMETIS library file (debug version)") 

mark_as_advanced(PARMETIS_LIBRARY_RELEASE
                 PARMETIS_LIBRARY_DEBUG)

select_library_configurations(PARMETIS)

set(PARMETIS_LIBRARIES ${PARMETIS_LIBRARY}) 

find_package_handle_standard_args(PARMETIS
                                  FOUND_VAR PARMETIS_FOUND
                                  REQUIRED_VARS PARMETIS_LIBRARIES)

# Set package properties if FeatureSummary was included
if(COMMAND set_package_properties)
    set_package_properties(Libedit PROPERTIES DESCRIPTION "ParMETIS is an MPI-based parallel library that implements a variety of algorithms for partitioning unstructured graphs, meshes, and for computing fill-reducing orderings of sparse matrices. "
                                              URL "http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview")
endif()
