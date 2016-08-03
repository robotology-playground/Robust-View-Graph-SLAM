project( "PwgOptimiser" )

cmake_minimum_required(VERSION 2.6)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# Find eigen.
find_path(Eigen_INCLUDE_DIR NAMES signature_of_eigen3_matrix_library PATHS $ENV{EIGEN_ROOT} /usr/include/eigen3 /usr/local/include/eigen3 )
if(NOT EXISTS ${Eigen_INCLUDE_DIR})
    message( "could not find Eigen - should be the path containing the Eigen 3 directory and /Eigen/Eigen, /Eigen/Core, etc..." )
else(NOT EXISTS ${Eigen_INCLUDE_DIR})
    include_directories(${Eigen_INCLUDE_DIR})
endif(NOT EXISTS ${Eigen_INCLUDE_DIR})

# Find snark.
#find_package(snark PATHS ${snark_DIR} /usr/local/CMakeFiles)


file(GLOB pwg_source *.cpp *.c)

# Filter out other sources
file(GLOB to_remove main_motion_comp* test* *old* main_piped*)
list(REMOVE_ITEM icp_source ${to_remove})

add_executable(PwgOptimiser ${pwg_source})
target_link_libraries(PwgOptimiser)

