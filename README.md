
- Low-level vision
- 
- [ ] get aligned point matches
    - [ ] get image
    - [ ] get features
        - [ ] get sift
        - [ ] get fast
        - [ ] get kaze
    - [ ] track features
        - [ ] track sift
        - [ ] track fast
        - [ ] track kaze
        - [ ] remove points_at_infinity
    - [ ] calibrate image_points
        - [ ] get intrinsics
        - [ ] remove lens distortion
- [ ] initialise linearisation point
    - [ ] essential matrix with RANSAC
        - [ ] bucketing.
        - [ ] five-Points Algorithm.
        - [ ] compute R and T.
        - [ ] triangulate inverse depth
        - [ ] resolve ambiguities.
    - [ ] kinematics
    - [ ] update visibilities
    - [ ] build the state vector (\mathbf{x}_c \mathbf{x}_f)^\top
- [ ] integrate Ceres-Solver.

- RecoverMoments
- 

- PwgOptimiser
-
- [x] initialise_info_matrix.
- [x] generate_constraints_info_Mviews.
- [x] update_info_matrix_Mviews.
- [x] constraints_addition_inverse_depth.
- [x] constraints_removal_inverse_depth.
- [x] compute_gate_inverse_depth_Mviews.
- [ ] solving and inversion (RecoverMoments).
    
- GraphOptimiser
- 
- [ ] assign_constraint_weight
- [ ] run_pose_graph_estimation
    - [ ] initialise
        - [ ] pose_generate_spanning_tree
        - [ ] pose_generate_sequential
        - [ ] initialise_info_matrix
        - [ ] generate_joint_info_matrix
    - [ ] optimise_constraint_graph
        - [ ] generate_constraint_info_pose
        - [ ] constraint_graph_add
        - [ ] constraint_graph_subtract
        - [ ] compute_gate_graph
    
Compile vlFeat
--------------
First, in the makefile, comment out the line #include make/matlab.mak. Then, "make". 
To compule our code with vlFeat, in the CMakeLists.txt, update the path to where you compiled vlFeat:

set(VLFEAT_INCLUDE_DIR /path/to/vlfeat)
message("-- Using VLFeat: ${VLFEAT_INCLUDE_DIR}")
include_directories(${VLFEAT_INCLUDE_DIR})
find_library(VLFEAT_LIB NAMES vl PATHS /path/to/vlfeat/bin/glnxa64)
if (EXISTS ${VLFEAT_LIB})
	message("-- VLFEAT libs: ${VLFEAT_LIB}")
endif(EXISTS ${VLFEAT_LIB})

Install libgsl: 
--------------
sudo apt-get update
sudo apt-get install libgsl-dev
