# Robust-View-Graph-SLAM
C++ implementation of view-graph SLAM using nonlinear least-squares

Code structure tree to be converted from Matlab to C++:

- batch_motion_and_map_inverse_depth
    - get_aligned_point_matches
        - get_bundle_images
            - demosaic (Bayer decoding, possibly using OpenCV)
        - initialise_keyframe
            - test_triangulate_inverse_depth
            - test_triangulate_jacobian_inverse_depth
            - get_scan_from_range
        - track_bundle_points
            - calcOpticalFlowPyrLK (done, already in OpenCV)
            - remove_points_at_infinity
        - calibrate_bundle_points
            - get_intrinsics
            - remove_lens_distortion
            - pextend
    - initialise_inverse_depth
    - optimise_constraint_image_inverse_depth
        - generate_constraints_info_Mviews          .... DONE !!!
        - initialise_info_matrix
        - generate_constraints_info_Mviews          .... DONE !!!
        - update_info_matrix_Mviews
        - constraints_addition_inverse_depth
        - constraints_removal_inverse_depth
        - compute_gate_inverse_depth_Mviews         .... DONE !!!
    
- assign_constraint_weight

- run_pose_graph_estimation
    - initialise
        - pose_generate_spanning_tree
        - pose_generate_sequential
        - initialise_info_matrix
        - generate_joint_info_matrix
    - optimise_constraint_graph
        - generate_constraint_info_pose
        - constraint_graph_add
        - constraint_graph_subtract
        - compute_gate_graph                        .... DONE !!!
        
- recover_moments
