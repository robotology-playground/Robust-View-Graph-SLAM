# Robust-View-Graph-SLAM
C++ implementation of view-graph SLAM using nonlinear least-squares

Code structure tree to be converted from Matlab to C++:

get_aligned_point_matches  
    get_bundle_images
        demosaic (Bayer decoding, possibly using OpenCV)
    initialise_keyframe
        test_triangulate_inverse_depth
        test_triangulate_jacobian_inverse_depth
        get_scan_from_range
    track_bundle_points
        calcOpticalFlowPyrLK (done, already in OpenCV)
        remove_points_at_infinity
    calibrate_bundle_points
        get_intrinsics
        remove_lens_distortion
        pextend
initialise_inverse_depth
optimise_constraint_image_inverse_depth
