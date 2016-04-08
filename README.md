# Robust-View-Graph-SLAM
C++ implementation of view-graph SLAM using nonlinear least-squares

Code structure tree to be converted from Matlab to C++:

get_aligned_point_matches  

    |___ get_bundle_images
    
            |___ demosaic (Bayer decoding, possibly using OpenCV)
            
    |___ initialise_keyframe
    
            |___ test_triangulate_inverse_depth
            
            |___ test_triangulate_jacobian_inverse_depth
            
            |___ get_scan_from_range
            
    |___ track_bundle_points
    
            |___ calcOpticalFlowPyrLK (done, already in OpenCV)
            
            |___ remove_points_at_infinity
            
    |___ calibrate_bundle_points
    
            |___ get_intrinsics
            
            |___ remove_lens_distortion
            
            |___ pextend
    
initialise_inverse_depth
    |___
    
optimise_constraint_image_inverse_depth
    |___
