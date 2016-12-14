function [u2,v2] = undistort_points(u,v,CameraParams)

% undistort_points - Function undistorts a set of input pixel coordinates

% Loop through points
npoints = length(u);
u2 = zeros(npoints,1);
v2 = zeros(npoints,1);
for i = 1:npoints
    
    % Get normalised coordinates
    xd = (u(i,1) - CameraParams.u0)/CameraParams.fu;
	yd = (v(i,1) - CameraParams.v0)/CameraParams.fv;
    
    % iteratively solve for undistorted normal coordinate
    x = xd;
	y = yd;
    for j = 1:20
        r_2 = x^2 + y^2;
		k_radial = 1 + CameraParams.K(1)*r_2 + CameraParams.K(2)*(r_2^2) + CameraParams.K(5)*(r_2^3);
		delta_x = 2*CameraParams.K(3)*x*y + CameraParams.K(4)*(r_2 + 2*x^2);
		delta_y = CameraParams.K(3)*(r_2 + 2*y^2) + 2*CameraParams.K(4)*x*y;
		x = (xd - delta_x)/k_radial;
		y = (yd - delta_y)/k_radial;
    end
    
    % Transform to Pixels
	u2(i,1) = CameraParams.fu*x + CameraParams.u0;
	v2(i,1) = CameraParams.fv*y + CameraParams.v0;
        
end

% end of undistort_points

