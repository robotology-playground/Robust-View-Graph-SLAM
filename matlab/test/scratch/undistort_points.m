function [u2,v2] = undistort_points(u,v,CameraParams)

% undistort_points - Function computes the undistorted pixel locations 

% Transform to normal coordinates
xd = (u - CameraParams.u0)./CameraParams.fu;
yd = (v - CameraParams.v0)./CameraParams.fv;

k = CameraParams.K;

% Loop through points
NumPoints = length(xd);
x = zeros(NumPoints,1);
y = zeros(NumPoints,1);
for i = 1:NumPoints
    
    % Iteratively Solve undistorted normal coordinates
    x(i,1) = xd(i,1);
	y(i,1) = yd(i,1);
    for j = 1:20
        r_2 = x(i,1)^2 + y(i,1)^2;
        k_radial = 1 + k(1)*r_2 + k(2)*r_2^2 + k(5)*r_2^3;
        delta_x = 2*k(3)*x(i,1)*y(i,1) + k(4)*(r_2 + 2*x(i,1)^2);
        delta_y = k(3)*(r_2 + 2*y(i,1)^2) + 2*k(4)*x(i,1)*y(i,1);
        x(i,1) = (xd(i,1) - delta_x)/k_radial;
        y(i,1) = (yd(i,1) - delta_y)/k_radial;
    end
    
end

% Transform back into pixel coordinates
u2 = CameraParams.fu.*x + CameraParams.u0;
v2 = CameraParams.fv.*y + CameraParams.v0;

% end of undistort_points

