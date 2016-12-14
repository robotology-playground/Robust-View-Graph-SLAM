function [H1,H2]=epipolar_rectification(F,pts1,pts2,cx,cy)


% Epipolar geometry
[e1,e2]=epipoles_from_F(F);
% we should have
% (FM' * e1 ~ 0) and (FM  * e2 ~ 0)

% find the homography Hprime that takes e1 to infinity then normalize things 
Hprime = map_to_infinity( e1, cx, cy );
e1_ = Hprime * e1;

% Normalize Hprime so that Hprime*eprime = (1,0,0)'
Hprime = Hprime / e1_(1);
e1_ = Hprime * e1;

% find the matching transform for Hprime that aligns image 2's rows with 
% the rectified image 1's rows and also causes minimum distortion between 
% matching points. 
% Get canonical camera matrices for F12 and compute H0, one possible
% rectification homography for image 2
[P,Pprime] = get_canonical_cameras(F);
M = Pprime(:,1:3);
H0 = Hprime * M;

% Now we need to find H so that the epipolar lines match
% each other, i.e., inv(H)' * l = inv(Hprime)' * lprime
% and the disparity is minimized, i.e.,
% min \sum_i d(H x_i, Hprime xprime_i)^2

% Now that we have homography H0 which will align our epipolar lines, 
% we also want to shear and translate along the X axis so that the rectified 
% corresponding points are as close together as possible in Euclidean distance. 
% The idea is to estimate a matrix HA which performs the translate and shear 
% using least squares.
% First, we tranform the points in image 1 according to Hprime and the 
% points in image 2 according to H0: 
% Transform data initially according to Hprime (img 1) and H0 (img 2)
hpts1 = makehomogeneous(pts1');
hpts1 = Hprime * hpts1;
npts1 = hnormalise(hpts1);
hpts2 = makehomogeneous(pts2');
hpts2 = H0 * hpts2;
npts2 = hnormalise(hpts2);
rmse_x = sqrt( mean( (npts1(1,:) - npts2(1,:) ).^2 ));
rmse_y = sqrt( mean( (npts1(2,:) - npts2(2,:) ).^2 ));
fprintf( 1, 'Before Ha, RMSE for corresponding points in Y: %g X: %g\n', ...
    rmse_y, rmse_x );

% Estimate [ a b c ; 0 1 0 ; 0 0 1 ] aligning H, Hprime
n = size(pts1',2);
A = [ npts2(1,:)', npts2(2,:)', ones(n,1) ];
b = npts1(1,:)';
abc = A\b;
HA = [ abc' ; 0 1 0 ; 0 0 1 ];
H = HA*H0;
hpts2 = makehomogeneous(pts2');
hpts2 = H * hpts2;
npts2 = hnormalise(hpts2);
rmse_x = sqrt( mean( (npts1(1,:) - npts2(1,:) ).^2 ));
rmse_y = sqrt( mean( (npts1(2,:) - npts2(2,:) ).^2 ));
fprintf( 1, 'After Ha, RMSE for corresponding points in Y: %g X: %g\n', ...
    rmse_y, rmse_x );

% Return the homographies as appropriate
H1 = Hprime;
H2 = H;