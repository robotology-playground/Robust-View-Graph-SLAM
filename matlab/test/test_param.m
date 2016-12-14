

% setup folders, parameters and read images
% Tariq Abuhashim - August 2014, iCub

function [K_left, K_right, d_left, d_right, Rc, Tc, imgsiz] = test_param( )


% intrinsic and extrinsic calibration parameters.

% % left camera is at the origin.
% % the camera co-ordinate system is defined as follows: x-left, y-up, z-forward
% % (K_left) and (K_right) are the left and right intrinsic parameters matrix
% % translation is (tc), and is in mm, while rotation matrix is (Rc). The 
% % camera distortion parameters are d_left and d_right
% K_left=[1246.561673900704,0,532.287936343425;
%     0,1247.092336797580,383.604069688392;
%     0,0,1];
% K_right=[1244.805048771648,0,506.635543659909;
%     0,1245.815687282771,390.704813552426;
%     0,0,1];
% d_left=[0 0 0 0 0];
% d_right=[0 0 0 0 0];
% Rc=[0.999985020769994   0.001356661749366  -0.005302612991045;
%     -0.001329315763428   0.999985817981778   0.005157204176297;
%     0.005309534370930  -0.005150078078314   0.999972642396056];
% tc=[-750.6666877284111;-3.857678457065100;3.0959446910902];

% % iCub01
% K_left = [232.921,0,162.202;
%           0,232.43,125.738;
%           0,0,1];
% K_right = [234.474,0,147.91;
%            0,234.169,121.982;
%            0,0,1];
% d_left = [0 0 0 0 0];
% d_right = [0 0 0 0 0];
% Rc = eye(3);
% Tc = [0;0;0];
% imgsiz = [240 320 3];

% iCub01 (high resolution undistorted)
K_left = [426.11508, 0,      336.14615;
          0,       432.00070, 241.68977;
          0,       0,      1];
K_right = [xxx,0,xxx;
           0,xxx,xxx;
           0,0,1];
      
d_left = [ -0.13342   0.07130   0.02253   -0.00805  0.00000 ];
d_right =  [ xxx   xxx   xxx   xxx  xxx ];
Rc = eye(3);
Tc = [0;0;0];
imgsiz = [240 320 3];

% iCub01 (high resolution distorted)