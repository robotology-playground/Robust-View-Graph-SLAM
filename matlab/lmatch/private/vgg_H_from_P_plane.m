% vgg_H_from_P_plane  Given 3D plane and camera, returns 4-by-3 matrix mapping image points to points on the plane.
%
% H = vgg_H_from_P_plane(A,P), where
%   A ... size (4,1), scene plane
%   P ... size (3,4), camera matrix
%   H ... size (4,3), matrix such that X = H*x, where x (size (3,1)) is an image point
%     and X is a scene point lying on the scene plane A.

% T. Werner, March 2003

function H = vgg_H_from_P_plane(A,P)

H = -[vgg_wedge([A; P(2,:); P(3,:)]),...
      vgg_wedge([P(1,:); A; P(3,:)]),...
      vgg_wedge([P(1,:); P(2,:); A])] /...
    (A*vgg_wedge(P));

return
