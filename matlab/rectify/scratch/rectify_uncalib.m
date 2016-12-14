function [newim1, newim2, b, H1, H2, newe] = rectify_uncalib(im1, im2, x1, x2, F12)

% rectify_uncalib:rectifies two images to achieve correspondences of scanlines.
%
%   [newim1, newim2] = rectify_uncalib(im1, im2, x1, x2, F12)
%   [newim1, newim2, b, H1, H2] = rectify_uncalib(im1, im2, x1, x2, F12)
%   [newim1, newim2, b, H1, H2, newe] = rectify_uncalib(im1, im2, x1, x2, F12)
%   rectifies two images, im1 and im2, using the fundamental matrix, F12, which
%   must satisfy the epipolar equation:
%                          (x1)
%      (x2, y2, 1) * F12 * (y1) = 0
%                          ( 1)
% 
%   The input arguments are:
%   - im1 and im2 should both be an m-by-n array of doubles (or uint8) for some
%     values m and n
%   - x1 and x2 should both be 3-by-n matrix, where each column of the matrix
%     is an image point in homogeneous coordinates.  Corresponding columns in
%     the two matrices contain corresponding points in the two images.
%   - F12 must be a 3-by-3 rank-2 matrix.  Fundamental matrices can be computed
%     using one of Torr's routines (available for download on his Microsoft home
%     page) or Zhang's home page.
%
%     ***Note: the fundamental matrix is assumed to be computed with the following
%              image coordinate systems in both images being adopted: the origins
%              of the image coordinate systems are at the centres of the images;
%              the x-axes point to the right, y-axes up.
%
%   The output arguments are:
%   - the two new rectified images, newim1 and newim2.
%   - (optional) the bounding box, b, of the form [minx,miny,maxx,maxy] which bound
%     the new images newim1 and newim2.
%   - (optional) H1 and H2 are the computed rectification transformations.
%   - (optional) newe is the new epipole in the second image after the rectification.
%     On return, newe is always set [1;0;0] if the horizontal scanlines correspond
%     or [0;1;0] if the vertical scanlines correspond.
%
%   The implementation here is based on that described in the paper:
%   Richard I. Hartley, "Theory and Practice of Projective Rectification"
%   International Journal of Computer Vision, vol 35, no 2, pages 115-127, 1999.
%

% compute the epipoles
%---------------------%
%[e1,e2]=e1_e2_from_F(F12);
e1 = null(F12);
e2 = null(F12');

% Check that both epipoles are outside the image boundary (a condition for
% Hartley's method)
%---------------------%
% not using check_epipoles_centred cause x1,x2 weren't centred when F12 
% was calculated
[m,n,p]=size(im1);
if check_epipoles(e1,e2,[m,n]); 
    error(['Check the epipoles; ',...
        'at least one of them is inside the image,', ... 
        'non-existing, or badly conditioned']);
end

% compute the two 3-by-4 projection matrices
% result 9.15 Hartley, page 256
%---------------------%
[P1,P2]=P1_P2_from_F(F12);

% check if points are homogenous coordinates
[x1,x2]=fix_input_points(x1,x2);

% Computes the rectification transformation matrices
% section 11.12.1 Hartley
%---------------------%
if nargout == 2
   [H1,H2] = rectify_transform(P1, P2, e1, e2, x1, x2);
elseif nargout == 5
   [H1,H2,newe,G2,R2] = rectify_transform(P1, P2, e1, e2, x1, x2);
end
   
% do rectification
%---------------------%
nx = size(im1,2)/2;  ny = size(im1,1)/2;
% look for the smallest image size that encloses all the mapped corners
corners1 = pflat(H1*[-nx -ny 1; nx -ny 1; -nx ny 1; nx ny 1]');

nx = size(im2,2)/2;  ny = size(im2,1)/2;
corners2 = pflat(H2*[-nx -ny 1; nx -ny 1; -nx ny 1; nx ny 1]');

corners = [corners1 corners2];
minx = floor(min(corners(1,:))-1);
miny = floor(min(corners(2,:))-1);
maxx = ceil(max(corners(1,:))+1);
maxy = ceil(max(corners(2,:))+1);

b = [minx,miny,maxx,maxy];
newim1 = rectify_im(im1, H1, b);
newim2 = rectify_im(im2, H2, b);

return

