function [F,varargout] = fundmatrix_nonlin(xy1xy2, idx, scale)

%FUNDMATRIX_NONLIN computes the fundamental matrix using the 7 point (non-linear) method.
%
%   F = fundmatrix_nonlin(xy1xy2, idx)
%   F = fundmatrix_nonlin(xy1xy2, idx, scale)
%   [F,errs] = fundmatrix_nonlin(xy1xy2, idx)
%   [F,errs] = fundmatrix_nonlin(xy1xy2, idx, scale)
%   [F,errs,avgerr] = fundmatrix_nonlin(xy1xy2, idx)
%   [F,errs,avgerr] = fundmatrix_nonlin(xy1xy2, idx, scale)
%
%   estimates the fundamental matrix using the iterative reweighted least squares
%   method.  The corresponding point coordinates should be stored in the large
%   6-by-n matrix of the following format:
%      xy1xy2 = [x1_1  x1_2 ... x1_n;
%                y1_1  y1_2 ... y1_n;
%                z1_1  z1_2 ... z1_n;
%                x2_1  x2_2 ... x2_n;
%                y2_1  y2_2 ... y2_n;
%                z2_1  z2_2 ... z2_n ]
%   where [x1_i, y1_i, z1_i] in the first three rows is the i-th feature point
%                            (in homogeneous coordinates) in the first image;
%         [x2_i, y2_i, z2_i] in the last three rows is the i-th feature point
%                            (in homogeneous coordinates) in the second image;
%         each column of the matrix stores a pair of corresponding points in
%         homogeneous coordinates.
%
%   The fundmental matrix, F, is a rank-2 3-by-3 matrix satisfying the epipolar
%   constraint:
%                            (x1_i)
%     (x2_i y2_i z2_i) * F * (y1_i) = 0
%                            (z1_i)
%
%   Input arguments:
%   - xy1xy2 should be a 6-by-n matrix where n >= 7.  It is recommended that
%     the coordinates of points in this matrix be relative to an image coordinate
%     system whose origin is at the centre (or principal point, if known) of
%     the image.
%   - idx must be an array of integers in the range [1,n] where n=size(xy1xy2,2).
%     Alternatively, set idx=[] if all the points in xy1xy2 are intended to be
%     used for the estimation of F.
%   - scale (optional) is the scale factor to be used in least squares to keep
%     things well behaved.  If not specified, a scale factor would be automatically
%     computed (Referencee: see Hartley's ICCV'95 paper).  If specified, a sensible
%     scale factor is half of the average of the image dimensions, eg. if the images
%     are of dimensions 400-by-600 then scale should be set to 500/2 = 250.
%
%   Output arguments:
%   - Unlike other estimation method for the fundamental matrix, this function
%     may produce 1 or 3 solutions.  F is thus a cell array of up to length 3,
%     with each cell element contains an estimated fundamental matrix.
%   - Similar to the first output argument, errs is a cell array where each cell
%     element contains the list of squared reprojection errors the feature points
%     from the epipolar lines) of all the feature point in xy1xy2.
%   - avgerr is also a cell array where each cell element contains the average
%     value of the corresponding list of errors stored in the second output argument
%     above.
%
%Created 2001.
%Revised September 2003.
%
%Copyright Du Huynh
%The University of Western Australia
%School of Computer Science and Software Engineering

if isempty(idx)
   idx = 1:size(xy1xy2,2);
elseif any(idx <= 0 | idx > size(xy1xy2,2))
   error('fundmatrix_nonlin: idx must be an array of integers in the range [1,n], where n=size(xy1xy2,2).');   
end

xy = xy1xy2(:,idx);

% check input arguments
if size(xy,1) ~= 6
   error('fundmatrix_nonlin: xy1xy2 must have 6 rows.');
elseif length(idx) < 7
   error('fundmatrix_nonlin: at least 7 pairs of corresponding points should be used.');
end

xy1 = pflat(xy(1:3,:));
xy2 = pflat(xy(4:6,:));

if nargin < 3
   % compute scale factor and shifting
   % this scaling is not quite the same as that described in the Hartley's
   % iccv95 paper (and others' papers), but the difference to the estimated
   % fundamental matrix would be insignificant for other ways of scaling.
   % It is more important that the data points are well scattered in the
   % images and are not degenerate (e.g. the image points are projections
   % of coplanar scene points)
   scale = sqrt(2) / ( max(max([xy1;xy2]))-min(min([xy1;xy2])) );
   T1 = [scale 0 -scale*mean(xy1(1,:));
      0 scale -scale*mean(xy1(2,:));
      0 0 1];
   T2 = [scale 0 -scale*mean(xy2(1,:));
      0 scale -scale*mean(xy2(2,:));
      0 0 1];
elseif scale == 0
   T1 = eye(3);  T2 = eye(3);
else
   T1 = diag([1/scale, 1/scale, 1]);
   T2 = T1;
end

% transform the corresponding points by T1 and T2
xy1new = T1*xy1;
xy2new = T2*xy2;

% A is a 7x9 matrix
A(:,1) = xy2new(1,:)' .* xy1new(1,:)';
A(:,2) = xy2new(1,:)' .* xy1new(2,:)';
A(:,3) = xy2new(1,:)' .* xy1new(3,:)';

A(:,4) = xy2new(2,:)' .* xy1new(1,:)';
A(:,5) = xy2new(2,:)' .* xy1new(2,:)';
A(:,6) = xy2new(2,:)' .* xy1new(3,:)';

A(:,7) = xy2new(3,:)' .* xy1new(1,:)';
A(:,8) = xy2new(3,:)' .* xy1new(2,:)';
A(:,9) = xy2new(3,:)' .* xy1new(3,:)';

% since A'*A is a 9x9 rank-7 matrix, the dimension of its null space is 2.
[U,S,V] = svd(A'*A, 0);
G = reshape(V(:,8), 3, 3)';
H = reshape(V(:,9), 3, 3)';

% undo the scaling
G = T2'*G*T1;  G = G/norm(G,2);
H = T2'*H*T1;  H = H/norm(H,2);

G11=G(1,1); G12=G(1,2); G13=G(1,3);
G21=G(2,1); G22=G(2,2); G23=G(2,3);
G31=G(3,1); G32=G(3,2); G33=G(3,3);

H11=H(1,1); H12=H(1,2); H13=H(1,3);
H21=H(2,1); H22=H(2,2); H23=H(2,3);
H31=H(3,1); H32=H(3,2); H33=H(3,3);

% the fundamental matrix, F, that we are after is a linear combination
% of G and H, i.e. F = alpha*G + (1-alpha)*H.  We impose the det(F)=0
% constraint to find alpha.  The condition det(F)=0 gives a polynomial of degree
% 3 for alpha.  So we would either get one real solution for alpha or three
% real solutions for alpha.  We need to deal with the latter case later.
% Below are the coefficients of the polynomial, obtained from Matlab's
% symbolic maths toolbox.

coeffs = zeros(4,1);
% coeffcient of the alpha^3 term
coeffs(1) = -H31*H12*H23-H11*H22*H33-G11*G22*H33+ ...
   H11*H23*H32-G11*H23*H32-H21*H13*H32+H21*H12*H33+ ...
   H31*H13*H22+G11*G23*H32+H11*G22*H33+G11*G22*G33- ...
   G11*H22*G33-G11*G23*G32-H11*G22*G33+G11*H23*G32+ ...
   G11*H22*H33+H11*H22*G33+H11*G23*G32-H11*G23*H32- ...
   H11*H23*G32-G21*G12*G33+G21*G12*H33+G21*H12*G33- ...
   G21*H12*H33+G21*G13*G32-G21*G13*H32-G21*H13*G32+ ...
   G21*H13*H32+H21*G12*G33-H21*G12*H33-H21*H12*G33- ...
   H21*G13*G32+H21*G13*H32+H21*H13*G32+G31*G12*G23- ...
   G31*G12*H23-G31*H12*G23+G31*H12*H23-G31*G13*G22+ ...
   G31*G13*H22+G31*H13*G22-G31*H13*H22-H31*G12*G23+ ...
   H31*G12*H23+H31*H12*G23+H31*G13*G22-H31*G13*H22- ...
   H31*H13*G22;
% coefficient of the alpha^2 term
coeffs(2) = 3*H31*H12*H23+3*H11*H22*H33+ ...
   G11*G22*H33-3*H11*H23*H32+2*G11*H23*H32+ ...
   3*H21*H13*H32-3*H21*H12*H33-3*H31*H13*H22- ...
   G11*G23*H32-2*H11*G22*H33+G11*H22*G33+ ...
   H11*G22*G33-G11*H23*G32-2*G11*H22*H33- ...
   2*H11*H22*G33-H11*G23*G32+2*H11*G23*H32+ ...
   2*H11*H23*G32-G21*G12*H33-G21*H12*G33+ ...
   2*G21*H12*H33+G21*G13*H32+G21*H13*G32- ...
   2*G21*H13*H32-H21*G12*G33+2*H21*G12*H33+ ...
   2*H21*H12*G33+H21*G13*G32-2*H21*G13*H32- ...
   2*H21*H13*G32+G31*G12*H23+G31*H12*G23- ...
   2*G31*H12*H23-G31*G13*H22-G31*H13*G22+ ...
   2*G31*H13*H22+H31*G12*G23-2*H31*G12*H23- ...
   2*H31*H12*G23-H31*G13*G22+2*H31*G13*H22+ ...
   2*H31*H13*G22;
% coefficient of the alpha term
coeffs(3) = G31*H12*H23+H11*H22*G33- ...
   G21*H12*H33-H31*G13*H22+G11*H22*H33+ ...
   3*H31*H13*H22-H11*G23*H32-3*H21*H13*H32- ...
   G31*H13*H22+3*H11*H23*H32+G21*H13*H32- ...
   3*H31*H12*H23-G11*H23*H32+3*H21*H12*H33- ...
   H31*H13*G22-3*H11*H22*H33+H21*H13*G32+ ...
   H21*G13*H32-H21*H12*G33-H21*G12*H33+ ...
   H31*H12*G23+H31*G12*H23-H11*H23*G32+ ...
   H11*G22*H33;
% coefficient of the constant term
coeffs(4) = H31*H12*H23+H11*H22*H33- ...
   H11*H23*H32+H21*H13*H32-H21*H12*H33-H31*H13*H22;

alphas = roots(coeffs);
no_solns = 0;  % number of real solutions
for i=1:3
   if isreal(alphas(i))
      no_solns = no_solns + 1;
      F{no_solns} = alphas(i)*G + (1-alphas(i))*H;
      F{no_solns} = F{no_solns} / norm(F{no_solns},2);
   end
end

if nargout >= 2
   xy1 = xy1xy2(1:3,:);
   xy2 = xy1xy2(4:6,:);
   for soln=1:no_solns
      % points in the first image cast epipolar lines in the second image under
      % the mapping of the fundamental matrix
      line2 = F{soln}*xy1;
      % similarly...
      line1 = F{soln}'*xy2;
      
      % squared reprojection errors (see Hartley & Zisserman's book, Eq (10.9), p.271)
      % Here, I keep the array of errors rather than summing them together as in Eq (10.9)
      errs{soln} = (sum(line2.*xy2,1)).^2 ./ ...
         (line2(1,:).^2 + line2(2,:).^2 + line1(1,:).^2 + line1(2,:).^2);
   end
   varargout{1} = errs;
end

if nargout >= 3
   % this is the average value of the reprojection errors, i.e. Eq (10.9) divided
   % by the number of corresponding points (see Hartley & Zisserman's book)
   for soln=1:no_solns
      avgerr{soln} = mean(errs{soln});
   end
   varargout{2} = avgerr;
end

return
