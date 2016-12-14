function val = bilinear_interpolate(im, rc)

%BILINEAR_INTERPOLATE performs bilinear interpolation.
%
%   val = bilinear_interpolate(im, rc) returns the bilinearly interpolated
%   pixel values at the given list of positions, rc.
%
%   Input arguments:
%   - im should be a m-by-n-by-3 (colour) image or a m-by-n (grey scale) image.
%   - rc should be a 2-by-k matrix where k is the number of positions
%     whose values are to be computed via bilinear interpolation.  The first
%     row of the matrix should contain the row-coordinates and the second row the
%     column-coordinates of these positions.  It is important that all the
%     (row,column)-positions stored in the matrix, rc, are within the boundary
%     of the image im.
%
%   Output argument:
%   - val is the output k-by-3 (if image im has 3 colour bands) or k-by-1
%     (if image im is a grey scale image) matrix containing the interpolated pixel
%     values.
%
%Created July 2002.
%Last modified September 2003.
%
%Copyright Du Huynh
%The University of Western Australia
%School of Computer Science and Software Engineering

nrows = size(im,1);  ncols = size(im,2);
% number of bands (3 for coloured images; 1 for gray scale images)
nobands = size(im,3);

% check that all the entries in matrix rc are within the boundary of image im.
if sum(rc(1,:) < 1 | rc(1,:) > nrows | rc(2,:) < 1 | rc(2,:) > ncols)
   error('bilinear_interpolate: elements of the rc matrix must be within the image boundary');
end


% The four corner points used in the bilinear interpolation:
%  c4      c3
%  +--------+
%  |        |
%  | o      |     (a point o which is stored in a column vector of rc
%  |        |     and the two corner points used for its bilinear interpoation)
%  +--------+
%  c1      c2
% The row-column coordinate system is used.  All the corner points c1, c2, c3,
% and c4 are two vectors, whose 1st components contains the row coordinates
% and 2nd components contains the column coordinates.

% note that we should have given a small margin for variable rc
% so that the four corner pixels surrounding rc are within the
% image boundary of im.  This condition should be enforced in the
% caller of this function.
c4 = (floor(rc));
c2 = (ceil(rc));
c3 = ([c4(1,:); c2(2,:)]);
c1 = ([c2(1,:); c4(2,:)]);

% d(diffRC_idx) = (rc(diffRC_idx) - c1) ./ (c4 - c1 + eps);
%
% the interpolation procedure above fails for those points in
% rc whose x- or y- component is a whole number (in which case,
% the respective components of these points in the c1 and c4 matrices
% would be the same.  the formula for d below would cause a division
% by zero problem.
sameC_idx = find(c2(2,:) == c4(2,:));
sameR_idx = find(c2(1,:) == c4(1,:));
diffRC_idx = find(c2(1,:) ~= c4(1,:) & c2(2,:) ~= c4(2,:));
% now the formula for d can be safely applied...
d(:,diffRC_idx) = (rc(:,diffRC_idx) - c4(:,diffRC_idx)) ./ ...
   (c2(:,diffRC_idx) - c4(:,diffRC_idx) + eps);
d(2,sameC_idx) = 0;
d(1,sameC_idx) = (rc(1,sameC_idx) - c4(1,sameC_idx)) ./ ...
   (c2(1,sameC_idx) - c4(1,sameC_idx) + eps);
d(1,sameR_idx) = 0;
d(2,sameR_idx) = (rc(2,sameR_idx) - c4(2,sameR_idx)) ./ ...
   (c2(2,sameR_idx) - c4(2,sameR_idx) + eps);

% convert c1, c2, c3, c4 into 1D array for fast retrieval of image
% intensity from im
c1 = (c1(2,:)-1)*nrows + c1(1,:);
c2 = (c2(2,:)-1)*nrows + c2(1,:);
c3 = (c3(2,:)-1)*nrows + c3(1,:);
c4 = (c4(2,:)-1)*nrows + c4(1,:);

im = reshape(im, nrows*ncols, nobands);
c1val = im(c1,:);
c2val = im(c2,:);
c3val = im(c3,:);
c4val = im(c4,:);

for i=1:nobands
   val(:,i) = (1-d(1,:)').*(1-d(2,:)').*double(c4val(:,i)) + ...
      d(1,:)'.*(1-d(2,:)').*double(c1val(:,i)) + ...
      (1-d(1,:)').*d(2,:)'.*double(c3val(:,i)) + ...
      d(1,:)'.*d(2,:)'.*double(c2val(:,i));
   val(:,i) = uint8(val(:,i));
end

return
