function newim = rectify_im(im, H, box)

%RECTIFY_IM applies the rectification transformation to a given image.
%
%   newim = rectify_im(im, H, box) applies the rectification transformation
%   H (a 3-by-3 matrix) to the given image, im.  The bounding box
%   which determines the size of the output image newim, should be of the
%   format: [minx,miny,maxx,maxy], where minx and miny can be negative
%   numbers as the origin of the image coordinate system is set to the
%   centre of the image when the transformation H is computed.
%
%   See also rectify_images
%
%   Created in July 2002.
%   Last modified September 2003.
%   Bug fixed November 2004.  (thanks to Tzu Yen Wong)
%
%   Copyright Du Huynh
%   The University of Western Australia
%   School of Computer Science and Software Engineering

minx = box(1); miny = box(2); maxx = box(3); maxy = box(4);
% the two matrices xx and yy returned by meshgrid are of the same dimension
[xx,yy] = meshgrid(minx:maxx, maxy:-1:miny);

% dimensions of the new (rectified) image
new_nrows = size(xx,1);  new_ncols = size(xx,2);

x = reshape(xx, 1, new_nrows*new_ncols); clear('xx');
y = reshape(yy, 1, new_nrows*new_ncols); clear('yy');

invH = inv(H);
len = length(x);
% We will encounter the "out of memory" problem if the image is large.
% So, to get around the problem, we do the following operation in
% several steps
mm = 50000;
idx=1:mm;
while (1)
   if idx(1) > len
      break;
   elseif idx(end) > len
      idx = idx(1:len-idx(1)+1);
   end
   newxy(:,idx) = invH*([x(idx); y(idx); ones(1,length(idx))]);
   newxy(:,idx) = pflat(newxy(:,idx));
   idx = idx + mm;
end
clear('idx');

% convert the x-y image coordinate system (origin at the image centre) to
% row-column coordinate system (origin at the top-left corner) before
% calling bilinear interpolation
nrows = size(im,1);  ncols = size(im,2);
newrc = [0 -1 nrows/2; 1 0 ncols/2; 0 0 1]*newxy;
clear('newxy');

% can't interpolate those points that fall outside the image im.  So discard them.
idx = find(newrc(1,:) >= 1 & newrc(1,:) <= nrows & ...
   newrc(2,:) >= 1 & newrc(2,:) <= ncols);
newrc = newrc(1:2,idx);
x = x(idx);
y = y(idx);

val = [];
len = size(newrc,2);
idx = 1:mm;
while (1)
   if idx(1) > len
      break;
   elseif idx(end) > len
      idx = idx(1:len-idx(1)+1);
   end
   val(idx,:) = bilinear_interpolate(im, newrc(:,idx));
   idx = idx + mm;
end

% compose the new image, newim, and the 2 arrays, x and y,
% into 1D array (note that x and y are still defined relative
% to the image origin at the centre of the image buffer.
% newim = zeros(new_ncols*new_nrows, size(im,3));
% covert into row-column coordinate system also
newim = zeros(new_nrows*new_ncols, 1, size(im,3));
rc = (x-minx)*new_nrows + (maxy-y) + 1;
if ~isempty(val)
   newim(rc,:) = val;
end

newim = reshape(uint8(newim), new_nrows, new_ncols, size(im,3));

return

