function [Iout, imageTransformed1, imageTransformed2] = TransformImagePair(I1, tform1, I2, tform2)
%cvexTransformImagePair Return the common area of transformed images
%  IOUT = cvexTransformImagePair(I1, TFORM1, I2, TFORM2)) transforms images
%  I1 and I2 using projective transformations TFORM1 and TFORM2, and
%  returns a composite image IOUT made by their common area. The
%  transformation of I1 is shown in red while that of I2 is shown in cyan.
%
%  Notes
%  -----
%  This is a helper function in support of examples and is subject to
%  change or removal in a future release.
%
%  See also estimateUncalibratedRectification, videorectification,
%     cvexTransformImagePair

%  Copyright  The MathWorks, Inc.
%  $Revision: 1.1.6.2 $Date: 2012/11/15 14:56:43 $

%--------------------------------------------------------------------------
% Compute the transformed location of image corners.
numRows = size(I1, 1);
numCols = size(I1, 2);
inPts = [1, 1; 1, numRows; numCols, numRows; numCols, 1];
outPts(1:4,1:2) = transformPointsForward(tform1, inPts);
numRows = size(I2, 1);
numCols = size(I2, 2);
inPts = [1, 1; 1, numRows; numCols, numRows; numCols, 1];
outPts(5:8,1:2) = transformPointsForward(tform2, inPts);

%--------------------------------------------------------------------------
% Compute the common rectangular area of the transformed images.
xSort   = sort(outPts(:,1));
ySort   = sort(outPts(:,2));
xLim(1) = ceil(xSort(4)) - 0.5;
xLim(2) = floor(xSort(5)) + 0.5;
yLim(1) = ceil(ySort(4)) - 0.5;
yLim(2) = floor(ySort(5)) + 0.5;
width   = xLim(2) - xLim(1) - 1;
height  = yLim(2) - yLim(1) - 1;
outputView = imref2d([height, width], xLim, yLim);

%--------------------------------------------------------------------------
% Generate a composite made by the common rectangular area of the
% transformed images.
imageTransformed1 = imwarp(I1, tform1, 'OutputView', outputView);
imageTransformed2 = imwarp(I2, tform2, 'OutputView', outputView);
Iout(:,:,1) = imageTransformed1;
Iout(:,:,2) = imageTransformed2;
Iout(:,:,3) = imageTransformed2;
