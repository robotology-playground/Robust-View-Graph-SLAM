
function [Irectified, Irect1, Irect2] = demo_rectification(file1, file2)

%file1 = 'parkinglot_left.png';
%file2 = 'parkinglot_right.png';

%--------------------------------------------------------------------------
% Read Stereo Image Pair
I1 = imread(file1);
if size(I1, 3) == 3
  I1 = rgb2gray(I1);
end

I2 = imread(file2);
if size(I2, 3) == 3
  I2 = rgb2gray(I2);
end

%--------------------------------------------------------------------------
% Try to rectify the images with up to 5 iterations. If an iteration fails
% to rectify the images, the function will change the parameters to use
% more memory and computational time, and to do rectification again.
numExperiment = 0;
while numExperiment < 5
  %------------------------------------------------------------------------
  % Step 1. Set Parameters
  % Set the parameters used in the algorithm. They are adjusted based on
  % the number of experiments. The effect of a larger Q is indicated after
  % the parameter.
  numExperiment = numExperiment + 1;
  Q = numExperiment;
  metricThreshold = .1 / Q;     % Leads to more interesting points
  matchThreshold = 50 / Q;         % Leads to more matching point pairs
  numTrials = 10000 * Q;          % Allow lower inlier ratio
  confidence = 100 - 0.01 / Q;    % Allow lower inlier ratio
  rangeThreshold = [5, 2] * Q;    % Allow larger disparity
  
  %------------------------------------------------------------------------
  % Step 2. Collect Interest Points from Each Image
  blobs1 = detectSURFFeatures(I1, 'MetricThreshold', metricThreshold);
  blobs2 = detectSURFFeatures(I2, 'MetricThreshold', metricThreshold);
  
  %------------------------------------------------------------------------
  % Step 3. Select Correspondences Between Points Based on SURF Features
  [features1, validBlobs1] = extractFeatures(I1, blobs1);
  [features2, validBlobs2] = extractFeatures(I2, blobs2);
  indexPairs = matchFeatures(features1, features2, 'Metric', 'SAD', ...
    'MatchThreshold', matchThreshold);
  
  % Retrieve locations of matched points for each image
  matchedPoints1 = validBlobs1(indexPairs(:,1),:);
  matchedPoints2 = validBlobs2(indexPairs(:,2),:);
  
  %------------------------------------------------------------------------
  % Step 4. Remove Outliers Using Epipolar Constraints
  [fMatrix, epipolarInliers, status] = estimateFundamentalMatrix(...
    matchedPoints1, matchedPoints2, 'Method', 'RANSAC', ...
    'NumTrials', numTrials, 'DistanceThreshold', 0.1, ...
    'Confidence', confidence);
  
  % If the function fails to find enough inliers or if either epipole is
  % inside the image, the images cannot be rectified. The function will
  % adjust the parameters and start another iteration.
  if status ~= 0 || isEpipoleInImage(fMatrix, size(I1)) ...
      || isEpipoleInImage(fMatrix', size(I2))
    continue;
  end
  
  inlierPoints1 = matchedPoints1(epipolarInliers, :);
  inlierPoints2 = matchedPoints2(epipolarInliers, :);
  
  %------------------------------------------------------------------------
  % Step 5. Compute Rectification Transformations
  [t1, t2] = estimateUncalibratedRectification(fMatrix, ...
    inlierPoints1.Location, inlierPoints2.Location, size(I2));
  tform1 = projective2d(t1);
  tform2 = projective2d(t2);

  %------------------------------------------------------------------------
  % Step 6. Check Quality Of Rectification And Generate Rectification
  % Composite
  matchingError = pointMatchingError(I1, tform1, inlierPoints1, ...
    I2, tform2, inlierPoints2, [7, 7], rangeThreshold);
  
  % maximumMatchingError is the maximum value of registration error.
  % Increase this parameter if the images are very different.
  maximumMatchingError = 0.5;
  if matchingError < maximumMatchingError
    [Irectified, Irect1, Irect2] = TransformImagePair(I1, tform1, I2, tform2);
    %figure, imshow(Irectified);
    %title('Rectified Stereo Images (Red - Left Image, Cyan - Right Image)');
    %figure, imshow(Irect1);
    %figure, imshow(Irect2);
    imwrite(Irect1,'left.png');
    imwrite(Irect2,'right.png');
    return;
  end
end