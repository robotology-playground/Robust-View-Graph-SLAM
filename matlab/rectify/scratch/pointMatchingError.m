%%=========================================================================
% Function pointMatchingError returns the matching error when the images
% (I1 and I2) and the points (pts1 and pts2) are transformed by the
% projective transformations (tform1 and tform2). Parameter block specifies
% the window size for computing the error. Parameter range specifies the
% maximum distance which the corresponding points can have.
%%=========================================================================
function matchingError = pointMatchingError(I1, tform1, pts1, ...
  I2, tform2, pts2, block, range)

points1 = transformPointsForward(tform1, pts1.Location);
points2 = transformPointsForward(tform2, pts2.Location);

outView = outputView(I1, tform1, I2, tform2);
J1 = imwarp(I1, tform1, 'OutputView', outView);
J2 = imwarp(I2, tform2, 'OutputView', outView);

htm = vision.TemplateMatcher('Metric', 'Sum of squared differences',...
  'OutputValue', 'Metric matrix');
count = 0;
matchingError = 0;
for idx = 1: size(points1, 1)
  p1 = round(points1(idx, :));
  [T1, flag1] = cropImage(J1, p1, block);
  p2 = round(points2(idx, :));
  [T2, flag2] = cropImage(J2, p2, block+range);
  if flag1 && flag2
    metricMatrix = step(htm, T2, T1);
    matchingError = matchingError + min(min(metricMatrix));
    count = count + 1;
  end
end

if count > 0
  matchingError = matchingError / count;
else
  matchingError = inf;
end