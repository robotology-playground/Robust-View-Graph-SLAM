function outView = outputView(I1, tform1, I2, tform2)

numRows = size(I1, 1);
numCols = size(I1, 2);
inPts = [1, 1; 1, numRows; numCols, numRows; numCols, 1];
outPts(1:4,1:2) = transformPointsForward(tform1, inPts);
numRows = size(I2, 1);
numCols = size(I2, 2);
inPts = [1, 1; 1, numRows; numCols, numRows; numCols, 1];
outPts(5:8,1:2) = transformPointsForward(tform2, inPts);

xSort   = sort(outPts(:,1));
ySort   = sort(outPts(:,2));
xLim(1) = ceil(xSort(4)) - 0.5;
xLim(2) = floor(xSort(5)) + 0.5;
yLim(1) = ceil(ySort(4)) - 0.5;
yLim(2) = floor(ySort(5)) + 0.5;
width   = xLim(2) - xLim(1) - 1;
height  = yLim(2) - yLim(1) - 1;
outView = imref2d([height, width], xLim, yLim);