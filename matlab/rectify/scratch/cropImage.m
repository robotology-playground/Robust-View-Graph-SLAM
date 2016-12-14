%%=========================================================================
% Function cropImage returns T, the sub-image of I locating at loc with
% size of 2*block+1. The function also returns flag, which is true when T
% is completely inside I and false otherwise.
%%=========================================================================

function [T, flag] = cropImage(I, loc, block)
numRows = size(I, 1);
numCols = size(I, 2);
if loc(1) > block(1) && loc(1) <= numCols-block(1) ...
    && loc(2) > block(2) && loc(2) <= numRows-block(2)
  T = I(loc(2)-block(2): loc(2)+block(2), loc(1)-block(1): loc(1)+block(1));
  flag = true;
else
  T = zeros(2*block+1);
  flag = false;
end