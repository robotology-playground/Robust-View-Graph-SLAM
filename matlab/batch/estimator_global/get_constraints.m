function [C, idx] = get_constraints(C, i, j)
edges = vertcat(C.edge);
if nargin < 3;
    if length(i) == 2;
        idx = find( edges(:,2) == i(1) | edges(:,2) == i(2) );%
    else
        idx = find( edges(:,2) < i+1 );%
    end
else
    idx = find( edges(:,1) > i-1 & edges(:,2) < j+1 );
end
% copy constraints
C = C(idx);
% re-number edges from 1 to j-i
for ii = 1:length(C)
    C(ii).edge = C(ii).edge - i + 1;
end
return