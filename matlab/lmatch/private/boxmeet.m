function c = boxmeet(a,b)

% c = boxmeet(a,b)  Intersection of boxes. If a, b are disjoint, c = [nan nan; ...]. 
%
% The corners must be sorted, use boxsortu.

c(:,1) = max(a(:,1),b(:,1));
c(:,2) = min(a(:,2),b(:,2));

if any(c(:,1) > c(:,2)) | any(isnan([a(:); b(:)]))
  c = ones(size(a,1),1)*[nan nan];
end

return
