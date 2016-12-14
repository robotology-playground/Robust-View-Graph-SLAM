% n2 = epibeam_search(EB,x1)  Finds segments in image 2 satisfying epipolar beam constraint wrt given segment in image 1.
%
% EB ... output of epibeam_init
% x1 ... double (3,2), query segment in image 1 (homogeneous)
% n2 ... double (1,:), segments in image 2 satisfying the constraint

function n2 = epibeam_search(EB,x1)

% n2 := lines in image 2 that have consistent orientation with line x1
p1 = det([x1 EB.e1]) > 0;
if p1
  n2 = find(~EB.p2);
else
  n2 = find(EB.p2);
end
if isempty(n2)
  return
end

% epipolars of line x1 in image 2
m2 = (EB.F*x1)';

% remove segments with end point above first epipolar 
n2( (-1)^p1*m2(1,:)*EB.y2(:,n2) < 0 ) = [];
if isempty(n2)
  return
end

% remove segments with start point under second epipolar
n2( (-1)^p1*m2(2,:)*EB.x2(:,n2) > 0 ) = [];
if isempty(n2)
  return
end

% remove segments that are behind any clipping plane
l2 = reshape( (-1)^p1*cros(x1)'*EB.H, [3 size(EB.H,2)/3] )';
n2( any(l2*EB.x2(:,n2) < 0, 1) ) = [];
n2( any(l2*EB.y2(:,n2) < 0, 1) ) = [];

return
