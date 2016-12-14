% EB = epibeam_init(x2,y2,P1,P2,A)  Initializes epibeam_search
%
% x2, y2 ... double(2,N2), end points of line sgments in image 2
% A ... double(?,4), clipping planes, typically plane at infty (optional)
% P1, P2 ... double(3,4), cameras, must be oriented consistently with each other and with A.
%  Changing orientation of P1, P2 can be done by :-
%   - either swapping signs
%   - or interchanging their rows (= mirroring the whole scene)
%
function EB = epibeam_init(x2,y2,P1,P2,A)

EB.F = vgg_F_from_P(P1,P2);
EB.F(:) = normx(EB.F(:));
EB.e1 = normx(P1*vgg_wedge(P2));
EB.e2 = normx(P2*vgg_wedge(P1));
EB.x2 = x2;
EB.y2 = y2;

if nargin < 5
  A = [];
end
EB.H = zeros(3,0);
for n = 1:size(A,1)
  EB.H = [EB.H P1*vgg_H_from_P_plane(A(n,:),P2)];
end

for n = 1:size(x2,2)
  EB.p2(n) = det([x2(:,n) y2(:,n) EB.e2])>0;
end

return
